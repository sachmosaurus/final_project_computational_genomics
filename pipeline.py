"""
pipeline.py — End-to-end demonstration of information loss in sequence alignment.

Covers all five modules from the project outline:
  1. Synthetic sequence generator (ambiguous cases with ground truth)
  2. Alignment runner (parameter sweep)
  3. Alignment classifier (CIGAR + alignment -> category)
  4. Flip detector (interpretation changes under parameter sweep)
  5. Output audit (what standard SAM output does and does not preserve)

Requires: pip install biopython
Usage:    python pipeline.py
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from itertools import islice
from typing import Iterator

from Bio import Align


# ============================================================
# MODULE 1 — SYNTHETIC SEQUENCE GENERATOR
# ============================================================
# Two classes of ambiguous cases, each with known ground truth:
#
#   A) HOMOPOLYMER ambiguity — a deletion from a run of identical bases.
#      The deletion is well-defined biologically, but the alignment has
#      multiple equally-optimal placements. (Demo 1 style.)
#
#   B) SUB-VS-INDEL ambiguity — a single mismatch that a gap-friendly
#      scoring scheme will re-interpret as a deletion+insertion pair.
#      (Demo 2 style.)

BASES = ("A", "C", "G", "T")


@dataclass
class Case:
    case_id: str
    ref: str
    read: str
    ground_truth: str           # "HOMOPOLYMER_DEL", "POINT_MUTATION"
    category: str               # "homopolymer_ambiguity", "sub_vs_indel"


def _rand_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def _mutate_one(seq: str, rng: random.Random) -> tuple[str, int]:
    """Change exactly one base in `seq` to a different base. Returns (new_seq, position)."""
    pos = rng.randrange(len(seq))
    others = [b for b in BASES if b != seq[pos]]
    return seq[:pos] + rng.choice(others) + seq[pos + 1:], pos


def generate_homopolymer_cases(n: int, rng: random.Random) -> list[Case]:
    """Read = reference with one base deleted from a homopolymer run.

    Multiple alignment placements of the gap will be equally optimal.
    """
    cases = []
    for i in range(n):
        flank_left = _rand_seq(rng.randint(45, 60), rng)
        homo_base = rng.choice(BASES)
        homo_len = rng.randint(6, 10)
        flank_right = _rand_seq(rng.randint(45, 60), rng)
        ref = flank_left + homo_base * homo_len + flank_right
        # Delete exactly one base from the homopolymer run
        read = flank_left + homo_base * (homo_len - 1) + flank_right
        cases.append(Case(
            case_id=f"homo_{i:03d}",
            ref=ref, read=read,
            ground_truth="HOMOPOLYMER_DEL",
            category="homopolymer_ambiguity",
        ))
    return cases


def generate_sub_vs_indel_cases(n: int, rng: random.Random) -> list[Case]:
    """Read = reference with exactly one mismatched base (same length).

    Gap-friendly parameters will re-interpret this as insertion+deletion.
    """
    cases = []
    for i in range(n):
        ref = _rand_seq(rng.randint(100, 130), rng)
        read, _ = _mutate_one(ref, rng)
        cases.append(Case(
            case_id=f"sub_{i:03d}",
            ref=ref, read=read,
            ground_truth="POINT_MUTATION",
            category="sub_vs_indel",
        ))
    return cases


def generate_corpus(n_each: int = 500, seed: int = 42) -> list[Case]:
    rng = random.Random(seed)
    return generate_homopolymer_cases(n_each, rng) + generate_sub_vs_indel_cases(n_each, rng)


# ============================================================
# MODULE 2 — ALIGNMENT RUNNER (parameter sweep)
# ============================================================

@dataclass(frozen=True)
class Config:
    name: str
    match: float
    mismatch: float
    open_gap: float
    extend_gap: float
    mode: str = "global"   # "global" or "local"


SWEEP: list[Config] = [
    # --- Global alignment (forces ends to align) ---
    Config("default",            match=2, mismatch=-4, open_gap=-6,  extend_gap=-1, mode="global"),
    Config("gap_friendly",       match=2, mismatch=-8, open_gap=-1,  extend_gap=-1, mode="global"),
    Config("mismatch_friendly",  match=2, mismatch=-1, open_gap=-10, extend_gap=-10, mode="global"),
    Config("balanced",           match=2, mismatch=-3, open_gap=-3,  extend_gap=-3, mode="global"),
    # --- Local alignment (can trim mismatched ends instead of forcing gaps) ---
    Config("local_default",      match=2, mismatch=-4, open_gap=-6,  extend_gap=-1, mode="local"),
    Config("local_gap_friendly", match=2, mismatch=-8, open_gap=-1,  extend_gap=-1, mode="local"),
]


def build_aligner(cfg: Config) -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = cfg.mode
    a.match_score = cfg.match
    a.mismatch_score = cfg.mismatch
    a.open_gap_score = cfg.open_gap
    a.extend_gap_score = cfg.extend_gap
    return a


def take_optimal(alignments, cap: int = 50) -> list:
    """Safely collect up to `cap` optimal alignments (the iterator can be huge)."""
    return list(islice(alignments, cap))


# ============================================================
# MODULE 3 — CIGAR EXTRACTION + CLASSIFIER
# ============================================================

def alignment_to_cigar(alignment) -> str:
    """CIGAR string using M for match/mismatch, I for insertion, D for deletion."""
    target, query = alignment[0], alignment[1]
    ops = []
    for t, q in zip(target, query):
        if t == "-":
            ops.append("I")
        elif q == "-":
            ops.append("D")
        else:
            ops.append("M")
    out, count = [], 1
    for i in range(1, len(ops)):
        if ops[i] == ops[i - 1]:
            count += 1
        else:
            out.append(f"{count}{ops[i - 1]}")
            count = 1
    out.append(f"{count}{ops[-1]}")
    return "".join(out)


# Maps each ground-truth category to the alignment category that correctly
# describes the underlying biological event. A classifier output that does NOT
# match this is a categorically wrong biological interpretation.
CORRECT_CATEGORY: dict[str, str] = {
    "POINT_MUTATION":   "SUB",
    "HOMOPOLYMER_DEL":  "INDEL",
}


def classify(alignment) -> str:
    """Categorize an alignment into {EXACT, SUB, INDEL, COMPLEX}.

    EXACT   -> no gaps, no mismatches
    SUB     -> mismatches only, no gaps
    INDEL   -> one or more gaps, no mismatches
    COMPLEX -> both gaps and mismatches
    """
    target, query = alignment[0], alignment[1]
    has_gap = False
    has_mismatch = False
    for t, q in zip(target, query):
        if t == "-" or q == "-":
            has_gap = True
        elif t != q:
            has_mismatch = True
    if not has_gap and not has_mismatch:
        return "EXACT"
    if has_gap and not has_mismatch:
        return "INDEL"
    if has_mismatch and not has_gap:
        return "SUB"
    return "COMPLEX"


# ============================================================
# MODULE 4 — FLIP DETECTOR
# ============================================================

@dataclass
class SweepResult:
    case: Case
    by_config: dict[str, dict] = field(default_factory=dict)
    # by_config[cfg_name] -> {"cigar": str, "category": str, "score": float,
    #                         "n_optimal": int}

    @property
    def categories(self) -> set[str]:
        return {v["category"] for v in self.by_config.values()}

    @property
    def cigars(self) -> set[str]:
        return {v["cigar"] for v in self.by_config.values()}

    @property
    def is_category_flip(self) -> bool:
        return len(self.categories) > 1

    @property
    def is_cigar_flip(self) -> bool:
        return len(self.cigars) > 1

    @property
    def has_inherent_ambiguity(self) -> bool:
        return any(v["n_optimal"] > 1 for v in self.by_config.values())

    def config_is_correct(self, cfg_name: str) -> bool:
        """True if this config's chosen CIGAR correctly describes the underlying event."""
        expected = CORRECT_CATEGORY[self.case.ground_truth]
        return self.by_config[cfg_name]["category"] == expected


def run_sweep(cases: list[Case], sweep: list[Config] = SWEEP) -> list[SweepResult]:
    results = []
    for case in cases:
        sr = SweepResult(case=case)
        for cfg in sweep:
            aligner = build_aligner(cfg)
            raw = aligner.align(case.ref, case.read)
            best = raw[0]
            # Count optimal alignments but cap at 100 — pathological homopolymer
            # cases can produce combinatorial blowups that hang the iterator.
            try:
                n_optimal = min(len(raw), 100)
            except OverflowError:
                n_optimal = 100
            sr.by_config[cfg.name] = {
                "cigar": alignment_to_cigar(best),
                "category": classify(best),
                "score": float(best.score),
                "n_optimal": n_optimal,
            }
        results.append(sr)
    return results


# ============================================================
# MODULE 5 — OUTPUT AUDIT (what SAM records, what it doesn't)
# ============================================================

def sam_line(read_id: str, ref_name: str, pos: int, cigar: str, seq: str) -> str:
    """Minimal SAM record. Standard SAM has ONE CIGAR field per primary alignment."""
    flag = 0
    mapq = 60
    rnext, pnext, tlen = "*", 0, 0
    qual = "*"
    return "\t".join([read_id, str(flag), ref_name, str(pos + 1),
                      str(mapq), cigar, rnext, str(pnext), str(tlen), seq, qual])


def write_sam(results: list[SweepResult], cfg_name: str, path: str) -> None:
    """Write a SAM file using one config's results. Deliberately mimics a real pipeline."""
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        fh.write("@SQ\tSN:synthetic_ref\tLN:1000\n")
        fh.write(f"@PG\tID:pipeline.py\tCN:midpoint_demo\tVN:0.1\tCL:config={cfg_name}\n")
        for sr in results:
            rec = sr.by_config[cfg_name]
            fh.write(sam_line(sr.case.case_id, "synthetic_ref", 0, rec["cigar"], sr.case.read) + "\n")


# ============================================================
# REPORTING
# ============================================================

def summarize(results: list[SweepResult]) -> None:
    n = len(results)
    n_cig_flip = sum(1 for r in results if r.is_cigar_flip)
    n_cat_flip = sum(1 for r in results if r.is_category_flip)
    n_inherent = sum(1 for r in results if r.has_inherent_ambiguity)

    print()
    print("=" * 66)
    print("OVERALL")
    print("=" * 66)
    print(f"  total cases:                               {n}")
    print(f"  cases with CIGAR flip across configs:      {n_cig_flip}  ({n_cig_flip/n:.0%})")
    print(f"  cases with CATEGORY flip (SUB<->INDEL):    {n_cat_flip}  ({n_cat_flip/n:.0%})")
    print(f"  cases with >1 equally-optimal alignment")
    print(f"  under at least one config (inherent amb.): {n_inherent}  ({n_inherent/n:.0%})")

    # Per-category breakdown
    print()
    print("=" * 66)
    print("BY GROUND-TRUTH CATEGORY")
    print("=" * 66)
    by_cat: dict[str, list[SweepResult]] = {}
    for r in results:
        by_cat.setdefault(r.case.category, []).append(r)
    for cat, group in by_cat.items():
        ncf = sum(1 for r in group if r.is_cigar_flip)
        ncat = sum(1 for r in group if r.is_category_flip)
        ninh = sum(1 for r in group if r.has_inherent_ambiguity)
        print(f"  [{cat}]  n={len(group)}")
        print(f"     CIGAR flips: {ncf}/{len(group)}   category flips: {ncat}/{len(group)}   "
              f"inherent amb.: {ninh}/{len(group)}")

    # A worked example of a category flip
    print()
    print("=" * 66)
    print("WORKED EXAMPLE — CATEGORY FLIP")
    print("=" * 66)
    examples = [r for r in results if r.is_category_flip]
    if examples:
        ex = examples[0]
        print(f"  case: {ex.case.case_id}  (ground truth: {ex.case.ground_truth})")
        print(f"  ref:  {ex.case.ref}")
        print(f"  read: {ex.case.read}")
        print()
        for cfg_name, rec in ex.by_config.items():
            print(f"    under '{cfg_name}':  category={rec['category']:7s}"
                  f"  CIGAR={rec['cigar']:20s}  score={rec['score']}")
        print()
        print(f"  -> Same input. Different biological claim per config.")
        print(f"     The rejected interpretations are not recorded anywhere in SAM.")

    # A worked example of inherent ambiguity
    print()
    print("=" * 66)
    print("WORKED EXAMPLE — INHERENT AMBIGUITY")
    print("=" * 66)
    inherent = [r for r in results if r.has_inherent_ambiguity]
    if inherent:
        ex = inherent[0]
        # Pick the config with the most optimal alignments
        cfg_name = max(ex.by_config, key=lambda c: ex.by_config[c]["n_optimal"])
        rec = ex.by_config[cfg_name]
        print(f"  case: {ex.case.case_id}  (ground truth: {ex.case.ground_truth})")
        print(f"  ref:  {ex.case.ref}")
        print(f"  read: {ex.case.read}")
        print(f"  under '{cfg_name}':  {rec['n_optimal']} equally-optimal alignments found")
        print(f"  -> Standard SAM output records ONE. The other "
              f"{rec['n_optimal']-1} leave no trace.")

    # Ground-truth recovery: which configs correctly recover the underlying event
    print()
    print("=" * 66)
    print("GROUND-TRUTH RECOVERY (is the reported category biologically correct?)")
    print("=" * 66)
    configs = list(results[0].by_config.keys())
    ground_truths = sorted({r.case.ground_truth for r in results})

    # Header row
    header = f"  {'config':<22s}" + "".join(f"{gt:>22s}" for gt in ground_truths) + f"{'overall':>14s}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for cfg in configs:
        row = f"  {cfg:<22s}"
        correct_total = 0
        total = 0
        for gt in ground_truths:
            subset = [r for r in results if r.case.ground_truth == gt]
            correct = sum(1 for r in subset if r.config_is_correct(cfg))
            correct_total += correct
            total += len(subset)
            row += f"{correct}/{len(subset)} ({correct/len(subset):.0%})".rjust(22)
        row += f"{correct_total}/{total} ({correct_total/total:.0%})".rjust(14)
        print(row)
    print()
    print("  A row with less than 100% in a ground-truth column means: that config's")
    print("  scoring model produces biologically incorrect interpretations at that rate.")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    corpus = generate_corpus(n_each=500, seed=42)
    print(f"Generated corpus: {len(corpus)} cases "
          f"({sum(1 for c in corpus if c.category == 'homopolymer_ambiguity')} homopolymer + "
          f"{sum(1 for c in corpus if c.category == 'sub_vs_indel')} sub-vs-indel)")

    print(f"Running sweep over {len(SWEEP)} configs: "
          f"{[c.name for c in SWEEP]}")
    results = run_sweep(corpus)

    summarize(results)

    # Module 5: write a SAM file under one config. Demonstrates that the
    # on-disk record is single-CIGAR-per-read — the alternatives are gone.
    sam_path = "synthetic.default.sam"
    write_sam(results, cfg_name="default", path=sam_path)
    print()
    print("=" * 66)
    print(f"WROTE: {sam_path}  (config=default)")
    print("=" * 66)
    print("  This is what a standard downstream pipeline would consume.")
    print("  Note: exactly one CIGAR per read. No field stores the")
    print("  alternative interpretations from other configs.")
