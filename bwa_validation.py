"""
bwa_validation.py — Validate the category-flip finding with BWA-MEM.

Asks: does the parameter-induced category flip (point mutation misreported as
indel) replicate when we swap Biopython's classical DP aligner for BWA-MEM's
seed-chain-extend?

If yes  →  the failure is format-structural; it survives a production aligner.
If no   →  the claim must be qualified to DP aligners specifically.

Strategy
--------
1. Generate the same 1000-case synthetic corpus (seed=42, same as Phase 1).
2. Write one concatenated reference FASTA (1000 sequences, one per case).
3. Index it once with `bwa index`.
4. For each parameter config, run `bwa mem` with the equivalent -A/-B/-O/-E
   flags and write the SAM.
5. Parse each SAM: classify each alignment as EXACT/SUB/INDEL/COMPLEX using
   the same logic as pipeline.py.
6. Verify each read mapped to its own reference segment (sanity check).
7. Print a comparison table: Biopython results vs BWA-MEM results per config.

Requirements
------------
    sudo apt install bwa          # or: conda install bwa
    pip install biopython         # for corpus generation

Usage
-----
    python bwa_validation.py

Output files (in ./bwa_output/)
--------------------------------
    refs.fa              — concatenated reference FASTA
    reads.fq             — all reads as FASTQ (dummy Q40 quality)
    refs.fa.{amb,ann,…}  — BWA index files
    <config_name>.sam    — one SAM per parameter config
    validation_report.txt — summary comparison table
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# 0.  Dependency check — fail early with a useful message
# ---------------------------------------------------------------------------

def _require(cmd: str, install_hint: str) -> None:
    if subprocess.run(["which", cmd], capture_output=True).returncode != 0:
        sys.exit(
            f"\nERROR: '{cmd}' not found on PATH.\n"
            f"Install with:  {install_hint}\n"
            f"Then re-run this script.\n"
        )

_require("bwa",    "sudo apt install bwa   (or: conda install -c bioconda bwa)")

# ---------------------------------------------------------------------------
# 1.  Import corpus generation from pipeline.py
# ---------------------------------------------------------------------------

try:
    from pipeline import generate_corpus, CORRECT_CATEGORY
except ModuleNotFoundError:
    sys.exit(
        "\nERROR: pipeline.py not found in the current directory.\n"
        "Run this script from the same folder as pipeline.py.\n"
    )

# ---------------------------------------------------------------------------
# 2.  BWA-MEM parameter configs
#     Mirrors the SWEEP in pipeline.py as closely as BWA-MEM flags allow.
#
#     Key flags:
#       -A  match score            (pipeline: match)
#       -B  mismatch penalty       (pipeline: abs(mismatch))
#       -O  gap-open penalty       (pipeline: abs(open_gap))
#       -E  gap-extend penalty     (pipeline: abs(extend_gap))
#       -L  clipping penalty       high value → discourages soft-clipping
#                                  (approximates global alignment end-behavior)
#
#     BWA-MEM is always local internally (Smith-Waterman), so "global" here
#     means: penalise clipping heavily so the aligner is forced to use the
#     full read rather than trimming ends.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BwaConfig:
    name: str
    A: int   # match
    B: int   # mismatch penalty (positive int; bwa negates internally)
    O: int   # gap-open penalty
    E: int   # gap-extend penalty
    clip_penalty: int  # -L: high = suppress soft-clipping (≈ global)

BWA_CONFIGS: list[BwaConfig] = [
    # ---- Global-like (high clip penalty suppresses end-trimming) ----
    BwaConfig("default",           A=2, B=4, O=6,  E=1,  clip_penalty=100),
    BwaConfig("gap_friendly",      A=2, B=8, O=1,  E=1,  clip_penalty=100),
    BwaConfig("mismatch_friendly", A=2, B=1, O=10, E=10, clip_penalty=100),
    BwaConfig("balanced",          A=2, B=3, O=3,  E=3,  clip_penalty=100),
    # ---- Local-like (default clip penalty allows end-trimming) ------
    BwaConfig("local_default",     A=2, B=4, O=6,  E=1,  clip_penalty=5),
    BwaConfig("local_gap_friendly",A=2, B=8, O=1,  E=1,  clip_penalty=5),
]

# ---------------------------------------------------------------------------
# 3.  CIGAR parser + classifier
#
#     The fundamental problem: SAM's M operation covers both matches and
#     mismatches.  A homopolymer-deletion CIGAR like "48M1D61M" has both M
#     and D, but the M regions are pure matches — it should be INDEL, not
#     COMPLEX.  COMPLEX requires actual mismatches co-occurring with gaps.
#
#     Fix: use the MD auxiliary tag that BWA-MEM always emits.  MD encodes
#     exactly which reference positions under M ops were mismatches, without
#     needing to re-fetch the reference.
#
#     MD syntax recap:
#       - Digits  → that many consecutive matches (no mismatch)
#       - Letter  → one mismatch (the reference base at that position)
#       - ^Letters → deleted reference bases (inside a D op)
#     So MD="109"    → 109 consecutive matches, no mismatches
#        MD="83A42"  → 83 matches, then ref=A (mismatch), then 42 matches
#        MD="48^A61" → 48 matches, 1-base deletion, 61 matches
# ---------------------------------------------------------------------------

def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """Return list of (length, op) tuples from a CIGAR string, e.g. '83M1D1I42M'."""
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]

def md_has_mismatches(md: str) -> bool:
    """Return True if the MD tag contains any substitution mismatches.

    Deletions in MD look like '^ACG' — we skip those.
    Any remaining letter (not a digit, not '^', not part of a deletion block)
    is a substitution mismatch at that reference position.
    """
    if not md:
        return False
    # Remove deletion blocks: ^[ACGT]+
    md_no_del = re.sub(r"\^[ACGT]+", "", md)
    # Any remaining letter is a mismatch
    return bool(re.search(r"[ACGT]", md_no_del))

def cigar_category(cigar: str, md: str = "") -> str:
    """
    Classify a SAM alignment into {EXACT, SUB, INDEL, COMPLEX, UNMAPPED}.

    Uses CIGAR for gap detection and MD tag for mismatch detection.
    This resolves SAM's M ambiguity: M can be match or mismatch, but MD
    tells us which — so we never need to re-fetch the reference.

    EXACT   : only M, no gaps, no mismatches (MD is all digits)
    SUB     : only M, no gaps, has mismatches (MD contains letters)
    INDEL   : has I/D, no mismatches in M regions (MD has no letters outside ^)
    COMPLEX : has I/D and also has mismatches in M regions
    UNMAPPED: no alignment
    """
    if cigar in ("*", "") or cigar is None:
        return "UNMAPPED"

    ops = parse_cigar(cigar)
    op_types = {op for _, op in ops}

    has_gap  = bool(op_types & {"I", "D"})
    has_clip = bool(op_types & {"S", "H"})
    has_m    = "M" in op_types

    if not has_m and not has_gap:
        return "UNMAPPED"

    # Use MD to determine whether M regions contain mismatches
    has_mismatch = md_has_mismatches(md) if md else False

    if has_clip and not has_gap and not has_mismatch:
        # Local alignment trimmed mismatched ends — variant silently dropped
        return "EXACT"

    if has_gap and has_mismatch:
        return "COMPLEX"
    if has_gap:
        return "INDEL"
    if has_mismatch:
        return "SUB"
    return "EXACT"

# ---------------------------------------------------------------------------
# 4.  Write corpus to FASTA / FASTQ
# ---------------------------------------------------------------------------

DUMMY_QUAL_CHAR = "I"  # Phred 40 — BWA-MEM doesn't use base qualities for scoring

def write_fasta(cases, path: Path) -> None:
    with open(path, "w") as fh:
        for case in cases:
            fh.write(f">{case.case_id}\n{case.ref}\n")

def write_fastq(cases, path: Path) -> None:
    with open(path, "w") as fh:
        for case in cases:
            qual = DUMMY_QUAL_CHAR * len(case.read)
            fh.write(f"@{case.case_id}\n{case.read}\n+\n{qual}\n")

# ---------------------------------------------------------------------------
# 5.  Run BWA-MEM
# ---------------------------------------------------------------------------

def bwa_index(ref_fa: Path) -> None:
    print(f"  Indexing {ref_fa.name} ...", end=" ", flush=True)
    r = subprocess.run(
        ["bwa", "index", str(ref_fa)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        sys.exit(f"\nERROR: bwa index failed:\n{r.stderr}")
    print("done.")

def bwa_mem(ref_fa: Path, reads_fq: Path, cfg: BwaConfig, out_sam: Path,
            threads: int = 4) -> None:
    cmd = [
        "bwa", "mem",
        f"-A{cfg.A}",
        f"-B{cfg.B}",
        f"-O{cfg.O}",
        f"-E{cfg.E}",
        f"-L{cfg.clip_penalty}",
        "-t", str(threads),
        # Suppress secondary alignments to keep output clean
        "-c", "1",
        str(ref_fa),
        str(reads_fq),
    ]
    print(f"  Running bwa mem [{cfg.name}] ...", end=" ", flush=True)
    with open(out_sam, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        sys.exit(f"\nERROR: bwa mem failed for config '{cfg.name}':\n{r.stderr}")
    print("done.")

# ---------------------------------------------------------------------------
# 6.  Parse SAM → {read_id: {cigar, category, rname}}
# ---------------------------------------------------------------------------

def parse_sam(sam_path: Path) -> dict[str, dict]:
    """Return {read_id: {cigar, category, rname, flag}} for primary alignments.

    Extracts the MD auxiliary tag so the classifier can distinguish matches
    from mismatches inside M operations.
    """
    records: dict[str, dict] = {}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            read_id, flag_str, rname, _, _, cigar = fields[:6]
            flag = int(flag_str)
            if flag & 0x100 or flag & 0x800:
                continue
            if read_id in records:
                continue
            md = ""
            for opt in fields[11:]:
                if opt.startswith("MD:Z:"):
                    md = opt[5:]
                    break
            records[read_id] = {
                "cigar":    cigar,
                "category": cigar_category(cigar, md),
                "rname":    rname,
                "flag":     flag,
            }
    return records

# ---------------------------------------------------------------------------
# 7.  Sanity check: did each read map to its own reference?
# ---------------------------------------------------------------------------

def check_mapping_fidelity(cases, sam_records: dict, config_name: str) -> float:
    """Return fraction of reads that mapped to their own reference segment."""
    total = len(cases)
    correct_rname = sum(
        1 for case in cases
        if sam_records.get(case.case_id, {}).get("rname") == case.case_id
    )
    fidelity = correct_rname / total if total else 0.0
    if fidelity < 0.95:
        print(
            f"  WARNING [{config_name}]: only {correct_rname}/{total} reads mapped "
            f"to their own reference segment ({fidelity:.1%}). "
            f"Results for mis-mapped reads may be unreliable."
        )
    return fidelity

# ---------------------------------------------------------------------------
# 8.  Compute ground-truth recovery table
# ---------------------------------------------------------------------------

def recovery_table(cases, sam_records: dict) -> dict[str, dict[str, int]]:
    """
    Returns {ground_truth: {category: count}}.
    Only counts reads that mapped to their own reference (fidelity check).
    """
    table: dict[str, dict[str, int]] = {}
    for case in cases:
        rec = sam_records.get(case.case_id)
        if rec is None or rec["rname"] != case.case_id:
            continue
        gt = case.ground_truth
        cat = rec["category"]
        table.setdefault(gt, {})
        table[gt][cat] = table[gt].get(cat, 0) + 1
    return table

# ---------------------------------------------------------------------------
# 9.  Print comparison report
# ---------------------------------------------------------------------------

BIOPYTHON_RESULTS = {
    # From pipeline.py run with seed=42, n=500 per class.
    # Format: {config_name: {ground_truth: {category: count}}}
    "default": {
        "POINT_MUTATION":  {"SUB": 500},
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
    "gap_friendly": {
        "POINT_MUTATION":  {"INDEL": 500},         # 100% FLIP
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
    "mismatch_friendly": {
        "POINT_MUTATION":  {"SUB": 500},
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
    "balanced": {
        "POINT_MUTATION":  {"SUB": 500},
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
    "local_default": {
        "POINT_MUTATION":  {"SUB": 473, "EXACT": 27},  # 27 silently dropped
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
    "local_gap_friendly": {
        "POINT_MUTATION":  {"INDEL": 500},             # 100% FLIP
        "HOMOPOLYMER_DEL": {"INDEL": 500},
    },
}

GROUND_TRUTHS  = ["POINT_MUTATION", "HOMOPOLYMER_DEL"]
CORRECT_CAT = {"POINT_MUTATION": "SUB", "HOMOPOLYMER_DEL": "INDEL"}

def correct_count(table: dict, gt: str, n: int = 500) -> tuple[int, int]:
    """Return (correct, total) for a given ground truth."""
    cat_counts = table.get(gt, {})
    total = sum(cat_counts.values())
    correct = cat_counts.get(CORRECT_CAT[gt], 0)
    return correct, total

def print_report(all_results: dict[str, dict], output_path: Path) -> None:
    lines = []

    header = (
        "\n" + "=" * 72 + "\n"
        "BWA-MEM VALIDATION REPORT\n"
        "Does the category-flip finding replicate on a seed-chain-extend aligner?\n"
        + "=" * 72
    )
    lines.append(header)

    col_w = 22
    for gt in GROUND_TRUTHS:
        lines.append(f"\n--- Ground truth: {gt} ---")
        lines.append(
            f"  {'config':<22s}  {'BWA-MEM correct':>16s}  "
            f"{'Biopython correct':>17s}  {'BWA-MEM category breakdown'}"
        )
        lines.append("  " + "-" * 70)

        for cfg in BWA_CONFIGS:
            bwa_table  = all_results.get(cfg.name, {})
            bwa_c, bwa_t = correct_count(bwa_table, gt)
            bio_data   = BIOPYTHON_RESULTS.get(cfg.name, {}).get(gt, {})
            bio_c      = bio_data.get(CORRECT_CAT[gt], 0)
            bio_t      = sum(bio_data.values()) or 500
            bwa_cats   = bwa_table.get(gt, {})
            cat_str    = "  ".join(f"{k}={v}" for k, v in sorted(bwa_cats.items()))
            bwa_label  = f"{bwa_c}/{bwa_t} ({bwa_c/bwa_t:.0%})" if bwa_t else "n/a"
            bio_label  = f"{bio_c}/{bio_t} ({bio_c/bio_t:.0%})"
            lines.append(
                f"  {cfg.name:<22s}  {bwa_label:>16s}  {bio_label:>17s}  {cat_str}"
            )

    # Key question answer
    lines.append("\n" + "=" * 72)
    lines.append("KEY QUESTION: does gap_friendly flip replicate in BWA-MEM?")
    lines.append("=" * 72)
    gf_bwa = all_results.get("gap_friendly", {})
    pm_cats = gf_bwa.get("POINT_MUTATION", {})
    total_pm = sum(pm_cats.values())
    indel_pm = pm_cats.get("INDEL", 0) + pm_cats.get("COMPLEX", 0)
    if total_pm > 0:
        flip_rate = indel_pm / total_pm
        if flip_rate > 0.8:
            verdict = (
                f"  YES — {indel_pm}/{total_pm} ({flip_rate:.0%}) point-mutation reads\n"
                f"  reported as INDEL/COMPLEX under gap_friendly BWA-MEM.\n"
                f"  The category flip is NOT a Biopython artifact.\n"
                f"  It survives a production seed-chain-extend aligner.\n"
                f"  The failure is format-structural."
            )
        elif flip_rate > 0.2:
            verdict = (
                f"  PARTIAL — {indel_pm}/{total_pm} ({flip_rate:.0%}) point-mutation reads\n"
                f"  reported as INDEL/COMPLEX under gap_friendly BWA-MEM.\n"
                f"  The flip rate is lower than Biopython (100%) but non-zero.\n"
                f"  BWA-MEM's heuristics partially suppress but do not eliminate the flip.\n"
                f"  The format-structural claim still holds."
            )
        else:
            verdict = (
                f"  NO (or weak) — only {indel_pm}/{total_pm} ({flip_rate:.0%}) point-mutation\n"
                f"  reads reported as INDEL/COMPLEX under gap_friendly BWA-MEM.\n"
                f"  BWA-MEM's heuristics suppress the gap-friendly flip significantly.\n"
                f"  The claim should be qualified: the flip is parameter-induced in DP\n"
                f"  aligners and weaker (but format-structural concern remains) in BWA-MEM."
            )
        lines.append(verdict)
    else:
        lines.append("  Could not determine — no POINT_MUTATION records found for gap_friendly.")

    lines.append("\n" + "=" * 72 + "\n")

    report = "\n".join(lines)
    print(report)
    output_path.write_text(report)
    print(f"\nReport saved to: {output_path}")

# ---------------------------------------------------------------------------
# 10.  Main
# ---------------------------------------------------------------------------

def main() -> None:
    out_dir = Path("bwa_output")
    out_dir.mkdir(exist_ok=True)

    # Generate corpus
    print("\nGenerating synthetic corpus (seed=42, n=500 per class) ...")
    cases = generate_corpus(n_each=500, seed=42)
    print(f"  {len(cases)} cases generated.")

    # Write FASTA and FASTQ
    ref_fa   = out_dir / "refs.fa"
    reads_fq = out_dir / "reads.fq"
    write_fasta(cases, ref_fa)
    write_fastq(cases, reads_fq)
    print(f"  Reference FASTA: {ref_fa}")
    print(f"  Reads FASTQ:     {reads_fq}")

    # Index reference once
    print("\nIndexing reference ...")
    bwa_index(ref_fa)

    # Run each config
    print("\nRunning BWA-MEM parameter sweep ...")
    all_results: dict[str, dict] = {}

    for cfg in BWA_CONFIGS:
        sam_path = out_dir / f"{cfg.name}.sam"
        bwa_mem(ref_fa, reads_fq, cfg, sam_path)
        records = parse_sam(sam_path)
        check_mapping_fidelity(cases, records, cfg.name)
        all_results[cfg.name] = recovery_table(cases, records)

    # Report
    report_path = out_dir / "validation_report.txt"
    print_report(all_results, report_path)


if __name__ == "__main__":
    main()
