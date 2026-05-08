"""
flip_demo.py — Minimum viable demonstration for the midpoint.

Thesis: alignment algorithms output ONE CIGAR per read, even when multiple
interpretations are equally valid or parameter-dependent. The alternatives
leave no trace.

Two demonstrations:
  1. Inherent ambiguity: for a homopolymer read, multiple alignments
     are equally optimal. The aligner picks one; the others are lost.
  2. Parameter-induced flip: the same (ref, read) pair yields a
     different "best" alignment under different parameters.

Requires: pip install biopython
"""

from Bio import Align


def make_aligner(match, mismatch, open_gap, extend_gap):
    """Configure a global pairwise aligner with explicit scoring."""
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score = match
    a.mismatch_score = mismatch
    a.open_gap_score = open_gap
    a.extend_gap_score = extend_gap
    return a


def alignment_to_cigar(alignment):
    """Convert a Biopython PairwiseAlignment to a CIGAR-like string.

    Uses M for match/mismatch (as SAM does), I for insertion to reference,
    D for deletion from reference.
    """
    target, query = alignment[0], alignment[1]
    ops = []
    for t, q in zip(target, query):
        if t == "-":
            ops.append("I")
        elif q == "-":
            ops.append("D")
        else:
            ops.append("M")
    # Run-length encode
    cigar = []
    count = 1
    for i in range(1, len(ops)):
        if ops[i] == ops[i - 1]:
            count += 1
        else:
            cigar.append(f"{count}{ops[i - 1]}")
            count = 1
    cigar.append(f"{count}{ops[-1]}")
    return "".join(cigar)


def show(alignment, label=""):
    if label:
        print(f"[{label}]  score={alignment.score}  CIGAR={alignment_to_cigar(alignment)}")
    print(alignment)


# ============================================================
# DEMO 1 — Inherent ambiguity in a homopolymer
# ============================================================
# A 4bp read aligned against a 5bp homopolymer. There are 5 equally
# optimal alignments (the gap can go in any of 5 positions). Every
# standard aligner picks exactly one and silently drops the other four.

print("=" * 62)
print("DEMO 1 — Inherent ambiguity (5 equally-optimal alignments)")
print("=" * 62)

ref_1, read_1 = "AAAAA", "AAAA"
aligner_1 = make_aligner(match=1, mismatch=-1, open_gap=-1, extend_gap=-1)
alignments = aligner_1.align(ref_1, read_1)

print(f"\nref:  {ref_1}")
print(f"read: {read_1}")
print(f"Optimal alignments found: {len(alignments)}   (all with score {alignments[0].score})\n")

for i, aln in enumerate(alignments, 1):
    print(f"  Alternative {i}:  CIGAR = {alignment_to_cigar(aln)}")
    for line in str(aln).splitlines():
        print(f"    {line}")
    print()

print("-> A standard SAM output would record ONE of these and discard the rest.\n")


# ============================================================
# DEMO 2 — Parameter-induced flip (sub vs. indel)
# ============================================================
# Same reference, same read, both 12 bp. One position differs.
# Under mismatch-friendly parameters the aligner reports ONE point mutation.
# Under gap-friendly parameters the aligner reports ONE deletion + ONE
# insertion — a fundamentally different biological claim.

print("=" * 62)
print("DEMO 2 — Parameter-induced flip (point mutation vs. indel pair)")
print("=" * 62)

ref_2  = "ACGTACGTACGT"
read_2 = "ACGTAAGTACGT"   # identical length; one position altered

print(f"\nref:  {ref_2}")
print(f"read: {read_2}\n")

configs = {
    "MISMATCH-FRIENDLY (mismatch cheap, gaps expensive)": dict(
        match=2, mismatch=-1, open_gap=-10, extend_gap=-10
    ),
    "GAP-FRIENDLY (mismatch expensive, gaps cheap)": dict(
        match=2, mismatch=-20, open_gap=-1, extend_gap=-1
    ),
}

for label, cfg in configs.items():
    aligner = make_aligner(**cfg)
    best = aligner.align(ref_2, read_2)[0]
    print(f"--- {label} ---")
    print(f"    params: {cfg}")
    print(f"    result: score={best.score}  CIGAR={alignment_to_cigar(best)}")
    for line in str(best).splitlines():
        print(f"    {line}")
    print()

print("-> Interpretation flipped from '1 point mutation' to '1 deletion + 1 insertion'.")
print("   Same input. Different biology reported. The rejected claim is not recorded.")
