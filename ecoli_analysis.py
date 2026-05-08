"""
ecoli_analysis.py — Real-data validation on E. coli K-12 MG1655.

Question: at real homopolymer loci in a biological genome, do different
alignment parameter configs produce different CIGARs for the same read?

If yes → the mechanism demonstrated on synthetic data fires on real data.
The rate of config disagreement at homopolymer loci vs. non-homopolymer
loci is the key comparison.

Requirements
------------
    sudo apt install bwa samtools
    Files in same directory:
        ecoli_k12.fna          — E. coli K-12 MG1655 reference
        ecoli_reads.fastq.gz   — real Illumina reads (SRR2584863 R1)

Usage
-----
    python ecoli_analysis.py
"""

from __future__ import annotations

import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# 0.  Dependency check
# ---------------------------------------------------------------------------

for cmd, hint in [
    ("bwa",      "sudo apt install bwa"),
    ("samtools", "sudo apt install samtools"),
]:
    if subprocess.run(["which", cmd], capture_output=True).returncode != 0:
        sys.exit(f"\nERROR: '{cmd}' not found. Install with: {hint}\n")

REF    = Path("ecoli_k12.fna")
READS  = Path("ecoli_reads.fastq.gz")
OUTDIR = Path("ecoli_output")

for p in [REF, READS]:
    if not p.exists():
        sys.exit(f"\nERROR: {p} not found. See download instructions.\n")

OUTDIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1.  Parameter configs  (same six as the synthetic sweep)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Cfg:
    name: str
    A: int; B: int; O: int; E: int; L: int

CONFIGS = [
    Cfg("default",           A=2, B=4, O=6,  E=1,  L=100),
    Cfg("gap_friendly",      A=2, B=8, O=1,  E=1,  L=100),
    Cfg("mismatch_friendly", A=2, B=1, O=10, E=10, L=100),
    Cfg("balanced",          A=2, B=3, O=3,  E=3,  L=100),
    Cfg("local_default",     A=2, B=4, O=6,  E=1,  L=5),
    Cfg("local_gap_friendly",A=2, B=8, O=1,  E=1,  L=5),
]

# ---------------------------------------------------------------------------
# 2.  Find homopolymer loci in the reference
#     Returns list of (chrom, start, end, base) for runs >= min_len
# ---------------------------------------------------------------------------

def find_homopolymers(fna: Path, min_len: int = 6) -> list[tuple[str, int, int, str]]:
    print(f"  Scanning {fna.name} for homopolymer runs (>= {min_len}bp) ...")
    loci: list[tuple[str, int, int, str]] = []
    chrom = ""
    seq_parts: list[str] = []

    def flush(chrom, seq_parts):
        seq = "".join(seq_parts).upper()
        for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % (
                min_len, min_len, min_len, min_len), seq):
            loci.append((chrom, m.start(), m.end(), m.group()[0]))

    with open(fna) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if chrom:
                    flush(chrom, seq_parts)
                chrom = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if chrom:
            flush(chrom, seq_parts)

    print(f"  Found {len(loci)} homopolymer loci.")
    return loci

# ---------------------------------------------------------------------------
# 3.  Build BWA index (once)
# ---------------------------------------------------------------------------

def bwa_index(ref: Path) -> None:
    idx = ref.with_suffix(ref.suffix + ".bwt")
    if idx.exists():
        print(f"  Index exists, skipping.")
        return
    print(f"  Indexing {ref.name} ...", end=" ", flush=True)
    r = subprocess.run(["bwa", "index", str(ref)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"\nERROR: bwa index failed:\n{r.stderr}")
    print("done.")

# ---------------------------------------------------------------------------
# 4.  Align reads for each config, sort + index BAM
# ---------------------------------------------------------------------------

def align(ref: Path, reads: Path, cfg: Cfg, outdir: Path,
          threads: int = 4) -> Path:
    bam = outdir / f"{cfg.name}.bam"
    if bam.exists():
        print(f"  [{cfg.name}] BAM exists, skipping alignment.")
        return bam
    print(f"  Aligning [{cfg.name}] ...", end=" ", flush=True)
    bwa_cmd = [
        "bwa", "mem",
        f"-A{cfg.A}", f"-B{cfg.B}", f"-O{cfg.O}", f"-E{cfg.E}", f"-L{cfg.L}",
        "-t", str(threads),
        str(ref), str(reads),
    ]
    sort_cmd = ["samtools", "sort", "-o", str(bam), "-"]
    bwa  = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    sort = subprocess.Popen(sort_cmd, stdin=bwa.stdout,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    bwa.stdout.close()
    sort.wait()
    if sort.returncode != 0:
        sys.exit(f"\nERROR: alignment failed for {cfg.name}")
    subprocess.run(["samtools", "index", str(bam)],
                   capture_output=True, check=True)
    print("done.")
    return bam

# ---------------------------------------------------------------------------
# 5.  Extract CIGARs at homopolymer loci
#     For each locus, collect {read_name: cigar} for reads overlapping it.
#     We extend the window by ±5bp to capture reads that start just outside
#     the run but span it.
# ---------------------------------------------------------------------------

def cigars_at_locus(bam: Path, chrom: str, start: int, end: int,
                    window: int = 5) -> dict[str, str]:
    """Return {read_name: cigar} for primary alignments overlapping the locus."""
    region = f"{chrom}:{max(1, start - window)}-{end + window}"
    r = subprocess.run(
        ["samtools", "view", "-F", "0x900", str(bam), region],
        capture_output=True, text=True,
    )
    result: dict[str, str] = {}
    for line in r.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        result[fields[0]] = fields[5]   # read_name → CIGAR
    return result

# ---------------------------------------------------------------------------
# 6.  Classify CIGAR + MD  (same logic as bwa_validation.py)
# ---------------------------------------------------------------------------

def md_has_mismatches(md: str) -> bool:
    if not md:
        return False
    return bool(re.search(r"[ACGT]", re.sub(r"\^[ACGT]+", "", md)))

def classify(cigar: str, md: str = "") -> str:
    if not cigar or cigar == "*":
        return "UNMAPPED"
    ops = {op for _, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)}
    has_gap  = bool(ops & {"I", "D"})
    has_clip = bool(ops & {"S", "H"})
    has_m    = "M" in ops
    if not has_m and not has_gap:
        return "UNMAPPED"
    has_mm = md_has_mismatches(md)
    if has_clip and not has_gap and not has_mm:
        return "EXACT"
    if has_gap and has_mm:
        return "COMPLEX"
    if has_gap:
        return "INDEL"
    if has_mm:
        return "SUB"
    return "EXACT"

# ---------------------------------------------------------------------------
# 7.  Compute config-disagreement rate at homopolymer vs. other loci
#
#     For a sample of loci, collect all reads that appear under >= 2 configs.
#     A read "disagrees" if at least two configs give it different CIGARs.
#     Disagreement rate = (reads with any config disagreement) / (reads seen
#     under >= 2 configs).
# ---------------------------------------------------------------------------

def disagreement_rate(
    bams: dict[str, Path],
    loci: list[tuple[str, int, int, str]],
    sample: int = 500,
) -> tuple[float, int, int]:
    """
    Sample up to `sample` loci, count config disagreements.
    Returns (rate, n_disagreed, n_total_reads_compared).
    """
    import random
    random.seed(42)
    sampled = random.sample(loci, min(sample, len(loci)))

    disagreed = 0
    total_compared = 0

    for chrom, start, end, base in sampled:
        # Collect CIGARs from each config for reads at this locus
        per_cfg: dict[str, dict[str, str]] = {}
        for cfg_name, bam in bams.items():
            per_cfg[cfg_name] = cigars_at_locus(bam, chrom, start, end)

        # Find reads present in at least 2 configs
        all_reads: set[str] = set()
        for cigars in per_cfg.values():
            all_reads.update(cigars.keys())

        for read in all_reads:
            seen = [per_cfg[c][read] for c in per_cfg if read in per_cfg[c]]
            if len(seen) < 2:
                continue
            total_compared += 1
            if len(set(seen)) > 1:
                disagreed += 1

    rate = disagreed / total_compared if total_compared > 0 else 0.0
    return rate, disagreed, total_compared

# ---------------------------------------------------------------------------
# 8.  Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("\n=== E. coli K-12 Real-Data Validation ===\n")

    # Find homopolymer loci
    homo_loci = find_homopolymers(REF, min_len=6)

    # Index reference
    print("\nIndexing reference ...")
    bwa_index(REF)

    # Align all configs
    print("\nAligning reads (6 configs) ...")
    bams: dict[str, Path] = {}
    for cfg in CONFIGS:
        bams[cfg.name] = align(REF, READS, cfg, OUTDIR)

    # Compute disagreement at homopolymer loci
    print(f"\nComputing config disagreement at homopolymer loci")
    print(f"  (sampling 500 loci from {len(homo_loci)} total) ...")
    homo_rate, homo_dis, homo_total = disagreement_rate(bams, homo_loci, sample=500)

    # For comparison: sample the same number of non-homopolymer loci
    # We approximate by picking random positions in the genome
    print(f"\nComputing config disagreement at random (non-targeted) loci ...")
    # Build random loci from the reference chromosome
    import random
    random.seed(99)
    ref_len = 4_641_652  # MG1655 genome length
    chrom_name = None
    with open(REF) as fh:
        for line in fh:
            if line.startswith(">"):
                chrom_name = line[1:].split()[0]
                break
    rand_loci = [
        (chrom_name, pos, pos + 8, "N")
        for pos in random.sample(range(100, ref_len - 100), 500)
    ]
    rand_rate, rand_dis, rand_total = disagreement_rate(bams, rand_loci, sample=500)

    # Report
    print("\n" + "=" * 66)
    print("REAL-DATA CONFIG DISAGREEMENT REPORT")
    print("E. coli K-12 MG1655  ·  SRR2584863 Illumina reads")
    print("=" * 66)
    print(f"\n  Homopolymer loci (>= 6bp runs):")
    print(f"    Total loci in genome:      {len(homo_loci):>6,}")
    print(f"    Loci sampled:              {min(500, len(homo_loci)):>6,}")
    print(f"    Reads compared:            {homo_total:>6,}")
    print(f"    Reads with config conflict:{homo_dis:>6,}  ({homo_rate:.1%})")
    print(f"\n  Random (non-targeted) loci:")
    print(f"    Loci sampled:              {min(500, len(rand_loci)):>6,}")
    print(f"    Reads compared:            {rand_total:>6,}")
    print(f"    Reads with config conflict:{rand_dis:>6,}  ({rand_rate:.1%})")
    print(f"\n  Enrichment at homopolymer loci vs. random:")
    if rand_rate > 0:
        enrichment = homo_rate / rand_rate
        print(f"    {enrichment:.1f}x higher disagreement rate at homopolymer loci")
    else:
        print(f"    (no disagreement at random loci — cannot compute enrichment)")
    print("\n" + "=" * 66)

    # Save report
    report = OUTDIR / "ecoli_report.txt"
    with open(report, "w") as fh:
        fh.write(f"E. coli K-12 config disagreement\n")
        fh.write(f"Homopolymer loci: {homo_dis}/{homo_total} reads disagree ({homo_rate:.1%})\n")
        fh.write(f"Random loci:      {rand_dis}/{rand_total} reads disagree ({rand_rate:.1%})\n")
        if rand_rate > 0:
            fh.write(f"Enrichment: {homo_rate/rand_rate:.1f}x\n")
    print(f"\nReport saved to: {report}\n")


if __name__ == "__main__":
    main()
