"""
vcf_comparison.py — End-to-end downstream harm demonstration.

Question: does the parameter-induced CIGAR instability propagate to variant
calls? If gap_friendly and default produce different VCFs from the same reads,
the format-level information loss has measurable clinical consequence.

Pipeline
--------
1. Verify BAMs from ecoli_analysis.py exist and are indexed
2. Install/check bcftools
3. For each config: bcftools mpileup | bcftools call → .vcf
4. Compare VCFs: count variants unique to each config, shared across all,
   and classify differences by variant type (SNP vs INDEL)
5. Focus comparison: default vs gap_friendly (the configs that produce the
   most dramatic CIGAR flip in BWA-MEM — 77% of point-mutation reads)

Key question for the report
----------------------------
Are there variants called under gap_friendly that are NOT called under default,
or vice versa? Are there loci where gap_friendly calls an INDEL and default
calls a SNP for the same position?

Requirements
------------
    sudo apt install bcftools samtools
    Existing ecoli_output/*.bam from ecoli_analysis.py

Usage
-----
    python vcf_comparison.py

Outputs (in vcf_output/)
------------------------
    {config}.vcf.gz         — compressed, indexed VCF per config
    comparison_report.txt   — variant overlap table + key differences
    default_only.vcf        — variants called only under default config
    gap_friendly_only.vcf   — variants called only under gap_friendly config
    type_flip.tsv           — positions where default=SNP, gap_friendly=INDEL
                              (the most direct demonstration of downstream harm)
"""

from __future__ import annotations

import subprocess
import sys
import os
from pathlib import Path

# ─── Dependency check ────────────────────────────────────────────────────────

def require(cmd: str, hint: str) -> None:
    if subprocess.run(["which", cmd], capture_output=True).returncode != 0:
        sys.exit(f"\nERROR: '{cmd}' not found.\nInstall with: {hint}\n")

require("bcftools", "sudo apt install bcftools")
require("samtools", "sudo apt install samtools")

# ─── Paths ────────────────────────────────────────────────────────────────────

BAM_DIR  = Path("ecoli_output")
REF      = Path("ecoli_k12.fna")
OUT_DIR  = Path("vcf_output")
OUT_DIR.mkdir(exist_ok=True)

CONFIGS = [
    "default",
    "gap_friendly",
    "mismatch_friendly",
    "balanced",
    "local_default",
    "local_gap_friendly",
]

for p in [BAM_DIR, REF]:
    if not p.exists():
        sys.exit(f"\nERROR: {p} not found. Run ecoli_analysis.py first.\n")

# ─── Step 1: Verify and index BAMs ────────────────────────────────────────────

print("\nVerifying BAMs ...")
for cfg in CONFIGS:
    bam = BAM_DIR / f"{cfg}.bam"
    bai = BAM_DIR / f"{cfg}.bam.bai"
    if not bam.exists():
        sys.exit(f"\nERROR: {bam} not found. Run ecoli_analysis.py first.\n")
    if not bai.exists():
        print(f"  Indexing {bam.name} ...", end=" ", flush=True)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        print("done.")
    else:
        print(f"  {bam.name} — index exists.")

# ─── Step 2: Call variants per config ─────────────────────────────────────────

print("\nCalling variants (bcftools mpileup | call) ...")

def call_variants(cfg: str, ref: Path, bam_dir: Path, out_dir: Path) -> Path:
    bam     = bam_dir / f"{cfg}.bam"
    vcf_gz  = out_dir / f"{cfg}.vcf.gz"

    if vcf_gz.exists():
        print(f"  [{cfg}] VCF exists, skipping.")
        return vcf_gz

    print(f"  [{cfg}] calling ...", end=" ", flush=True)

    # mpileup → call pipeline
    mpileup_cmd = [
        "bcftools", "mpileup",
        "-f", str(ref),
        "-q", "20",          # min mapping quality
        "-Q", "20",          # min base quality
        "--annotate", "FORMAT/AD,FORMAT/DP",
        str(bam),
    ]
    call_cmd = [
        "bcftools", "call",
        "-mv",               # multiallelic caller, output only variant sites
        "--ploidy", "1",     # E. coli is haploid
        "-Oz",               # compressed VCF output
        "-o", str(vcf_gz),
    ]

    mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)
    call    = subprocess.Popen(call_cmd, stdin=mpileup.stdout,
                                stderr=subprocess.DEVNULL)
    mpileup.stdout.close()
    call.wait()

    if call.returncode != 0:
        sys.exit(f"\nERROR: bcftools call failed for {cfg}")

    # Index the VCF
    subprocess.run(["bcftools", "index", str(vcf_gz)],
                   check=True, capture_output=True)
    print("done.")
    return vcf_gz

vcf_files: dict[str, Path] = {}
for cfg in CONFIGS:
    vcf_files[cfg] = call_variants(cfg, REF, BAM_DIR, OUT_DIR)

# ─── Step 3: Count variants per config ────────────────────────────────────────

def count_variants(vcf_gz: Path) -> dict[str, int]:
    """Return {total, SNP, INDEL} counts for a VCF."""
    r = subprocess.run(
        ["bcftools", "stats", str(vcf_gz)],
        capture_output=True, text=True, check=True,
    )
    counts = {"total": 0, "SNP": 0, "INDEL": 0}
    for line in r.stdout.splitlines():
        if line.startswith("SN") and "number of SNPs:" in line:
            counts["SNP"] = int(line.split()[-1])
        elif line.startswith("SN") and "number of indels:" in line:
            counts["INDEL"] = int(line.split()[-1])
    counts["total"] = counts["SNP"] + counts["INDEL"]
    return counts

print("\nCounting variants per config ...")
var_counts: dict[str, dict] = {}
for cfg in CONFIGS:
    var_counts[cfg] = count_variants(vcf_files[cfg])
    c = var_counts[cfg]
    print(f"  {cfg:<22s}  total={c['total']:>5}  SNP={c['SNP']:>5}  INDEL={c['INDEL']:>5}")

# ─── Step 4: Pairwise comparison — default vs gap_friendly ────────────────────

print("\nComparing default vs gap_friendly ...")

def vcf_isec(vcf_a: Path, vcf_b: Path, name_a: str, name_b: str,
             out_dir: Path) -> dict:
    """
    Use bcftools isec to find variants private to each VCF and shared.
    Returns counts and writes private VCFs.
    """
    isec_dir = out_dir / f"isec_{name_a}_vs_{name_b}"
    isec_dir.mkdir(exist_ok=True)

    r = subprocess.run(
        ["bcftools", "isec", "-p", str(isec_dir),
         str(vcf_a), str(vcf_b)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(f"  WARNING: isec failed: {r.stderr[:200]}")
        return {}

    # isec output: 0000.vcf = private to A, 0001.vcf = private to B,
    #              0002.vcf = in both
    def count_lines(f: Path) -> int:
        if not f.exists():
            return 0
        return sum(1 for l in open(f) if not l.startswith("#") and l.strip())

    private_a = count_lines(isec_dir / "0000.vcf")
    private_b = count_lines(isec_dir / "0001.vcf")
    shared    = count_lines(isec_dir / "0002.vcf")

    return {
        f"private_{name_a}": private_a,
        f"private_{name_b}": private_b,
        "shared":             shared,
    }

isec = vcf_isec(
    vcf_files["default"],
    vcf_files["gap_friendly"],
    "default", "gap_friendly",
    OUT_DIR,
)

# ─── Step 5: Find type-flip positions ─────────────────────────────────────────
# Positions where default calls SNP and gap_friendly calls INDEL (or vice versa)
# These are the direct consequence of the CIGAR flip on variant calling.

print("\nFinding SNP↔INDEL type flips ...")

def get_variant_types(vcf_gz: Path) -> dict[str, str]:
    """Return {chrom:pos: 'SNP'|'INDEL'} for all variants in VCF."""
    r = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_gz)],
        capture_output=True, text=True, check=True,
    )
    result: dict[str, str] = {}
    for line in r.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, _, ref, alt = fields[:5]
        key = f"{chrom}:{pos}"
        # Simple type classification: INDEL if any allele differs in length
        is_indel = any(len(a) != len(ref) for a in alt.split(",") if a != ".")
        result[key] = "INDEL" if is_indel else "SNP"
    return result

types_default     = get_variant_types(vcf_files["default"])
types_gap_friendly = get_variant_types(vcf_files["gap_friendly"])

# Sites present in both but with different types
type_flips: list[tuple[str, str, str]] = []
shared_sites = set(types_default.keys()) & set(types_gap_friendly.keys())
for site in sorted(shared_sites):
    td = types_default[site]
    tg = types_gap_friendly[site]
    if td != tg:
        type_flips.append((site, td, tg))

# Write type-flip TSV
flip_tsv = OUT_DIR / "type_flip.tsv"
with open(flip_tsv, "w") as fh:
    fh.write("position\tdefault_type\tgap_friendly_type\n")
    for site, td, tg in type_flips:
        fh.write(f"{site}\t{td}\t{tg}\n")

# ─── Step 6: Report ───────────────────────────────────────────────────────────

report_lines = []
report_lines.append("\n" + "=" * 70)
report_lines.append("VCF COMPARISON REPORT")
report_lines.append("Does CIGAR instability propagate to variant calls?")
report_lines.append("=" * 70)

report_lines.append("\n--- Variants called per configuration ---")
report_lines.append(f"  {'config':<22s}  {'total':>7}  {'SNP':>7}  {'INDEL':>7}")
report_lines.append("  " + "-" * 44)
for cfg in CONFIGS:
    c = var_counts[cfg]
    report_lines.append(
        f"  {cfg:<22s}  {c['total']:>7}  {c['SNP']:>7}  {c['INDEL']:>7}"
    )

if isec:
    report_lines.append("\n--- default vs gap_friendly overlap ---")
    report_lines.append(f"  Private to default:      {isec.get('private_default', 'n/a'):>6}")
    report_lines.append(f"  Private to gap_friendly: {isec.get('private_gap_friendly', 'n/a'):>6}")
    report_lines.append(f"  Shared:                  {isec.get('shared', 'n/a'):>6}")

report_lines.append(f"\n--- SNP↔INDEL type flips at shared sites ---")
report_lines.append(f"  Positions with type disagreement: {len(type_flips)}")
if type_flips:
    report_lines.append(f"  First 10 examples:")
    for site, td, tg in type_flips[:10]:
        report_lines.append(f"    {site:<30s}  default={td:<5s}  gap_friendly={tg}")

report_lines.append("\n--- Key finding ---")
total_default     = var_counts["default"]["total"]
total_gap         = var_counts["gap_friendly"]["total"]
snp_default       = var_counts["default"]["SNP"]
snp_gap           = var_counts["gap_friendly"]["SNP"]
indel_default     = var_counts["default"]["INDEL"]
indel_gap         = var_counts["gap_friendly"]["INDEL"]

report_lines.append(
    f"  INDEL calls: default={indel_default}, gap_friendly={indel_gap} "
    f"({'+'if indel_gap > indel_default else ''}{indel_gap - indel_default})"
)
report_lines.append(
    f"  SNP calls:   default={snp_default}, gap_friendly={snp_gap} "
    f"({'+'if snp_gap > snp_default else ''}{snp_gap - snp_default})"
)
report_lines.append(
    f"  Type flips at shared sites: {len(type_flips)}"
)
if len(type_flips) > 0:
    report_lines.append(
        f"\n  YES — CIGAR instability propagates to variant calls.\n"
        f"  {len(type_flips)} position(s) where the same read data produces\n"
        f"  a SNP under default parameters and an INDEL under gap-friendly\n"
        f"  parameters (or vice versa). These represent direct clinical\n"
        f"  consequences of the format-level information loss."
    )
elif isec.get("private_default", 0) + isec.get("private_gap_friendly", 0) > 0:
    total_private = isec.get("private_default", 0) + isec.get("private_gap_friendly", 0)
    report_lines.append(
        f"\n  YES — {total_private} variant(s) called under one configuration\n"
        f"  but not the other. The CIGAR instability causes some variants\n"
        f"  to appear or disappear depending on scoring parameters."
    )
else:
    report_lines.append(
        f"\n  No differences detected at the VCF level.\n"
        f"  The variant caller may be robust to the CIGAR flip in this\n"
        f"  specific dataset, or coverage is insufficient to call variants\n"
        f"  at the affected sites. See type_flip.tsv."
    )

report_lines.append("\n" + "=" * 70 + "\n")

report = "\n".join(report_lines)
print(report)

report_path = OUT_DIR / "comparison_report.txt"
report_path.write_text(report)
print(f"Report saved to: {report_path}")
print(f"Type-flip TSV:   {flip_tsv}")

if __name__ == "__main__":
    pass
