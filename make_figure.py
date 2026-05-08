"""
make_figure.py — Generate the midpoint-presentation figure.

Two-panel stacked bar chart: for each ground-truth class (homopolymer_del,
point_mutation), show the category distribution (EXACT/SUB/INDEL/COMPLEX)
produced by each of the four parameter configs.

The visual argument:
  - For point mutations, gap_friendly flips ~100% of reads from SUB to INDEL/COMPLEX
  - For homopolymer deletions, all configs agree on category (INDEL) — but the
    "inherent ambiguity" panel reminds us that multiple equally-optimal
    alignments exist at that category.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from pipeline import generate_corpus, run_sweep, SWEEP, CORRECT_CATEGORY

# Navy/teal palette matching Victor's outline deck.
NAVY    = "#0A1E3C"
TEAL    = "#14B8C6"
SAND    = "#F5E6C8"
CORAL   = "#E85A4F"
AMBER   = "#F59E0B"
MUTED   = "#CBD5E1"

CATEGORY_ORDER  = ["EXACT", "SUB", "INDEL", "COMPLEX"]
CATEGORY_COLORS = {
    "EXACT":   AMBER,   # for point-mutation inputs, EXACT = silently dropped variant
    "SUB":     TEAL,    # "point mutation" reading
    "INDEL":   NAVY,    # "deletion/insertion" reading
    "COMPLEX": CORAL,   # mixed — usually a mis-interpretation artifact
}

mpl.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Inter", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        10,
    "axes.titlesize":   12,
    "axes.titleweight": "bold",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})


def tally(results, ground_truth: str) -> dict[str, dict[str, int]]:
    """For each config, count how many cases of the given ground truth fell
    into each category. Returns {config_name: {category: count}}."""
    out: dict[str, dict[str, int]] = {}
    for cfg in SWEEP:
        counts = {c: 0 for c in CATEGORY_ORDER}
        for r in results:
            if r.case.ground_truth != ground_truth:
                continue
            cat = r.by_config[cfg.name]["category"]
            counts[cat] = counts.get(cat, 0) + 1
        out[cfg.name] = counts
    return out


def draw_panel(ax, counts_by_cfg: dict[str, dict[str, int]],
               title: str, correct_category: str, n_cases: int) -> None:
    configs = list(counts_by_cfg.keys())
    x = range(len(configs))
    bottoms = [0] * len(configs)

    for cat in CATEGORY_ORDER:
        heights = [counts_by_cfg[cfg][cat] for cfg in configs]
        ax.bar(
            x, heights, bottom=bottoms,
            color=CATEGORY_COLORS[cat], label=cat,
            edgecolor="white", linewidth=1.0,
        )
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    ax.set_xticks(x)
    ax.set_xticklabels(configs, rotation=30, fontsize=9, ha="right")
    ax.set_ylim(0, n_cases * 1.15)
    ax.set_ylabel(f"cases classified (of {n_cases})", fontsize=9)
    ax.set_title(title, loc="left", pad=12)

    # Annotate each bar with correctness rate against ground truth
    for i, cfg in enumerate(configs):
        correct = counts_by_cfg[cfg].get(correct_category, 0)
        pct = correct / n_cases
        color = "#0A1E3C" if pct == 1.0 else ("#E85A4F" if pct < 0.5 else "#B45309")
        label = f"{correct}/{n_cases}"
        pct_label = f"({pct:.0%})"
        ax.text(
            i, n_cases * 1.08,
            label, ha="center", fontsize=8,
            color=color, fontweight="bold",
        )
        ax.text(
            i, n_cases * 1.02,
            pct_label, ha="center", fontsize=8,
            color=color,
        )

    ax.axhline(n_cases, color="#94A3B8", linewidth=0.6, linestyle="--", zorder=0)


def draw_ambiguity_panel(ax, results, ground_truth: str, title: str) -> None:
    """For each config, show the mean number of equally-optimal alignments on the given
    ground-truth class. Highlights that even 'correct' configs pick
    arbitrarily from multiple equally-optimal alignments."""
    subset = [r for r in results if r.case.ground_truth == ground_truth]
    configs = [cfg.name for cfg in SWEEP]
    x = range(len(configs))

    means = []
    for cfg in configs:
        ns = [r.by_config[cfg]["n_optimal"] for r in subset]
        means.append(sum(ns) / len(ns))

    ax.bar(x, means, color=TEAL, edgecolor="white", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(configs, rotation=30, fontsize=9, ha="right")
    ax.set_ylim(0, max(means) * 1.3)
    ax.set_ylabel("mean equally-optimal alignments per read", fontsize=9)
    ax.set_title(title, loc="left", pad=12)
    ax.axhline(1, color="#94A3B8", linewidth=0.6, linestyle="--", zorder=0)
    ax.text(
        len(configs) - 0.5, 1.05,
        "1 = no ambiguity", fontsize=8, color="#64748B",
        ha="right", va="bottom",
    )

    for i, m in enumerate(means):
        ax.text(
            i, m + max(means) * 0.04,
            f"{m:.1f}", ha="center", fontsize=10,
            color=NAVY, fontweight="bold",
        )


def main():
    corpus = generate_corpus(n_each=500, seed=42)
    results = run_sweep(corpus)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    fig.suptitle(
        "Two distinct failure modes of single-CIGAR alignment output  (n=500 per class)",
        fontsize=12, fontweight="bold", y=1.02,
    )

    # Panel A — the category flip
    panel_a = tally(results, "POINT_MUTATION")
    n_a = sum(panel_a[list(panel_a.keys())[0]].values())
    draw_panel(
        axes[0], panel_a,
        title="A.  Reported biology for point-mutation reads",
        correct_category=CORRECT_CATEGORY["POINT_MUTATION"],
        n_cases=n_a,
    )

    # Panel B — inherent ambiguity on homopolymer cases
    draw_ambiguity_panel(
        axes[1], results,
        ground_truth="HOMOPOLYMER_DEL",
        title="B.  Equally-optimal alternatives per homopolymer read",
    )

    # Shared legend (panel A only)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="lower center", bbox_to_anchor=(0.28, -0.12),
        ncol=4, frameon=False, fontsize=9,
        title="panel A — classification of reported CIGAR",
    )

    fig.tight_layout()
    out_path = "figure_category_distribution.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
