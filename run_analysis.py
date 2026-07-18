"""
run_analysis.py

CLI entry point: parse a PSSM ASN.1 text file, run the conservation
analysis, and write out a CSV, two plots, and a markdown report.

Usage:
    python3 run_analysis.py <path_to_pssm.asn> <output_dir>
"""

import sys
import os
import csv
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pssm_parser import parse_pssm_asn_text, STANDARD_AA
from pssm_analyzer import analyze_pssm, summary_narrative, predict_mutation, FUNCTIONAL_HINTS


def write_csv(results, out_path):
    fields = ["position", "wild_residue", "wild_score", "information_content_bits",
              "conservation", "best_residue", "best_score", "second_best_residue",
              "second_best_score", "worst_residue", "worst_score"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in results:
            writer.writerow(r)


def plot_conservation_profile(results, out_path, title):
    positions = [r["position"] for r in results]
    ic = [r["information_content_bits"] for r in results]

    fig, ax = plt.subplots(figsize=(14, 3.5))
    ax.fill_between(positions, ic, step="mid", alpha=0.35, color="#3b6fa0")
    ax.plot(positions, ic, linewidth=0.8, color="#1f4e79")
    ax.axhline(2.0, color="crimson", linestyle="--", linewidth=0.8, label="highly conserved threshold")
    ax.axhline(1.0, color="orange", linestyle="--", linewidth=0.8, label="moderately conserved threshold")
    ax.set_xlabel("Position")
    ax.set_ylabel("Information content (bits)")
    ax.set_title(f"Conservation profile -- {title}")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlim(min(positions), max(positions))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_score_heatmap(parsed, out_path, title, max_positions=261):
    n = min(parsed.num_columns, max_positions)
    matrix = np.array([[parsed.scores[c][aa] for c in range(n)] for aa in STANDARD_AA])

    fig, ax = plt.subplots(figsize=(max(10, n * 0.05), 6))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-8, vmax=8)
    ax.set_yticks(range(len(STANDARD_AA)))
    ax.set_yticklabels(STANDARD_AA, fontsize=7)
    ax.set_xlabel("Position")
    ax.set_title(f"PSSM score heatmap -- {title}")
    fig.colorbar(im, ax=ax, shrink=0.7, label="log-odds score")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_report(parsed, results, out_path, csv_name, profile_png, heatmap_png):
    narrative = summary_narrative(parsed, results)
    top10 = sorted(results, key=lambda r: r["information_content_bits"], reverse=True)[:10]

    lines = []
    lines.append(f"# PSSM Interpretation Report\n")
    lines.append(f"**Query:** {parsed.query_id}  ")
    lines.append(f"**Title:** {parsed.query_title}  ")
    lines.append(f"**Length:** {parsed.num_columns} residues  ")
    lines.append(f"**Karlin-Altschul kappa:** {parsed.kappa}\n")
    lines.append("## Summary\n")
    lines.append(narrative + "\n")
    lines.append(f"![Conservation profile]({os.path.basename(profile_png)})\n")
    lines.append(f"![Score heatmap]({os.path.basename(heatmap_png)})\n")
    lines.append("## Top 10 most conserved positions\n")
    lines.append("| Position | Wild residue | IC (bits) | Conservation | Best alt. | Worst alt. |")
    lines.append("|---|---|---|---|---|---|")
    for r in top10:
        lines.append(
            f"| {r['position']} | {r['wild_residue']} | {r['information_content_bits']} "
            f"| {r['conservation']} | {r['best_residue']} ({r['best_score']}) "
            f"| {r['worst_residue']} ({r['worst_score']}) |"
        )
    lines.append(f"\nFull per-position data: `{csv_name}`\n")
    lines.append("## Interpretation notes\n")
    lines.append(FUNCTIONAL_HINTS + "\n")
    lines.append(
        "**Scope note:** this report is derived entirely from the PSSM's own scores "
        "and frequency data. It does not (yet) cross-reference structure, domain "
        "databases (Pfam/CDD), or literature -- those require separate, verifiable "
        "integrations and are intentionally out of scope for this version so nothing "
        "here overstates what the profile alone can support."
    )

    with open(out_path, "w") as f:
        f.write("\n".join(lines))


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 run_analysis.py <pssm.asn> <output_dir>")
        sys.exit(1)

    pssm_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    parsed = parse_pssm_asn_text(pssm_path)
    results = analyze_pssm(parsed)

    title = parsed.query_title or parsed.query_id

    csv_path = os.path.join(out_dir, "pssm_analysis.csv")
    json_path = os.path.join(out_dir, "pssm_analysis.json")
    profile_png = os.path.join(out_dir, "conservation_profile.png")
    heatmap_png = os.path.join(out_dir, "score_heatmap.png")
    report_path = os.path.join(out_dir, "report.md")

    write_csv(results, csv_path)
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    plot_conservation_profile(results, profile_png, title)
    plot_score_heatmap(parsed, heatmap_png, title)
    write_report(parsed, results, report_path, os.path.basename(csv_path), profile_png, heatmap_png)

    print(f"Wrote: {csv_path}\n       {json_path}\n       {profile_png}\n       {heatmap_png}\n       {report_path}")


if __name__ == "__main__":
    main()
