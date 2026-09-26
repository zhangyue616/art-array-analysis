#!/usr/bin/env python3
"""Describe RT/partner conservation in the fixed ten Type-II ART loci.

This analysis is deliberately descriptive.  The ten loci form two related
lineages rather than ten independent experimental samples, and no activity is
inferred from catalytic-motif or sequence differences.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "p6_protein" / "python_packages"))
sys.path.insert(0, str(ROOT / "tools" / "python312_packages"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Align import PairwiseAligner  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402


LOCI_TSV = ROOT / "data" / "processed" / "p2_loci" / "fixed_10_loci.tsv"
RT_FASTA = (
    ROOT / "data" / "processed" / "p1_adjust_assignment" / "rt_sequences.faa"
)
PARTNER_FASTA = (
    ROOT
    / "data"
    / "processed"
    / "p1_adjust_assignment"
    / "direct_downstream_partner_sequences.faa"
)
GENPEPT_DIR = ROOT / "data" / "raw" / "p6_protein" / "fixed_proteins"
OUT_DIR = ROOT / "data" / "processed" / "p6_protein"
FIG_DIR = ROOT / "figures" / "p6_protein"

ORDER = [
    "SA1",
    "MarsHill",
    "Madawaska",
    "LY01",
    "S6",
    "PALS_2",
    "UFV_DC4",
    "AH12",
    "Machias",
    "PB50",
]
GROUP = {label: ("seven_locus_reference" if i < 7 else "PB50_related") for i, label in enumerate(ORDER)}
MARS_HILL_CORE_START_1BASED = 211


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    return aligner


ALIGNER = make_aligner()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_label_sequences(path: Path) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    sequences: dict[str, str] = {}
    metadata: dict[str, dict[str, str]] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        bits = rec.id.split("|")
        if len(bits) < 4:
            raise ValueError(f"Unexpected FASTA header: {rec.id}")
        label, accession, protein_id = bits[:3]
        sequences[label] = str(rec.seq).upper()
        metadata[label] = {
            "accession": accession,
            "protein_id": protein_id,
            "source_header": rec.id,
        }
    return sequences, metadata


def strict_yxdd(sequence: str) -> list[tuple[int, str]]:
    hits = []
    for i in range(max(0, len(sequence) - 350), len(sequence) - 3):
        token = sequence[i : i + 4]
        if token[0] == "Y" and token[2:] == "DD":
            hits.append((i + 1, token))
    return hits


def alignment_stats(seq_a: str, seq_b: str) -> dict[str, float | int]:
    aln = ALIGNER.align(seq_a, seq_b)[0]
    idx = np.asarray(aln.indices)
    both = (idx[0] >= 0) & (idx[1] >= 0)
    matches = sum(seq_a[i] == seq_b[j] for i, j in zip(idx[0, both], idx[1, both]))
    columns = int(idx.shape[1])
    aligned_residues = int(both.sum())
    return {
        "matches": int(matches),
        "aligned_residues": aligned_residues,
        "alignment_columns": columns,
        "identity_per_alignment_column": matches / columns if columns else math.nan,
        "identity_per_aligned_residue": matches / aligned_residues if aligned_residues else math.nan,
        "coverage_a": aligned_residues / len(seq_a) if seq_a else math.nan,
        "coverage_b": aligned_residues / len(seq_b) if seq_b else math.nan,
        "score": float(aln.score),
    }


def transfer_reference_boundary(query: str, reference: str, ref_boundary_0: int) -> tuple[int, str]:
    """Map a between-residue reference boundary onto query coordinates."""
    aln = ALIGNER.align(query, reference)[0]
    idx = np.asarray(aln.indices)
    # The first query residue aligned at or after the reference boundary.
    candidates = np.where((idx[1] >= ref_boundary_0) & (idx[0] >= 0))[0]
    if candidates.size:
        q_boundary = int(idx[0, candidates[0]])
        status = "aligned_first_query_residue_at_or_after_reference_boundary"
        return q_boundary, status
    # Defensive fallback: last mapped query residue plus one.
    mapped = np.where(idx[0] >= 0)[0]
    if not mapped.size:
        raise ValueError("Boundary transfer has no mapped query residues")
    return int(idx[0, mapped[-1]]) + 1, "fallback_after_last_mapped_query_residue"


def parse_regions(protein_id: str) -> list[dict[str, str | int]]:
    path = GENPEPT_DIR / f"{protein_id}.gp"
    rec = SeqIO.read(str(path), "genbank")
    rows = []
    for feature in rec.features:
        if feature.type != "Region":
            continue
        rows.append(
            {
                "region_start_1based": int(feature.location.start) + 1,
                "region_end_1based": int(feature.location.end),
                "region_name": "; ".join(feature.qualifiers.get("region_name", [])),
                "region_note": "; ".join(feature.qualifiers.get("note", [])),
                "region_db_xref": "; ".join(feature.qualifiers.get("db_xref", [])),
            }
        )
    return rows


def relation(label_a: str, label_b: str) -> str:
    ga, gb = GROUP[label_a], GROUP[label_b]
    if ga == gb == "seven_locus_reference":
        return "within_seven_locus_reference"
    if ga == gb == "PB50_related":
        return "within_PB50_related"
    return "cross_group"


def fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def median_range(values: list[float]) -> tuple[float, float, float]:
    return statistics.median(values), min(values), max(values)


def plot_figure(
    features: list[dict],
    matrices: dict[str, np.ndarray],
    out_pdf: Path,
    out_png: Path,
) -> None:
    width_in = 170.0 / 25.4
    height_in = 180.0 / 25.4
    palette = {"seven_locus_reference": "#3775BA", "PB50_related": "#9A4D8E"}
    labels = [r["label"] for r in features]
    positions = np.arange(len(labels))
    heat_values = np.concatenate([m[np.triu_indices_from(m, 1)] for m in matrices.values()])
    vmin = float(np.nanmin(heat_values))
    vmax = 1.0

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    ):
        fig = plt.figure(figsize=(width_in, height_in), layout="constrained")
        grid = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0])
        ax_arch = fig.add_subplot(grid[0, :])
        heat_axes = [fig.add_subplot(grid[1, i]) for i in range(3)]

        ntd = np.array([int(r["mapped_ntd_length_aa"]) for r in features])
        core = np.array([int(r["mapped_core_length_aa"]) for r in features])
        total = ntd + core
        colors = [palette[r["subgroup"]] for r in features]
        ax_arch.barh(positions, ntd, color="#D7DDE5", edgecolor="#555555", linewidth=0.35)
        ax_arch.barh(positions, core, left=ntd, color=colors, edgecolor="#555555", linewidth=0.35)
        for y, row, x_total in zip(positions, features, total):
            motif_x = int(row["c_terminal_yxdd_position_1based"]) - 1
            ax_arch.plot(motif_x, y, marker="|", markersize=9, markeredgewidth=1.2, color="#111111")
            ax_arch.text(x_total + 4, y, f"{row['c_terminal_yxdd']} | P {row['partner_length_aa']} aa", va="center", ha="left", fontsize=8)
        ax_arch.axhline(6.5, color="#4A4A4A", lw=0.7, ls="--")
        ax_arch.text(ntd[0] / 2, 0, "NTD", ha="center", va="center", fontsize=8)
        ax_arch.text(ntd[0] + core[0] / 2, 0, "RT core\nreference", ha="center", va="center", fontsize=8, color="white")
        ax_arch.text(ntd[7] + core[7] / 2, 7, "RT core\nPB50-related", ha="center", va="center", fontsize=8, color="white")
        ax_arch.set_yticks(positions, labels)
        ax_arch.invert_yaxis()
        ax_arch.set_xlim(0, max(total) + 78)
        ax_arch.set_xlabel("RT residue position; P denotes direct downstream partner length")
        ax_arch.set_title("a  RT architecture transferred from the MarsHill core boundary", loc="left", fontweight="bold")
        ax_arch.spines[["top", "right"]].set_visible(False)
        ax_arch.grid(axis="x", color="#E5E5E5", lw=0.45)

        cmap = plt.get_cmap("viridis")
        norm = Normalize(vmin=vmin, vmax=vmax)
        heat_specs = [
            ("RT N-terminal segment", "rt_ntd_identity"),
            ("RT core segment", "rt_core_identity"),
            ("Partner protein", "partner_identity"),
        ]
        image = None
        for panel, (title, key), ax in zip("bcd", heat_specs, heat_axes):
            image = ax.imshow(matrices[key], cmap=cmap, norm=norm, interpolation="nearest")
            ax.axhline(6.5, color="white", lw=1.0)
            ax.axvline(6.5, color="white", lw=1.0)
            ax.set_xticks(range(10), labels, rotation=90)
            ax.set_yticks(range(10), labels if panel == "b" else [""] * 10)
            ax.set_title(f"{panel}  {title}", loc="left", fontweight="bold")
            ax.tick_params(length=0, pad=1.5)
        assert image is not None
        cbar = fig.colorbar(image, ax=heat_axes, orientation="horizontal", fraction=0.05, pad=0.11, aspect=35)
        cbar.set_label("Global pairwise identity per alignment column")
        cbar.ax.tick_params(labelsize=8)

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, dpi=300)
        fig.savefig(out_png, dpi=300)
        plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    locus_rows = {r["label"]: r for r in read_tsv(LOCI_TSV)}
    rt, rt_meta = load_label_sequences(RT_FASTA)
    partner, partner_meta = load_label_sequences(PARTNER_FASTA)

    if set(rt) != set(ORDER) or set(partner) != set(ORDER) or set(locus_rows) != set(ORDER):
        raise ValueError("The fixed ten labels do not agree across inputs")

    mars_hill_rt = rt["MarsHill"]
    ref_boundary_0 = MARS_HILL_CORE_START_1BASED - 1
    features: list[dict] = []
    region_rows: list[dict] = []
    boundaries: dict[str, int] = {}
    for label in ORDER:
        qseq = rt[label]
        q_boundary, boundary_status = transfer_reference_boundary(qseq, mars_hill_rt, ref_boundary_0)
        boundaries[label] = q_boundary
        motif_hits = strict_yxdd(qseq)
        if len(motif_hits) != 1:
            raise ValueError(f"{label}: expected one C-terminal YxDD, found {motif_hits}")
        motif_pos, motif = motif_hits[0]
        rt_regions = parse_regions(rt_meta[label]["protein_id"])
        partner_regions = parse_regions(partner_meta[label]["protein_id"])
        for role, protein_id, rows in (
            ("RT", rt_meta[label]["protein_id"], rt_regions),
            ("partner", partner_meta[label]["protein_id"], partner_regions),
        ):
            if rows:
                for row in rows:
                    region_rows.append({"label": label, "role": role, "protein_id": protein_id, **row})
            else:
                region_rows.append(
                    {
                        "label": label,
                        "role": role,
                        "protein_id": protein_id,
                        "region_start_1based": "",
                        "region_end_1based": "",
                        "region_name": "",
                        "region_note": "no Region feature in retrieved GenPept record",
                        "region_db_xref": "",
                    }
                )
        lr = locus_rows[label]
        features.append(
            {
                "locus_order": lr["locus_order"],
                "label": label,
                "subgroup": GROUP[label],
                "genome_accession": lr["accession"],
                "assignment_status": lr["assignment_status"],
                "rt_protein_id": rt_meta[label]["protein_id"],
                "rt_length_aa": len(qseq),
                "partner_protein_id": partner_meta[label]["protein_id"],
                "partner_length_aa": len(partner[label]),
                "rt_partner_gap_nt": lr["rt_partner_gap_nt"],
                "mapped_ntd_length_aa": q_boundary,
                "mapped_core_start_1based": q_boundary + 1,
                "mapped_core_length_aa": len(qseq) - q_boundary,
                "boundary_transfer_status": boundary_status,
                "boundary_reference": f"MarsHill residue {MARS_HILL_CORE_START_1BASED}; source-defined core start",
                "c_terminal_yxdd": motif,
                "c_terminal_yxdd_position_1based": motif_pos,
                "rt_genpept_region_count": len(rt_regions),
                "partner_genpept_region_count": len(partner_regions),
                "interpretation_limit": "sequence feature only; no catalytic activity or causal effect inferred",
            }
        )

    write_tsv(OUT_DIR / "protein_feature_summary.tsv", features)
    write_tsv(OUT_DIR / "genpept_region_features.tsv", region_rows)

    pairs: list[dict] = []
    matrix_keys = ["rt_full_identity", "rt_ntd_identity", "rt_core_identity", "partner_identity"]
    matrices = {key: np.eye(len(ORDER), dtype=float) for key in matrix_keys}
    for i, label_a in enumerate(ORDER):
        for j in range(i + 1, len(ORDER)):
            label_b = ORDER[j]
            parts = {
                "rt_full": (rt[label_a], rt[label_b]),
                "rt_ntd": (rt[label_a][: boundaries[label_a]], rt[label_b][: boundaries[label_b]]),
                "rt_core": (rt[label_a][boundaries[label_a] :], rt[label_b][boundaries[label_b] :]),
                "partner": (partner[label_a], partner[label_b]),
            }
            stats = {name: alignment_stats(*seqs) for name, seqs in parts.items()}
            row = {
                "label_a": label_a,
                "label_b": label_b,
                "relation": relation(label_a, label_b),
            }
            for part_name, result in stats.items():
                for metric, value in result.items():
                    row[f"{part_name}_{metric}"] = fmt(value) if isinstance(value, float) else value
                key = f"{part_name}_identity"
                identity = float(result["identity_per_alignment_column"])
                matrices[key][i, j] = identity
                matrices[key][j, i] = identity
            row["independence_note"] = "descriptive pair; pairs share loci and are not independent samples"
            pairs.append(row)
    write_tsv(OUT_DIR / "protein_pairwise_identity.tsv", pairs)

    summary_rows: list[dict] = []
    for rel in ("within_seven_locus_reference", "within_PB50_related", "cross_group"):
        selected = [r for r in pairs if r["relation"] == rel]
        for metric in ("rt_full", "rt_ntd", "rt_core", "partner"):
            values = [float(r[f"{metric}_identity_per_alignment_column"]) for r in selected]
            med, low, high = median_range(values)
            summary_rows.append(
                {
                    "comparison_set": rel,
                    "metric": f"{metric}_global_identity_per_alignment_column",
                    "dependent_pair_count": len(values),
                    "median": fmt(med),
                    "minimum": fmt(low),
                    "maximum": fmt(high),
                    "independent_unit_note": "pairwise values are dependent because each locus occurs in multiple pairs",
                }
            )
    write_tsv(OUT_DIR / "group_pairwise_summary.tsv", summary_rows)

    feature_group_rows: list[dict] = []
    for subgroup in ("seven_locus_reference", "PB50_related"):
        selected = [r for r in features if r["subgroup"] == subgroup]
        for metric in ("rt_length_aa", "mapped_ntd_length_aa", "mapped_core_length_aa", "partner_length_aa"):
            values = [int(r[metric]) for r in selected]
            feature_group_rows.append(
                {
                    "subgroup": subgroup,
                    "locus_count": len(selected),
                    "metric": metric,
                    "median": fmt(statistics.median(values), 1),
                    "minimum": min(values),
                    "maximum": max(values),
                    "independent_lineage_note": "loci are closely related and do not establish ten independent evolutionary replicates",
                }
            )
        motifs = Counter(r["c_terminal_yxdd"] for r in selected)
        feature_group_rows.append(
            {
                "subgroup": subgroup,
                "locus_count": len(selected),
                "metric": "c_terminal_yxdd_counts",
                "median": "",
                "minimum": "",
                "maximum": ";".join(f"{k}:{v}" for k, v in sorted(motifs.items())),
                "independent_lineage_note": "motif x-residue split is descriptive; activity consequence is untested",
            }
        )
    write_tsv(OUT_DIR / "group_feature_summary.tsv", feature_group_rows)

    plot_figure(
        features,
        {
            "rt_ntd_identity": matrices["rt_ntd_identity"],
            "rt_core_identity": matrices["rt_core_identity"],
            "partner_identity": matrices["partner_identity"],
        },
        FIG_DIR / "figure1_typeII_protein_conservation.pdf",
        FIG_DIR / "figure1_typeII_protein_conservation.png",
    )

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "biopython": __import__("Bio").__version__,
        "matplotlib": matplotlib.__version__,
        "parameters": {
            "aligner": {
                "mode": ALIGNER.mode,
                "match_score": ALIGNER.match_score,
                "mismatch_score": ALIGNER.mismatch_score,
                "open_gap_score": ALIGNER.open_gap_score,
                "extend_gap_score": ALIGNER.extend_gap_score,
            },
            "identity_denominator": "global alignment columns including gap columns",
            "boundary_reference": f"MarsHill residue {MARS_HILL_CORE_START_1BASED}",
            "yxdd_scan": "single YxDD in final 350 aa",
        },
        "figure": {
            "width_mm": 170,
            "height_mm": 180,
            "minimum_source_font_pt": 8,
            "png_dpi": 300,
            "pdf_vector": True,
            "heatmap_color_limits": "full observed off-diagonal minimum through 1.0",
        },
        "validation": {
            "fixed_locus_count": len(features),
            "pair_count": len(pairs),
            "expected_pair_count": 45,
            "all_have_one_c_terminal_yxdd": all(len(strict_yxdd(rt[x])) == 1 for x in ORDER),
            "partner_records_without_region_feature": sum(int(r["partner_genpept_region_count"]) == 0 for r in features),
            "pdf_exists": (FIG_DIR / "figure1_typeII_protein_conservation.pdf").exists(),
            "png_exists": (FIG_DIR / "figure1_typeII_protein_conservation.png").exists(),
            "status": "PASS" if len(features) == 10 and len(pairs) == 45 else "FAIL",
        },
        "outputs": [
            "data/processed/p6_protein/protein_feature_summary.tsv",
            "data/processed/p6_protein/genpept_region_features.tsv",
            "data/processed/p6_protein/protein_pairwise_identity.tsv",
            "data/processed/p6_protein/group_pairwise_summary.tsv",
            "data/processed/p6_protein/group_feature_summary.tsv",
            "figures/p6_protein/figure1_typeII_protein_conservation.pdf",
            "figures/p6_protein/figure1_typeII_protein_conservation.png",
        ],
    }
    (OUT_DIR / "protein_analysis_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(receipt["validation"], indent=2))


if __name__ == "__main__":
    main()
