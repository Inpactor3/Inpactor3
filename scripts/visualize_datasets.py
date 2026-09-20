"""
Genera visualizaciones del corpus fusionado y de la auditoria de fuentes.

Salida: docs/plots/*.png
  1) lineage_distribution.png  — barras horizontales por linaje
  2) top_species.png           — top 20 especies mas representadas
  3) length_histogram.png      — histograma de longitud de LTR-RTs
  4) lineage_x_species_heatmap.png — heatmap species x linaje (top-15 x top-15)
  5) resumen.png               — panel 2x2 con las anteriores

Uso:
    python scripts/visualize_datasets.py
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "data" / "corpus" / "inpactor3_v0.fasta"
MANIFEST = ROOT / "data" / "corpus" / "inpactor3_v0.manifest.jsonl"
OUT = ROOT / "docs" / "plots"
OUT.mkdir(parents=True, exist_ok=True)

# Colores consistentes por superfamilia
COPIA_LIN = {"ALE", "ANGELA", "BIANCA", "IKEROS", "IVANA", "ORYCO", "OSSER", "RETROFIT", "SIRE", "TAR", "TORK"}
GYPSY_LIN = {"ATHILA", "TAT", "DEL", "TEKAY", "CRM", "GALADRIEL", "REINA", "CHROMOVIRUS"}
COLOR_COPIA = "#3b82f6"   # azul
COLOR_GYPSY = "#ef4444"   # rojo
COLOR_OTHER = "#94a3b8"


def lineage_color(lin: str) -> str:
    if lin in COPIA_LIN:
        return COLOR_COPIA
    if lin in GYPSY_LIN:
        return COLOR_GYPSY
    return COLOR_OTHER


def read_manifest() -> list[dict]:
    if not MANIFEST.exists():
        raise SystemExit(f"No existe {MANIFEST}. Corre scripts/merge_corpus.py primero.")
    with open(MANIFEST) as f:
        return [json.loads(l) for l in f]


def plot_lineages(records: list[dict]) -> None:
    counts = Counter(r["canonical_lineage"] for r in records)
    lins, vals = zip(*sorted(counts.items(), key=lambda x: x[1]))
    colors = [lineage_color(l) for l in lins]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(lins, vals, color=colors, edgecolor="white", linewidth=0.5)
    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.005, i, f"{v:,}  ({v/sum(vals)*100:.1f}%)",
                va="center", fontsize=9)
    ax.set_xlabel("Número de LTR-retrotransposones")
    ax.set_title(f"Distribución por linaje canónico — {sum(vals):,} secuencias")
    ax.set_xlim(0, max(vals) * 1.18)
    # Leyenda
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color=COLOR_COPIA, label="Copia (RLC)"),
        Patch(color=COLOR_GYPSY, label="Gypsy (RLG)"),
    ], loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "lineage_distribution.png", dpi=140)
    plt.close(fig)
    print("[ok]", OUT / "lineage_distribution.png")


def plot_top_species(records: list[dict], top: int = 20) -> None:
    counts = Counter(r["species"] for r in records)
    species_top = counts.most_common(top)
    sps, vals = zip(*species_top[::-1])

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(sps, vals, color="#059669", edgecolor="white", linewidth=0.5)
    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.005, i, f"{v:,}", va="center", fontsize=9)
    ax.set_xlabel("Número de LTR-RTs contribuidos")
    ax.set_title(f"Top {top} especies del corpus")
    ax.set_xlim(0, max(vals) * 1.15)
    fig.tight_layout()
    fig.savefig(OUT / "top_species.png", dpi=140)
    plt.close(fig)
    print("[ok]", OUT / "top_species.png")


def plot_length_histogram(records: list[dict]) -> None:
    lens = np.array([r["length"] for r in records])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(lens, bins=60, color="#8b5cf6", edgecolor="white")
    ax.axvline(np.median(lens), color="#dc2626", linestyle="--", label=f"mediana = {int(np.median(lens))} bp")
    ax.axvline(np.mean(lens), color="#0891b2", linestyle="--", label=f"media = {int(np.mean(lens))} bp")
    ax.set_xlabel("Longitud (bp)")
    ax.set_ylabel("Frecuencia")
    ax.set_title(f"Distribución de longitud de LTR-RTs (n = {len(lens):,})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "length_histogram.png", dpi=140)
    plt.close(fig)
    print("[ok]", OUT / "length_histogram.png")


def plot_heatmap(records: list[dict], top_species: int = 15, top_lineages: int = 13) -> None:
    sp_counts = Counter(r["species"] for r in records)
    lin_counts = Counter(r["canonical_lineage"] for r in records)
    top_sp = [s for s, _ in sp_counts.most_common(top_species)]
    top_lin = [l for l, _ in lin_counts.most_common(top_lineages)]

    matrix = np.zeros((len(top_sp), len(top_lin)), dtype=int)
    sp_idx = {s: i for i, s in enumerate(top_sp)}
    lin_idx = {l: j for j, l in enumerate(top_lin)}
    for r in records:
        s, l = r["species"], r["canonical_lineage"]
        if s in sp_idx and l in lin_idx:
            matrix[sp_idx[s], lin_idx[l]] += 1

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(top_lin)))
    ax.set_xticklabels(top_lin, rotation=45, ha="right")
    ax.set_yticks(range(len(top_sp)))
    ax.set_yticklabels(top_sp)
    for i in range(len(top_sp)):
        for j in range(len(top_lin)):
            v = matrix[i, j]
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center",
                        color="white" if v > matrix.max() * 0.5 else "black", fontsize=8)
    ax.set_title("Especies × Linajes — conteos en el corpus")
    fig.colorbar(im, ax=ax, label="nº LTR-RTs")
    fig.tight_layout()
    fig.savefig(OUT / "lineage_x_species_heatmap.png", dpi=140)
    plt.close(fig)
    print("[ok]", OUT / "lineage_x_species_heatmap.png")


def plot_summary(records: list[dict]) -> None:
    """Panel 2x2 combinando los graficos principales."""
    lens = np.array([r["length"] for r in records])
    lin_counts = Counter(r["canonical_lineage"] for r in records)
    sp_counts = Counter(r["species"] for r in records)

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    fig.suptitle(
        f"Corpus Inpactor3 v0 — {len(records):,} secuencias · "
        f"{len(sp_counts)} especies · {len(lin_counts)} linajes canónicos",
        fontsize=15, fontweight="bold"
    )

    # (0,0) linajes
    lins, vals = zip(*sorted(lin_counts.items(), key=lambda x: x[1]))
    axes[0, 0].barh(lins, vals, color=[lineage_color(l) for l in lins])
    axes[0, 0].set_title("Distribución por linaje canónico")
    axes[0, 0].set_xlabel("nº LTR-RTs")

    # (0,1) top especies
    top_sp = sp_counts.most_common(15)
    sps, svals = zip(*top_sp[::-1])
    axes[0, 1].barh(sps, svals, color="#059669")
    axes[0, 1].set_title("Top 15 especies")
    axes[0, 1].set_xlabel("nº LTR-RTs")

    # (1,0) histograma de longitud
    axes[1, 0].hist(lens, bins=60, color="#8b5cf6", edgecolor="white")
    axes[1, 0].axvline(np.median(lens), color="#dc2626", linestyle="--",
                       label=f"mediana {int(np.median(lens))} bp")
    axes[1, 0].axvline(np.mean(lens), color="#0891b2", linestyle="--",
                       label=f"media {int(np.mean(lens))} bp")
    axes[1, 0].set_title("Longitud de los LTR-RTs")
    axes[1, 0].set_xlabel("bp")
    axes[1, 0].legend()

    # (1,1) heatmap
    top_lin = [l for l, _ in lin_counts.most_common(13)]
    top_sp = [s for s, _ in sp_counts.most_common(15)]
    matrix = np.zeros((len(top_sp), len(top_lin)), dtype=int)
    sp_idx = {s: i for i, s in enumerate(top_sp)}
    lin_idx = {l: j for j, l in enumerate(top_lin)}
    for r in records:
        if r["species"] in sp_idx and r["canonical_lineage"] in lin_idx:
            matrix[sp_idx[r["species"]], lin_idx[r["canonical_lineage"]]] += 1
    im = axes[1, 1].imshow(matrix, cmap="YlOrRd", aspect="auto")
    axes[1, 1].set_xticks(range(len(top_lin)))
    axes[1, 1].set_xticklabels(top_lin, rotation=45, ha="right", fontsize=8)
    axes[1, 1].set_yticks(range(len(top_sp)))
    axes[1, 1].set_yticklabels(top_sp, fontsize=8)
    axes[1, 1].set_title("Especies × Linajes (heatmap)")
    fig.colorbar(im, ax=axes[1, 1], label="nº LTR-RTs", fraction=0.04)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "resumen.png", dpi=140)
    plt.close(fig)
    print("[ok]", OUT / "resumen.png")


def main() -> int:
    print(f"[read] {MANIFEST}")
    records = read_manifest()
    print(f"[ok  ] {len(records):,} registros")
    plot_lineages(records)
    plot_top_species(records)
    plot_length_histogram(records)
    plot_heatmap(records)
    plot_summary(records)
    print(f"\n[done] gráficos en {OUT.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
