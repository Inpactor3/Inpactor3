"""
Calcula qué porcentaje de un genoma fue anotado por Inpactor3 como LTR-RT.

Lee un archivo `Inpactor3_predictions.tab` y el FASTA del genoma original;
computa:
  - bases totales del genoma
  - bases cubiertas por predicciones (unión de intervalos, no doble conteo)
  - % del genoma cubierto (global y por linaje)
  - compara contra el valor esperado biológicamente (--expected)

Uso:
    python scripts/genome_coverage.py \\
        --genome data/genomes/Arabidopsis_thaliana.TAIR10.dna_sm.chromosome.1.fa \\
        --predictions results/arabidopsis/Inpactor3_predictions.tab \\
        --expected 10.0 \\
        --report results/arabidopsis/coverage_report.md
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def scaffold_lengths(fasta: Path) -> dict[str, int]:
    lens: dict[str, int] = {}
    cur, n = None, 0
    with open(fasta) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur is not None:
                    lens[cur] = n
                cur = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
        if cur is not None:
            lens[cur] = n
    return lens


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Une intervalos que se solapan para no doble-contar."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def load_predictions(tab: Path):
    """Lee .tab de 8 columnas. Devuelve {scaffold: [(start, end, linaje)]}."""
    by_scaf: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with open(tab) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue
            scaf, s, e, _len, lin = cols[0], int(cols[1]), int(cols[2]), int(cols[3]), cols[4]
            by_scaf[scaf].append((s, e, lin))
    return by_scaf


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, type=Path)
    ap.add_argument("--predictions", required=True, type=Path)
    ap.add_argument("--expected", type=float, default=None,
                    help="% LTR-RT esperado (ej. 10.0 para Arabidopsis)")
    ap.add_argument("--report", type=Path, required=True)
    args = ap.parse_args()

    lens = scaffold_lengths(args.genome)
    total_bp = sum(lens.values())
    print(f"[genoma] {args.genome.name} · {total_bp:,} bp · {len(lens)} scaffold(s)")

    preds = load_predictions(args.predictions)
    n_pred = sum(len(v) for v in preds.values())
    print(f"[preds ] {n_pred:,} predicciones en {len(preds)} scaffold(s)")

    # Cobertura global (unión de intervalos)
    covered_bp = 0
    for scaf, hits in preds.items():
        merged = merge_intervals([(h[0], h[1]) for h in hits])
        covered_bp += sum(e - s for s, e in merged)
    pct = 100 * covered_bp / max(total_bp, 1)

    # Cobertura por linaje
    by_lin_cov: dict[str, int] = {}
    lin_counts: dict[str, int] = defaultdict(int)
    lin_bp_raw: dict[str, int] = defaultdict(int)
    for scaf, hits in preds.items():
        by_lin: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for s, e, lin in hits:
            by_lin[lin].append((s, e))
            lin_counts[lin] += 1
            lin_bp_raw[lin] += e - s
        for lin, intervals in by_lin.items():
            merged = merge_intervals(intervals)
            by_lin_cov[lin] = by_lin_cov.get(lin, 0) + sum(e - s for s, e in merged)

    # Reporte
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w") as f:
        f.write(f"# Cobertura genómica — {args.genome.name}\n\n")
        f.write(f"- Bases del genoma: **{total_bp:,}**\n")
        f.write(f"- Predicciones totales: **{n_pred:,}**\n")
        f.write(f"- Bases cubiertas (unión): **{covered_bp:,}**\n")
        f.write(f"- **% del genoma anotado como LTR-RT: {pct:.2f}%**\n")
        if args.expected is not None:
            diff = pct - args.expected
            veredicto = ("✅ dentro del rango" if abs(diff) < 5 else
                         "⚠️  fuera del rango esperado")
            f.write(f"- % esperado biológicamente: **{args.expected:.1f}%**\n")
            f.write(f"- Diferencia: **{diff:+.2f} puntos porcentuales**  {veredicto}\n")
        f.write("\n## Cobertura por linaje\n\n")
        f.write("| Linaje | Predicciones | Bases cubiertas | % del genoma |\n")
        f.write("|---|---:|---:|---:|\n")
        for lin, cov in sorted(by_lin_cov.items(), key=lambda x: -x[1]):
            f.write(f"| {lin} | {lin_counts[lin]:,} | {cov:,} | "
                    f"{100*cov/total_bp:.2f}% |\n")

    print(f"\n[out ] {args.report}")
    print(f"\n{'='*50}")
    print(f"  % del genoma anotado como LTR-RT: {pct:.2f}%")
    if args.expected is not None:
        diff = pct - args.expected
        print(f"  Esperado: {args.expected:.1f}%  (diferencia: {diff:+.2f} pp)")
    print("=" * 50)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
