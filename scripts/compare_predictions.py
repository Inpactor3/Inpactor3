"""
Compara predicciones de Inpactor2 e Inpactor3 sobre el MISMO genoma.

Ambos programas deben haber producido su archivo `.tab` con formato:
    seqid  start  end  length  lineage  detect_prob  filter_prob  class_prob

Métricas:
  - Emparejamiento por IoU 1D ≥ umbral (default 0.5)
  - Por linaje: precisión, recall, F1
  - Confusión global de linajes
  - Concordancia agregada (¿los dos ven los mismos elementos?)

Uso:
    python scripts/compare_predictions.py \
        --inpactor2 results/inpactor2/Inpactor2_predictions.tab \
        --inpactor3 results/inpactor3/Inpactor3_predictions.tab \
        --iou 0.5 \
        --report results/comparison.md
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def load_tab(path: Path) -> list[dict]:
    """Lee un .tab de 8 columnas y devuelve lista de dicts."""
    hits = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            hits.append({
                "seqid": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "length": int(cols[3]),
                "lineage": cols[4].upper(),
                "det_prob": float(cols[5]) if len(cols) > 5 and cols[5] != "-" else None,
                "filter_prob": None if len(cols) < 7 or cols[6] == "-" else float(cols[6]),
                "cls_prob": float(cols[7]) if len(cols) > 7 and cols[7] != "-" else None,
            })
    return hits


def iou_1d(a: dict, b: dict) -> float:
    if a["seqid"] != b["seqid"]:
        return 0.0
    s = max(a["start"], b["start"])
    e = min(a["end"], b["end"])
    inter = max(0, e - s)
    union = (a["end"] - a["start"]) + (b["end"] - b["start"]) - inter
    return inter / union if union > 0 else 0.0


def match(refs: list[dict], preds: list[dict], iou_thr: float) -> tuple[list[tuple[int, int, float]], set[int], set[int]]:
    """Empareja greedy por mayor IoU. Devuelve (matches, unmatched_ref, unmatched_pred)."""
    # índice por seqid para no comparar todo contra todo
    ref_by_sq = defaultdict(list)
    for i, r in enumerate(refs):
        ref_by_sq[r["seqid"]].append(i)

    matches: list[tuple[int, int, float]] = []
    used_ref: set[int] = set()
    used_pred: set[int] = set()

    for j, p in enumerate(preds):
        best_i, best_iou = -1, iou_thr
        for i in ref_by_sq.get(p["seqid"], []):
            if i in used_ref:
                continue
            v = iou_1d(refs[i], p)
            if v >= best_iou:
                best_iou = v
                best_i = i
        if best_i >= 0:
            matches.append((best_i, j, best_iou))
            used_ref.add(best_i)
            used_pred.add(j)

    unmatched_ref = set(range(len(refs))) - used_ref
    unmatched_pred = set(range(len(preds))) - used_pred
    return matches, unmatched_ref, unmatched_pred


def per_lineage_prf(refs, preds, matches, unmatched_ref, unmatched_pred) -> dict:
    """Calcula precision/recall/F1 por linaje. TP requiere IoU>=thr Y linaje coincidente."""
    lineages = sorted({r["lineage"] for r in refs} | {p["lineage"] for p in preds})
    tp = defaultdict(int)
    fp = defaultdict(int)
    fn = defaultdict(int)

    matched_ref = set()
    matched_pred = set()
    for i, j, _ in matches:
        if refs[i]["lineage"] == preds[j]["lineage"]:
            tp[refs[i]["lineage"]] += 1
        else:
            fp[preds[j]["lineage"]] += 1  # linaje mal predicho
            fn[refs[i]["lineage"]] += 1
        matched_ref.add(i)
        matched_pred.add(j)

    for i in unmatched_ref:
        fn[refs[i]["lineage"]] += 1
    for j in unmatched_pred:
        fp[preds[j]["lineage"]] += 1

    metrics = {}
    for lin in lineages:
        p = tp[lin] / (tp[lin] + fp[lin]) if (tp[lin] + fp[lin]) else 0.0
        r = tp[lin] / (tp[lin] + fn[lin]) if (tp[lin] + fn[lin]) else 0.0
        f1 = 2 * p * r / (p + r) if (p + r) else 0.0
        metrics[lin] = {"TP": tp[lin], "FP": fp[lin], "FN": fn[lin],
                        "P": p, "R": r, "F1": f1}
    return metrics


def summary(name: str, hits: list[dict]) -> str:
    n = len(hits)
    lins = defaultdict(int)
    for h in hits:
        lins[h["lineage"]] += 1
    lengths = [h["length"] for h in hits]
    mean_l = sum(lengths) / n if n else 0
    lines = [f"### {name}",
             f"- Predicciones totales: **{n}**",
             f"- Longitud media: {mean_l:.0f} bp",
             f"- Linajes distintos: {len(lins)}",
             "- Top linajes:"]
    for lin, c in sorted(lins.items(), key=lambda x: -x[1])[:8]:
        lines.append(f"    - {lin:12s} {c:>6}")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inpactor2", required=True, type=Path)
    ap.add_argument("--inpactor3", required=True, type=Path)
    ap.add_argument("--iou", type=float, default=0.5)
    ap.add_argument("--report", type=Path, default=Path("results/comparison.md"))
    ap.add_argument("--reference", choices=["inpactor2", "inpactor3"], default="inpactor2",
                    help="Cuál se toma como 'verdad'; el otro se evalúa contra él.")
    args = ap.parse_args()

    inp2 = load_tab(args.inpactor2)
    inp3 = load_tab(args.inpactor3)
    print(f"[inp2] {len(inp2)} predicciones")
    print(f"[inp3] {len(inp3)} predicciones")

    refs = inp2 if args.reference == "inpactor2" else inp3
    preds = inp3 if args.reference == "inpactor2" else inp2
    ref_name = args.reference
    pred_name = "inpactor3" if ref_name == "inpactor2" else "inpactor2"

    matches, un_r, un_p = match(refs, preds, args.iou)
    metrics = per_lineage_prf(refs, preds, matches, un_r, un_p)

    # Confusión global de linaje entre las cajas que sí solaparon
    conf = defaultdict(lambda: defaultdict(int))
    for i, j, _ in matches:
        conf[refs[i]["lineage"]][preds[j]["lineage"]] += 1

    # ---- reporte ----
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w") as f:
        f.write(f"# Comparación Inpactor2 vs Inpactor3\n\n")
        f.write(f"- Referencia: **{ref_name}** ({len(refs)} cajas)\n")
        f.write(f"- Evaluado : **{pred_name}** ({len(preds)} cajas)\n")
        f.write(f"- Umbral IoU 1D: **{args.iou}**\n")
        f.write(f"- Cajas que coinciden posicionalmente: **{len(matches)}**\n")
        f.write(f"- Sólo en referencia (FN): **{len(un_r)}**\n")
        f.write(f"- Sólo en {pred_name} (FP): **{len(un_p)}**\n\n")

        f.write("## Resumen por herramienta\n\n")
        f.write(summary("Inpactor2", inp2) + "\n\n")
        f.write(summary("Inpactor3", inp3) + "\n\n")

        f.write("## Métricas por linaje (TP requiere IoU y linaje correctos)\n\n")
        f.write("| Linaje | TP | FP | FN | Precisión | Recall | F1 |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for lin in sorted(metrics):
            m = metrics[lin]
            f.write(
                f"| {lin} | {m['TP']} | {m['FP']} | {m['FN']} | "
                f"{m['P']:.3f} | {m['R']:.3f} | {m['F1']:.3f} |\n"
            )

        macro_p = sum(m["P"] for m in metrics.values()) / max(len(metrics), 1)
        macro_r = sum(m["R"] for m in metrics.values()) / max(len(metrics), 1)
        macro_f1 = sum(m["F1"] for m in metrics.values()) / max(len(metrics), 1)
        f.write(f"\n**Macro promedio**: P={macro_p:.3f}  R={macro_r:.3f}  F1={macro_f1:.3f}\n\n")

        f.write("## Confusión de linaje (filas = referencia, columnas = predicho)\n\n")
        all_lin = sorted(set(conf.keys()) | {k for c in conf.values() for k in c})
        f.write("| ref \\ pred | " + " | ".join(all_lin) + " |\n")
        f.write("|" + "---|" * (len(all_lin) + 1) + "\n")
        for r in all_lin:
            row = [str(conf[r][c]) for c in all_lin]
            f.write(f"| **{r}** | " + " | ".join(row) + " |\n")

    print(f"\n[out] {args.report}")
    print(f"\nMétricas macro: P={macro_p:.3f}  R={macro_r:.3f}  F1={macro_f1:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
