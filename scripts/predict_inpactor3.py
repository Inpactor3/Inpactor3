"""
Inferencia con Inpactor3 sobre un genoma o FASTA de scaffolds.

Produce dos archivos con el MISMO formato que Inpactor2 para permitir
comparación directa:

  Inpactor3_predictions.tab    TSV de 8 columnas idéntico a Inpactor2:
      seqid  start  end  length  lineage  detect_prob  filter_prob  class_prob

  Inpactor3_library.fasta      FASTA con las secuencias detectadas.

Diferencias operacionales con Inpactor2:
  - Un solo forward por ventana (vs 3 redes en cascada).
  - `filter_prob` se rellena con "-" (no hay red Filter separada).
  - `detect_prob` es el objectness sigmoide, `class_prob` es softmax del linaje.

Uso:
    python scripts/predict_inpactor3.py \
        --genome test_genome.fasta \
        --checkpoint models/inpactor3_demo.pt \
        --output-dir results/inpactor3/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inpactor3.dataset import BASE2IDX, one_hot
from inpactor3.model import Yoro1D, decode
from inpactor3.nms import nms

WINDOW = 50_000
CELL = 100
STRIDE = 40_000  # solape 10 kb entre ventanas para no cortar elementos


def iter_fasta(path: Path):
    head, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(seq)
                head, seq = line[1:].split()[0], []
            else:
                seq.append(line)
        if head is not None:
            yield head, "".join(seq)


def scan_scaffold(model, seqid: str, seq: str, id2lin: dict, thr: float, device: str) -> list[tuple]:
    """Desliza ventanas de 50 kb sobre un scaffold y agrega detecciones."""
    all_boxes: list[tuple] = []
    L = len(seq)
    if L < WINDOW:
        # padding con N si es más corto
        seq = seq + "N" * (WINDOW - L)
        starts = [0]
    else:
        starts = list(range(0, L - WINDOW + 1, STRIDE))
        if starts[-1] + WINDOW < L:
            starts.append(L - WINDOW)

    for w_start in starts:
        window_seq = seq[w_start : w_start + WINDOW]
        x = torch.from_numpy(one_hot(window_seq)).unsqueeze(0).to(device)
        with torch.no_grad():
            pred = model(x)
        # objectness y clase también nos sirven aparte
        obj = torch.sigmoid(pred[0, :, 0]).cpu().numpy()
        cls_probs = torch.softmax(pred[0, :, 3:], dim=-1).cpu().numpy()
        cls_conf = cls_probs.max(axis=-1)

        boxes = decode(pred, CELL, thr=thr)[0]  # lista para batch=1
        for start, end, cls_id, score in boxes:
            abs_start = int(w_start + start)
            abs_end = int(w_start + end)
            if abs_end > L:
                abs_end = L
            cell = int(start // CELL)
            c_conf = float(cls_conf[cell]) if cell < len(cls_conf) else float(score)
            all_boxes.append(
                (seqid, abs_start, abs_end, id2lin.get(cls_id, "unknown"), score, c_conf, window_seq[int(start):int(end)])
            )

    # NMS sobre TODAS las cajas del scaffold (entre ventanas solapadas)
    if not all_boxes:
        return []
    # transformar a formato de nms 1D (start, end, cls, score) por seqid
    for_nms = [(b[1], b[2], b[3], b[4]) for b in all_boxes]
    kept = nms(for_nms, iou_thr=0.3)
    kept_set = set((k[0], k[1], k[2]) for k in kept)
    return [b for b in all_boxes if (b[1], b[2], b[3]) in kept_set]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True, type=Path, help="FASTA de scaffolds")
    ap.add_argument("--checkpoint", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--det-thr", type=float, default=0.5)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    ckpt = torch.load(args.checkpoint, map_location=args.device, weights_only=False)
    lin2id = ckpt["lin2id"]
    id2lin = {v: k for k, v in lin2id.items()}
    id2lin[0] = "background"
    n_classes = max(id2lin) + 1

    model = Yoro1D(n_classes=n_classes, cell_size=CELL).to(args.device)
    model.load_state_dict(ckpt["model"])
    model.eval()
    print(f"[model] {args.checkpoint} · {n_classes} clases · device={args.device}")

    tab = args.output_dir / "Inpactor3_predictions.tab"
    fa = args.output_dir / "Inpactor3_library.fasta"
    n_total = 0

    with open(tab, "w") as ftab, open(fa, "w") as ffa:
        for seqid, seq in iter_fasta(args.genome):
            print(f"[scan] {seqid} · {len(seq):,} bp")
            hits = scan_scaffold(model, seqid, seq, id2lin, args.det_thr, args.device)
            for seqid_out, start, end, lin, det_p, cls_p, subseq in hits:
                length = end - start
                # formato IDÉNTICO a Inpactor2_predictions.tab:
                # seqid  start  end  length  lineage  detect_prob  filter_prob  class_prob
                ftab.write(
                    f"{seqid_out}\t{start}\t{end}\t{length}\t{lin}\t"
                    f"{det_p:.4f}\t-\t{cls_p:.4f}\n"
                )
                ffa.write(f">{seqid_out}#{start}#{end}#{lin}#{det_p:.4f}\n")
                for k in range(0, len(subseq), 80):
                    ffa.write(subseq[k : k + 80] + "\n")
                n_total += 1

    print(f"\n[ok  ] {n_total} predicciones")
    print(f"[out ] {tab}")
    print(f"[out ] {fa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
