"""
Demo de entrenamiento Inpactor3 (YORO 1D) usando InpactorDB real.

Requiere haber ejecutado antes:
    python scripts/audit_datasets.py --data-dir data/raw

Uso:
    python scripts/train_demo.py --fasta data/raw/InpactorDB_non_redundant.fasta \
        --epochs 3 --batch 4 --limit-ltrs 2000 --windows 256

Muestra loss por época y decodifica una ventana de validación con NMS,
comparando cajas predichas vs verdaderas.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch
from torch.utils.data import DataLoader

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inpactor3.dataset import LTRWindowDataset, collate, load_ltrs
from inpactor3.model import Yoro1D, decode, detection_loss
from inpactor3.nms import nms


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--windows", type=int, default=256)
    ap.add_argument("--val-windows", type=int, default=32)
    ap.add_argument("--window-size", type=int, default=50_000)
    ap.add_argument("--cell-size", type=int, default=100)
    ap.add_argument("--max-per-window", type=int, default=3)
    ap.add_argument("--limit-ltrs", type=int, default=2000)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--det-thr", type=float, default=0.3, help="Umbral de objectness en decode")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    print(f"[data] cargando LTR-RTs de {args.fasta} (limit={args.limit_ltrs})")
    ltrs, lineages = load_ltrs(args.fasta, limit=args.limit_ltrs)
    print(f"[data] {len(ltrs)} LTR-RTs · {len(lineages)} linajes: {lineages[:8]}{'...' if len(lineages)>8 else ''}")
    if not ltrs:
        print("[error] FASTA vacío o no filtró nada; ejecuta primero scripts/audit_datasets.py")
        return 1

    common = dict(
        ltrs=ltrs,
        lineages=lineages,
        window_size=args.window_size,
        cell_size=args.cell_size,
        max_per_window=args.max_per_window,
    )
    train_ds = LTRWindowDataset(n_windows=args.windows, seed=1, **common)
    val_ds = LTRWindowDataset(n_windows=args.val_windows, seed=999, **common)
    dl = DataLoader(train_ds, batch_size=args.batch, shuffle=True, collate_fn=collate)
    vdl = DataLoader(val_ds, batch_size=args.batch, shuffle=False, collate_fn=collate)

    model = Yoro1D(n_classes=train_ds.n_classes, cell_size=args.cell_size).to(args.device)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)
    print(f"[model] Yoro1D · {sum(p.numel() for p in model.parameters())/1e6:.2f} M params · device={args.device}")

    for ep in range(1, args.epochs + 1):
        model.train()
        tot = 0.0
        parts = {"obj": 0.0, "off": 0.0, "len": 0.0, "cls": 0.0}
        n = 0
        for batch in dl:
            x = batch["x"].to(args.device)
            tgt = {k: batch[k].to(args.device) for k in ("y_obj", "y_off", "y_len", "y_cls")}
            pred = model(x)
            loss, comp = detection_loss(pred, tgt)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item()
            for k, v in comp.items():
                parts[k] += v
            n += 1
        avg = {k: v / n for k, v in parts.items()}
        print(
            f"[ep {ep:02d}] loss={tot/n:.4f}  obj={avg['obj']:.3f} off={avg['off']:.3f} "
            f"len={avg['len']:.3f} cls={avg['cls']:.3f}"
        )

    # Demostración de inferencia en 1 ventana
    model.eval()
    id2lin = {v: k for k, v in train_ds.lin2id.items()}
    id2lin[0] = "background"
    with torch.no_grad():
        batch = next(iter(vdl))
        x = batch["x"].to(args.device)
        pred = model(x)
        raw = decode(pred, args.cell_size, thr=args.det_thr)

    for i, (boxes_pred, boxes_true) in enumerate(zip(raw, batch["boxes"])):
        kept = nms(boxes_pred, iou_thr=0.3)
        print(f"\n--- ventana {i} ---")
        print(f"  verdaderas ({len(boxes_true)}):")
        for s, e, c in boxes_true:
            print(f"    [{s:>6d}, {e:>6d})  cls={id2lin.get(c, c)}")
        print(f"  predichas post-NMS ({len(kept)}):")
        for s, e, c, sc in kept[:10]:
            print(f"    [{s:>8.0f}, {e:>8.0f})  cls={id2lin.get(c, c):15s}  score={sc:.2f}")
        if i >= 1:
            break

    out = Path("models/inpactor3_demo.pt")
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"model": model.state_dict(), "lin2id": train_ds.lin2id}, out)
    print(f"\n[save] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
