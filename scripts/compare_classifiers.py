"""
Comparación directa Inpactor2 vs Inpactor3 sobre secuencias PanTEon LTR.

Carga cada clasificador con sus pesos preentrenados y predice, para cada
secuencia PanTEon del subset, la etiqueta. Compara contra la superfamilia
verdadera del header.

Inpactor2 pipeline:
    seq → one-hot (5, 50000) → kmer_extractor(2D CNN) → StandardScaler
       → PCA → Inpactor_Class.hdf5 → lineage (20 clases InpactorDB)

Inpactor3 pipeline:
    seq → one-hot (4, 50000) → Yoro1D → argmax por celda mayor
       → lineage canónico Wicker

Ambas etiquetas se agregan a superfamilia (COPIA / GYPSY / OTHER) para
comparar contra la etiqueta PanTEon.

Uso:
    python scripts/compare_classifiers.py \
        --fasta test_panteon/subset_panteon_ltr.fasta \
        --inpactor2-dir /home/cami/Desktop/ProyectoIntegrador/Inpactor2 \
        --inpactor3-ckpt models/inpactor3_demo.pt \
        --report test_panteon/classifier_comparison.md
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

# --- setup rutas ---
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


LANGU = ["A", "C", "G", "T", "N"]
TOTAL_WIN = 50_000

# Mapeo de lineage (Inpactor2 InpactorDB) → superfamilia
LINEAGE_TO_SF = {
    "Negative": "OTHER",
    "RLC/ALE/RETROFIT": "COPIA", "RLC/ANGELA": "COPIA", "RLC/BIANCA": "COPIA",
    "RLC/IKEROS": "COPIA", "RLC/IVANA/ORYCO": "COPIA", "RLC/TAR/TORK": "COPIA",
    "RLC/TORK/TAR": "COPIA", "RLC/SIRE": "COPIA",
    "RLG/CRM": "GYPSY", "RLG/GALADRIEL": "GYPSY", "RLG/REINA": "GYPSY",
    "RLG/TEKAY/DEL": "GYPSY", "RLG/ATHILA": "GYPSY", "RLG/TAT": "GYPSY",
}
# Mapeo de lineage Inpactor3 canónico → superfamilia
INP3_TO_SF = {
    "ALE": "COPIA", "ANGELA": "COPIA", "BIANCA": "COPIA", "IKEROS": "COPIA",
    "IVANA": "COPIA", "ORYCO": "COPIA", "OSSER": "COPIA", "RETROFIT": "COPIA",
    "SIRE": "COPIA", "TAR": "COPIA", "TORK": "COPIA",
    "TAT": "GYPSY", "DEL": "GYPSY", "TEKAY": "GYPSY", "ATHILA": "GYPSY",
    "CRM": "GYPSY", "GALADRIEL": "GYPSY", "REINA": "GYPSY",
    "background": "OTHER", "unknown": "OTHER",
}


def fasta2onehot5(seq: str) -> np.ndarray:
    """Codifica secuencia como (5, TOTAL_WIN) — formato Inpactor2."""
    arr = np.zeros((5, TOTAL_WIN), dtype=bool)
    for i, nt in enumerate(seq[:TOTAL_WIN]):
        idx = LANGU.index(nt.upper()) if nt.upper() in LANGU else 4
        arr[idx, i] = True
    return arr


def fasta2onehot4(seq: str) -> np.ndarray:
    """Codifica secuencia como (4, TOTAL_WIN) — formato Inpactor3."""
    arr = np.zeros((4, TOTAL_WIN), dtype=np.float32)
    for i, nt in enumerate(seq[:TOTAL_WIN]):
        j = "ACGT".find(nt.upper())
        if j >= 0:
            arr[j, i] = 1.0
    return arr


def parse_panteon_sf(h: str) -> tuple[str, str]:
    if "#" not in h:
        return "unknown", "unknown"
    _, rest = h.split("#", 1)
    if "@" in rest:
        cls, sp = rest.split("@", 1)
        species = sp.strip().replace(" ", "_")
    else:
        cls, species = rest.strip(), "unknown"
    parts = cls.strip().split("/")
    if len(parts) >= 3 and parts[1] == "LTR":
        return parts[2], species
    return "NONLTR", species


def iter_fasta(path: Path):
    head, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(seq)
                head, seq = line[1:], []
            else:
                seq.append(line)
        if head is not None:
            yield head, "".join(seq)


def _build_kmer_extractor(inp2_dir: Path, input_shape):
    """Reproduce Inpactor2_utils.kmer_extractor_model compatible con Keras 3.x.

    Usa layer.set_weights(...) en vez del argumento `weights=` obsoleto.
    """
    import tensorflow as tf
    weights = np.load(inp2_dir / "Models" / "Weights_SL.npy", allow_pickle=True)
    from tensorflow.keras import layers as L
    inputs = tf.keras.Input(shape=input_shape, name="input_1")
    layer_outs = []
    sum_axis2 = L.Lambda(lambda t: tf.reduce_sum(t, axis=-2))
    for k in range(1, 7):
        W = weights[(k - 1) * 2]
        b = weights[(k - 1) * 2 + 1]
        n_filters = [4, 16, 64, 256, 1024, 4096][k - 1]
        conv = L.Conv2D(
            n_filters, (5, k), strides=(1, 1), activation="relu",
            use_bias=True, name=f"k_{k}",
        )
        y = conv(inputs)
        conv.set_weights([W, b])
        y = sum_axis2(y)
        layer_outs.append(y)
    concat = L.Concatenate(axis=2)(layer_outs)
    outputs = L.Flatten()(concat)
    model = tf.keras.Model(inputs=inputs, outputs=outputs)
    for layer in model.layers:
        layer.trainable = False
    return model


def run_inpactor2(seqs: list[np.ndarray], inp2_dir: Path) -> tuple[list[str], list[float]]:
    """Corre el clasificador de Inpactor2 sobre secuencias one-hot (5, W)."""
    sys.path.insert(0, str(inp2_dir))
    # Forzar el backend legacy (Keras 2) para cargar el .hdf5 antiguo
    os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")
    import tensorflow as tf
    try:
        import tf_keras
        tf.keras = tf_keras
    except ImportError:
        pass
    from joblib import load

    def f1_m(y_true, y_pred):  # placeholder para load_model
        return tf.constant(0.0)

    batch = np.stack(seqs).astype(np.float32)  # (N, 5, W)
    # Keras 3 espera channels_last o input con 4 dims para Conv2D
    batch = batch[..., np.newaxis]  # (N, 5, W, 1)
    print(f"[inp2] batch shape: {batch.shape}")

    # kmer extraction
    print("[inp2] extrayendo k-mers (esto tarda un poco por secuencia)...")
    kext = _build_kmer_extractor(inp2_dir, input_shape=batch.shape[1:])
    kmer_counts = kext.predict(batch, batch_size=1, verbose=0)

    scaler = load(inp2_dir / "Models" / "std_scaler.bin")
    features = scaler.transform(kmer_counts)
    pca = load(inp2_dir / "Models" / "std_pca.bin")
    features_pca = pca.transform(features)

    print("[inp2] clasificando...")
    model = tf.keras.models.load_model(
        inp2_dir / "Models" / "Inpactor_Class.hdf5",
        custom_objects={"f1_m": f1_m},
    )
    probs = model.predict(features_pca, verbose=0)
    ids = probs.argmax(axis=-1)

    lineages_names_dic = {
        0: "Negative", 1: "RLC/ALE/RETROFIT", 3: "RLC/ANGELA", 4: "RLC/BIANCA",
        8: "RLC/IKEROS", 9: "RLC/IVANA/ORYCO", 11: "RLC/TAR/TORK", 12: "RLC/TORK/TAR",
        13: "RLC/SIRE", 14: "RLG/CRM", 16: "RLG/GALADRIEL", 17: "RLG/REINA",
        18: "RLG/TEKAY/DEL", 19: "RLG/ATHILA", 20: "RLG/TAT",
    }
    labels = [lineages_names_dic.get(int(i), f"class_{int(i)}") for i in ids]
    scores = [float(probs[i, ids[i]]) for i in range(len(ids))]
    return labels, scores


def run_inpactor3(seqs: list[np.ndarray], ckpt_path: Path) -> tuple[list[str], list[float]]:
    """Corre Yoro1D sobre secuencias one-hot (4, W)."""
    import torch
    from inpactor3.model import Yoro1D

    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    lin2id = ckpt["lin2id"]
    id2lin = {v: k for k, v in lin2id.items()}
    id2lin[0] = "background"
    n_classes = max(id2lin) + 1
    model = Yoro1D(n_classes=n_classes, cell_size=100)
    model.load_state_dict(ckpt["model"])
    model.eval()

    labels: list[str] = []
    scores: list[float] = []
    print(f"[inp3] procesando {len(seqs)} secuencias...")
    with torch.no_grad():
        for arr in seqs:
            x = torch.from_numpy(arr).unsqueeze(0)  # (1, 4, W)
            pred = model(x)
            obj = torch.sigmoid(pred[0, :, 0])
            cls_probs = torch.softmax(pred[0, :, 3:], dim=-1)
            # Escoge celda de mayor objectness × prob de clase
            score_cell = obj.unsqueeze(-1) * cls_probs
            best_cell = score_cell.max(dim=-1).values.argmax().item()
            best_cls = cls_probs[best_cell].argmax().item()
            labels.append(id2lin.get(best_cls, "unknown"))
            scores.append(float(score_cell[best_cell, best_cls]))
    return labels, scores


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--inpactor2-dir", required=True, type=Path)
    ap.add_argument("--inpactor3-ckpt", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    args = ap.parse_args()

    # ---- leer secuencias ----
    entries: list[tuple[str, str, str, str]] = []  # (id, seq, gt_sf, species)
    for h, seq in iter_fasta(args.fasta):
        sf, sp = parse_panteon_sf(h)
        if sf == "NONLTR":
            continue
        seq_pad = (seq + "N" * TOTAL_WIN)[:TOTAL_WIN]
        entries.append((h.split()[0], seq_pad, sf, sp))
    print(f"[read] {len(entries)} secuencias LTR de PanTEon")

    if not entries:
        print("[err] sin secuencias LTR válidas")
        return 1

    # ---- Inpactor2 ----
    seqs5 = [fasta2onehot5(e[1]) for e in entries]
    inp2_lin, inp2_score = run_inpactor2(seqs5, args.inpactor2_dir)

    # ---- Inpactor3 ----
    seqs4 = [fasta2onehot4(e[1]) for e in entries]
    inp3_lin, inp3_score = run_inpactor3(seqs4, args.inpactor3_ckpt)

    # ---- comparar ----
    rows = []
    for (id_, _, gt, sp), l2, s2, l3, s3 in zip(entries, inp2_lin, inp2_score, inp3_lin, inp3_score):
        sf2 = LINEAGE_TO_SF.get(l2, "OTHER")
        sf3 = INP3_TO_SF.get(l3, "OTHER")
        rows.append({
            "id": id_, "species": sp, "gt_sf": gt,
            "inp2_lin": l2, "inp2_sf": sf2, "inp2_score": s2,
            "inp3_lin": l3, "inp3_sf": sf3, "inp3_score": s3,
        })

    # ---- métricas ----
    acc2 = sum(1 for r in rows if r["gt_sf"] == r["inp2_sf"]) / len(rows)
    acc3 = sum(1 for r in rows if r["gt_sf"] == r["inp3_sf"]) / len(rows)
    agree = sum(1 for r in rows if r["inp2_sf"] == r["inp3_sf"]) / len(rows)

    conf_i2 = defaultdict(lambda: defaultdict(int))
    conf_i3 = defaultdict(lambda: defaultdict(int))
    for r in rows:
        conf_i2[r["gt_sf"]][r["inp2_sf"]] += 1
        conf_i3[r["gt_sf"]][r["inp3_sf"]] += 1

    # ---- reporte ----
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w") as f:
        f.write("# Comparación Inpactor2 vs Inpactor3 sobre PanTEon LTR\n\n")
        f.write(f"- Secuencias evaluadas: **{len(rows)}**\n")
        f.write(f"- Accuracy Inpactor2 (superfamilia): **{acc2:.3f}**\n")
        f.write(f"- Accuracy Inpactor3 (superfamilia): **{acc3:.3f}**\n")
        f.write(f"- Concordancia entre ambos: **{agree:.3f}**\n\n")

        f.write("## Tabla completa\n\n")
        f.write("| seqid | especie | verdad | Inpactor2 lin | Inpactor2 SF | Inpactor3 lin | Inpactor3 SF |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            match2 = "✓" if r["gt_sf"] == r["inp2_sf"] else "✗"
            match3 = "✓" if r["gt_sf"] == r["inp3_sf"] else "✗"
            f.write(
                f"| {r['id'][:20]} | {r['species'][:20]} | **{r['gt_sf']}** | "
                f"{r['inp2_lin']} | {r['inp2_sf']} {match2} | "
                f"{r['inp3_lin']} | {r['inp3_sf']} {match3} |\n"
            )

        f.write("\n## Confusión Inpactor2 (filas=verdad, columnas=predicho)\n\n")
        cols = sorted(set(k for c in conf_i2.values() for k in c) | set(conf_i2.keys()))
        f.write("| ↓verdad \\ pred→ | " + " | ".join(cols) + " |\n")
        f.write("|" + "---|" * (len(cols) + 1) + "\n")
        for r in sorted(conf_i2.keys()):
            row = [str(conf_i2[r][c]) for c in cols]
            f.write(f"| **{r}** | " + " | ".join(row) + " |\n")

        f.write("\n## Confusión Inpactor3 (filas=verdad, columnas=predicho)\n\n")
        cols = sorted(set(k for c in conf_i3.values() for k in c) | set(conf_i3.keys()))
        f.write("| ↓verdad \\ pred→ | " + " | ".join(cols) + " |\n")
        f.write("|" + "---|" * (len(cols) + 1) + "\n")
        for r in sorted(conf_i3.keys()):
            row = [str(conf_i3[r][c]) for c in cols]
            f.write(f"| **{r}** | " + " | ".join(row) + " |\n")

    print(f"\n[out] {args.report}")
    print(f"Accuracy Inpactor2: {acc2:.3f}")
    print(f"Accuracy Inpactor3: {acc3:.3f}")
    print(f"Concordancia:       {agree:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
