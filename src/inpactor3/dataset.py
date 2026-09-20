"""
Dataset YORO-style para Inpactor3.

Toma un FASTA real de InpactorDB (o PanTEon, GyDB, etc.) con LTR-RTs
etiquetados por linaje y construye ventanas genómicas sintéticas de tamaño
fijo (por defecto 50 000 bp), plantando 0-K elementos en posiciones
aleatorias sobre "flancos" de DNA aleatorio.

Cada ventana produce:
  - x: tensor one-hot (4, W) de la ventana
  - y_obj:   (S,)    presencia de objeto por celda-ancla
  - y_off:   (S,)    offset del inicio dentro de la celda, [0, 1)
  - y_len:   (S,)    log(longitud / cell_size)
  - y_cls:   (S,)    id de linaje (0 = background, 1..C = linajes)
  - boxes:   lista de (start, end, cls) verdadera (para eval)

Convención de anclas: una celda de tamaño `cell_size` (100 bp como en YORO);
el objeto se asigna a la celda que contiene su punto medio.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

BASES = "ACGT"
BASE2IDX = {b: i for i, b in enumerate(BASES)}


def one_hot(seq: str) -> np.ndarray:
    """(4, L) float32. N y otros → todo cero."""
    L = len(seq)
    x = np.zeros((4, L), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        j = BASE2IDX.get(b)
        if j is not None:
            x[j, i] = 1.0
    return x


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choices(BASES, k=n))


def parse_header(h: str) -> str:
    """
    Soporta dos formatos:
      1) InpactorDB V5 real: `SUPERFAM-LINEAGE-Family-Species-Source-LENbp-ID`
         → devuelve LINEAGE (posiblemente compuesto TORK/TAR → TORK).
      2) Corpus canónico Inpactor3: `id|src=..|lin=CANON|sp=..|len=..`
         → devuelve CANON directamente.
    """
    if "|lin=" in h:
        for tag in h.split("|"):
            if tag.startswith("lin="):
                return tag[4:]
    parts = h.split("-")
    if len(parts) >= 2:
        lin = parts[1]
        return lin.split("/")[0]
    return "unknown"


@dataclass
class LTR:
    seq: str
    lineage: str


def load_ltrs(
    fasta: Path, min_len: int = 1000, max_len: int = 20000, limit: int | None = None
) -> tuple[list[LTR], list[str]]:
    """Carga LTR-RTs de un FASTA (InpactorDB-style) filtrando por longitud."""
    ltrs: list[LTR] = []
    cur_head = None
    cur_seq: list[str] = []

    def flush():
        if cur_head is None:
            return
        seq = "".join(cur_seq)
        if min_len <= len(seq) <= max_len:
            ltrs.append(LTR(seq=seq, lineage=parse_header(cur_head)))

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                if limit and len(ltrs) >= limit:
                    break
                cur_head = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        flush()

    lineages = sorted({l.lineage for l in ltrs})
    return ltrs, lineages


class LTRWindowDataset(Dataset):
    def __init__(
        self,
        ltrs: list[LTR],
        lineages: list[str],
        window_size: int = 50_000,
        cell_size: int = 100,
        max_per_window: int = 3,
        n_windows: int = 512,
        seed: int = 0,
    ):
        self.ltrs = ltrs
        self.window_size = window_size
        self.cell_size = cell_size
        self.n_cells = window_size // cell_size
        self.max_per_window = max_per_window
        self.n_windows = n_windows
        # class 0 = background
        self.lin2id = {l: i + 1 for i, l in enumerate(lineages)}
        self.n_classes = len(lineages) + 1
        self.seed = seed

    def __len__(self) -> int:
        return self.n_windows

    def _sample_window(self, rng: random.Random):
        k = rng.randint(1, self.max_per_window)
        # elige LTR-RTs que quepan
        chosen: list[LTR] = []
        for _ in range(k * 4):
            e = rng.choice(self.ltrs)
            if len(e.seq) < self.window_size - 2 * self.cell_size:
                chosen.append(e)
            if len(chosen) == k:
                break
        # coloca sin solaparse
        placements: list[tuple[int, int, LTR]] = []
        for e in chosen:
            for _ in range(20):
                start = rng.randint(0, self.window_size - len(e.seq) - 1)
                end = start + len(e.seq)
                if all(end <= s or start >= t for s, t, _ in placements):
                    placements.append((start, end, e))
                    break
        placements.sort()

        # construye secuencia: flancos aleatorios entre LTR-RTs
        parts: list[str] = []
        cursor = 0
        for s, t, e in placements:
            parts.append(random_dna(s - cursor, rng))
            parts.append(e.seq)
            cursor = t
        parts.append(random_dna(self.window_size - cursor, rng))
        seq = "".join(parts)[: self.window_size]
        return seq, placements

    def __getitem__(self, idx: int):
        rng = random.Random(self.seed * 1_000_003 + idx)
        seq, placements = self._sample_window(rng)
        x = torch.from_numpy(one_hot(seq))  # (4, W)

        S = self.n_cells
        y_obj = torch.zeros(S, dtype=torch.float32)
        y_off = torch.zeros(S, dtype=torch.float32)
        y_len = torch.zeros(S, dtype=torch.float32)
        y_cls = torch.zeros(S, dtype=torch.long)  # 0 = bg
        boxes = []
        for s, t, e in placements:
            midpoint = (s + t) // 2
            cell = midpoint // self.cell_size
            if cell >= S:
                continue
            cell_start = cell * self.cell_size
            y_obj[cell] = 1.0
            y_off[cell] = (s - cell_start) / self.cell_size  # puede ser negativo si s<cell_start
            y_len[cell] = float(np.log(max((t - s), 1) / self.cell_size))
            y_cls[cell] = self.lin2id[e.lineage]
            boxes.append((s, t, self.lin2id[e.lineage]))

        return {
            "x": x,
            "y_obj": y_obj,
            "y_off": y_off,
            "y_len": y_len,
            "y_cls": y_cls,
            "boxes": boxes,
        }


def collate(batch):
    out = {
        "x": torch.stack([b["x"] for b in batch]),
        "y_obj": torch.stack([b["y_obj"] for b in batch]),
        "y_off": torch.stack([b["y_off"] for b in batch]),
        "y_len": torch.stack([b["y_len"] for b in batch]),
        "y_cls": torch.stack([b["y_cls"] for b in batch]),
        "boxes": [b["boxes"] for b in batch],
    }
    return out
