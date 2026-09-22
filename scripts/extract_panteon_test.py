"""
Extrae un subset LTR de PanTEon como GENOMA de prueba para ambas herramientas.

Toma N LTR-RTs reales de PanTEon (uno por superfamilia balanceado) y los
planta en scaffolds de ~200 kb con flancos de ADN aleatorio, produciendo:

  test_panteon/genome.fasta                  scaffolds sintéticos con LTRs REALES
  test_panteon/ground_truth.tab              posiciones y superfamilia verdaderas
  test_panteon/subset_panteon_ltr.fasta      los LTRs originales sin cortar

Ambas herramientas (Inpactor2 e Inpactor3) reciben genome.fasta.
"""
from __future__ import annotations

import argparse
import random
from collections import defaultdict
from pathlib import Path

PANTEON = Path("data/raw/PanTEon_Database_v1.6.2.fasta")
LTR_SUPERFAMS = ["COPIA", "GYPSY", "ERV", "LARD", "BELPAO", "TRIM"]
BASES = "ACGT"


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


def parse_header(h: str) -> tuple[str, str]:
    """Devuelve (superfamily, species) del header PanTEon."""
    if "#" not in h:
        return "unknown", "unknown"
    _, rest = h.split("#", 1)
    if "@" in rest:
        classif, sp = rest.split("@", 1)
        species = sp.strip().replace(" ", "_")
    else:
        classif, species = rest.strip(), "unknown"
    parts = classif.strip().split("/")
    if len(parts) >= 3 and parts[1] == "LTR":
        return parts[2], species
    return "NONLTR", species


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-superfam", type=int, default=5,
                    help="Cuántos LTRs sacar de cada superfamilia (30 en total).")
    ap.add_argument("--scaffold-len", type=int, default=200_000)
    ap.add_argument("--n-scaffolds", type=int, default=3)
    ap.add_argument("--min-len", type=int, default=2000, help="filtro mínimo de longitud")
    ap.add_argument("--max-len", type=int, default=20000)
    ap.add_argument("--out-dir", type=Path, default=Path("test_panteon"))
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[read] escaneando {PANTEON}")
    by_sf: dict[str, list[tuple[str, str, str]]] = defaultdict(list)  # sf -> [(head, seq, species)]
    for h, seq in iter_fasta(PANTEON):
        sf, sp = parse_header(h)
        if sf not in LTR_SUPERFAMS:
            continue
        L = len(seq)
        if L < args.min_len or L > args.max_len:
            continue
        by_sf[sf].append((h, seq, sp))
        # sample fair: no acumular más del necesario
        if all(len(by_sf[s]) >= args.n_per_superfam * 3 for s in LTR_SUPERFAMS):
            break

    for sf in LTR_SUPERFAMS:
        print(f"  {sf:8s} disponibles: {len(by_sf[sf])}")

    # elegir n_per_superfam de cada
    chosen: list[tuple[str, str, str, str]] = []  # (sf, head, seq, species)
    for sf in LTR_SUPERFAMS:
        pool = by_sf[sf][:]
        rng.shuffle(pool)
        for h, seq, sp in pool[: args.n_per_superfam]:
            chosen.append((sf, h, seq, sp))
    rng.shuffle(chosen)
    print(f"\n[chosen] {len(chosen)} LTRs (balanceado por superfamilia)")

    # guardar subset FASTA sin flanco para referencia
    sub_fa = args.out_dir / "subset_panteon_ltr.fasta"
    with open(sub_fa, "w") as f:
        for sf, h, seq, sp in chosen:
            f.write(f">{h}\n")
            for k in range(0, len(seq), 80):
                f.write(seq[k : k + 80] + "\n")
    print(f"[out ] {sub_fa}")

    # construir scaffolds
    scaffolds: list[tuple[str, str]] = []
    truth: list[tuple[str, int, int, str, str]] = []  # (scaffold, start, end, sf, species)
    idx = 0
    for i in range(1, args.n_scaffolds + 1):
        name = f"scaffold_{i}"
        parts: list[str] = []
        cursor = 0
        while cursor < args.scaffold_len - 30_000 and idx < len(chosen):
            sf, h, seq, sp = chosen[idx]
            flank = rng.randint(5_000, 25_000)
            if cursor + flank + len(seq) > args.scaffold_len:
                break
            parts.append("".join(rng.choices(BASES, k=flank)))
            cursor += flank
            start = cursor
            parts.append(seq)
            cursor += len(seq)
            truth.append((name, start, cursor, sf, sp))
            idx += 1
        parts.append("".join(rng.choices(BASES, k=args.scaffold_len - cursor)))
        scaffolds.append((name, "".join(parts)[: args.scaffold_len]))

    fa = args.out_dir / "genome.fasta"
    with open(fa, "w") as f:
        for name, seq in scaffolds:
            f.write(f">{name}\n")
            for k in range(0, len(seq), 80):
                f.write(seq[k : k + 80] + "\n")
    total_bp = sum(len(s) for _, s in scaffolds)
    print(f"[out ] {fa}  ({total_bp:,} bp · {len(scaffolds)} scaffolds)")

    gt = args.out_dir / "ground_truth.tab"
    with open(gt, "w") as f:
        for seqid, s, e, sf, sp in truth:
            f.write(f"{seqid}\t{s}\t{e}\t{e-s}\t{sf}\t1.0000\t-\t1.0000\n")
    print(f"[out ] {gt}  ({len(truth)} LTR-RTs plantados)")

    print("\n[distribución plantada]")
    from collections import Counter
    for sf, n in Counter(t[3] for t in truth).most_common():
        print(f"  {sf:8s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
