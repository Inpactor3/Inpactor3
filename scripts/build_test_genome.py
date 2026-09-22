"""
Construye un genoma sintético de PRUEBA usando ESPECIES HELD-OUT del corpus.

Divide el corpus fusionado por especie: 80% de especies van a entrenamiento
(vía scripts/train_demo.py con --exclude-species) y 20% se usan aquí para
construir el genoma de test. Así no hay data leakage por especie.

Salidas:
  test_genome/genome.fasta            scaffolds sintéticos con LTR-RTs plantados
  test_genome/ground_truth.tab        formato Inpactor2 (referencia)
  test_genome/train_species.txt       especies para entrenamiento
  test_genome/test_species.txt        especies reservadas para test

Uso:
    python scripts/build_test_genome.py \
        --corpus data/corpus/inpactor3_panteon.fasta \
        --output-dir test_genome/panteon \
        --n-scaffolds 5 --scaffold-len 200000 --seed 42
"""
from __future__ import annotations

import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from inpactor3.dataset import LTR, load_ltrs, random_dna  # noqa: E402


def load_by_species(fasta: Path) -> tuple[dict[str, list[LTR]], list[str]]:
    """Agrupa LTR-RTs por especie (leyendo header canónico con sp=)."""
    ltrs, lineages = load_ltrs(fasta, limit=None)
    by_sp: dict[str, list[LTR]] = defaultdict(list)
    # necesitamos la especie por LTR; releemos el header para extraerla
    with open(fasta) as f:
        idx = 0
        for line in f:
            if line.startswith(">"):
                sp = "unknown"
                for tag in line[1:].split("|"):
                    if tag.startswith("sp="):
                        sp = tag[3:].strip()
                        break
                if idx < len(ltrs):
                    by_sp[sp].append(ltrs[idx])
                    idx += 1
    return by_sp, lineages


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", type=Path,
                    default=ROOT / "data" / "corpus" / "inpactor3_v0.fasta")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "test_genome")
    ap.add_argument("--n-scaffolds", type=int, default=3)
    ap.add_argument("--scaffold-len", type=int, default=200_000)
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[read] {args.corpus}")
    by_sp, lineages = load_by_species(args.corpus)
    print(f"[ok  ] {sum(len(v) for v in by_sp.values())} LTR-RTs · "
          f"{len(by_sp)} especies · {len(lineages)} linajes")

    species = sorted(by_sp.keys())
    rng.shuffle(species)
    n_test = max(1, int(len(species) * args.test_frac))
    test_sp = set(species[:n_test])
    train_sp = set(species[n_test:])
    test_pool = [ltr for sp in test_sp for ltr in by_sp[sp]]
    print(f"[split] train={len(train_sp)} especies · "
          f"test={len(test_sp)} especies · {len(test_pool)} LTR-RTs candidatos")

    with open(args.output_dir / "train_species.txt", "w") as f:
        f.write("\n".join(sorted(train_sp)))
    with open(args.output_dir / "test_species.txt", "w") as f:
        f.write("\n".join(sorted(test_sp)))

    # ---- fabricar scaffolds usando SOLO especies held-out ----
    scaffolds: list[tuple[str, str]] = []
    truth: list[tuple[str, int, int, str]] = []
    for i in range(1, args.n_scaffolds + 1):
        name = f"chr{i}"
        parts: list[str] = []
        cursor = 0
        while cursor < args.scaffold_len - 30_000:
            e = rng.choice(test_pool)
            flank = rng.randint(5_000, 25_000)
            if cursor + flank + len(e.seq) > args.scaffold_len:
                break
            parts.append(random_dna(flank, rng))
            cursor += flank
            start = cursor
            parts.append(e.seq)
            cursor += len(e.seq)
            truth.append((name, start, cursor, e.lineage))
        parts.append(random_dna(args.scaffold_len - cursor, rng))
        scaffolds.append((name, "".join(parts)[: args.scaffold_len]))

    fa = args.output_dir / "genome.fasta"
    with open(fa, "w") as f:
        for name, seq in scaffolds:
            f.write(f">{name}\n")
            for k in range(0, len(seq), 80):
                f.write(seq[k : k + 80] + "\n")
    total_bp = sum(len(s) for _, s in scaffolds)
    print(f"[out ] {fa}  ({total_bp:,} bp · {len(scaffolds)} scaffolds)")

    gt = args.output_dir / "ground_truth.tab"
    with open(gt, "w") as f:
        for seqid, s, e, lin in truth:
            f.write(f"{seqid}\t{s}\t{e}\t{e-s}\t{lin}\t1.0000\t-\t1.0000\n")
    print(f"[out ] {gt}  ({len(truth)} LTR-RTs verdaderos)")

    print("\n[distribución de linajes plantados]")
    for lin, n in Counter(t[3] for t in truth).most_common():
        print(f"  {lin:12s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
