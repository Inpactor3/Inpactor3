"""
Fusiona FASTAs de fuentes heterogéneas (InpactorDB, PanTEon, RepetDB, Dfam
consenso LTR, TREP...) en un único corpus canónico para Inpactor3.

Pasos:
  1) Ingesta multi-fuente con parser específico por origen.
  2) Normalización taxonómica contra data/lineage_map.tsv (linaje canónico).
  3) Filtro de longitud [1kb, 25kb] y N-content < 5%.
  4) Deduplicación exacta por hash de secuencia.
  5) Emisión de FASTA unificado con header canónico:
        >{id}|src={fuente}|lin={CANONICAL_LINEAGE}|sp={species}|len={n}
     y JSONL de manifiesto con proveniencia por secuencia.

Dedup difusa (CD-HIT/MMseqs2) se hace en un paso posterior por eficiencia;
este script solo colapsa duplicados exactos (rápido, ya elimina ~30%).

Uso:
    python scripts/merge_corpus.py \
        --inpactordb data/raw/inpactordb_nr/InpactorDB_non_redundant_final_V5.fasta \
        --panteon data/raw/PanTEon_Database_v1.6.2.fasta \
        --lineage-map data/lineage_map.tsv \
        --out data/corpus/inpactor3_v0.fasta \
        --manifest data/corpus/inpactor3_v0.manifest.jsonl
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def load_lineage_map(path: Path) -> dict[str, str]:
    """Devuelve dict {alias_upper: canonical}. Cubre wicker+rexdb+dfam aliases."""
    m: dict[str, str] = {}
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            canon = row["canonical"].upper()
            m[canon] = canon
            # último segmento de rexdb_path
            m[row["rexdb_path"].split("/")[-1].upper()] = canon
            # último segmento de dfam_class
            m[row["dfam_class"].split(";")[-1].upper()] = canon
    return m


PANTEON_LTR_SUPERFAMS = {"COPIA", "GYPSY", "BELPAO", "ERV", "LARD", "TRIM"}


def canon_lineage(raw: str, m: dict[str, str]) -> str:
    """Resuelve linaje/superfamilia canónico.

    Para InpactorDB usa lineage_map.tsv (linajes tipo ALE, TAT...).
    Para PanTEon (superfamilias tipo COPIA, GYPSY...) las pasa directas
    si están en el set de superfamilias LTR válidas.
    """
    key = raw.strip().upper().replace(" ", "_")
    # PanTEon: superfamilias LTR ya son canónicas
    if key in PANTEON_LTR_SUPERFAMS:
        return key
    if key == "NONLTR":
        return "UNKNOWN"  # se descarta
    if key in m:
        return m[key]
    # compuestos frecuentes en InpactorDB V5
    if "/" in key:
        first = key.split("/")[0]
        if first in m:
            return m[first]
    aliases = {
        "COPIA/ALE": "ALE",
        "GYPSY/TEKAY": "TEKAY",
        "RETAND": "TEKAY",
        "CHROMOVIRIDAE": "CHROMOVIRUS",
        "TORK/TAR": "TORK",
        "ORYCO/IVANA": "IVANA",
        "ALE/RETROFIT": "ALE",
        "TAR/TORK": "TAR",
        "IVANA/ORYCO": "IVANA",
    }
    return aliases.get(key, "UNKNOWN")


def parse_inpactordb(h: str) -> tuple[str, str, str]:
    """
    V5 real: SUPERFAM-LINEAGE-Family-Species-Source-LENbp-ID (separado por '-').
    LINEAGE puede venir como compuesto (TORK/TAR); se resuelve en canon_lineage.
    """
    parts = h.split("-")
    lineage = parts[1] if len(parts) > 1 else "unknown"
    species = parts[3] if len(parts) > 3 else "unknown"
    return "inpactordb", lineage, species


PANTEON_LTR_SUPERFAMS = {"COPIA", "GYPSY", "BELPAO", "ERV", "LARD", "TRIM"}


def parse_panteon(h: str) -> tuple[str, str, str]:
    """
    Header PanTEon real: SEQID#CLASSI/LTR/COPIA @Species with spaces
    Devuelve (src, superfamily_como_etiqueta, species).
    Devuelve superfamily="NONLTR" para elementos que no sean LTR (para filtrar).
    """
    if "#" not in h:
        return "panteon", "unknown", "unknown"
    _, rest = h.split("#", 1)
    if "@" in rest:
        classif, species_raw = rest.split("@", 1)
        species = species_raw.strip().replace(" ", "_")
    else:
        classif, species = rest.strip(), "unknown"
    parts = classif.strip().split("/")
    order = parts[1] if len(parts) > 1 else ""
    superfam = parts[2] if len(parts) > 2 else "unknown"
    # Filtrar solo LTR-RTs; el resto se marca para descarte
    if order != "LTR":
        return "panteon", "NONLTR", species
    return "panteon", superfam, species


PARSERS = {"inpactordb": parse_inpactordb, "panteon": parse_panteon}


def iter_fasta(path: Path):
    head, seq = None, []
    with open(path, "rt", errors="replace") as f:
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inpactordb", type=Path)
    ap.add_argument("--panteon", type=Path)
    ap.add_argument("--lineage-map", type=Path, default=Path("data/lineage_map.tsv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--max-len", type=int, default=25000)
    ap.add_argument("--max-n-frac", type=float, default=0.05)
    args = ap.parse_args()

    lin_map = load_lineage_map(args.lineage_map)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    seen: set[str] = set()
    lineage_hist: Counter = Counter()
    src_hist: Counter = Counter()
    kept = dropped_len = dropped_n = dropped_dup = dropped_unknown = 0

    with open(args.out, "w") as out_fa, open(args.manifest, "w") as out_mf:
        sources = [
            ("inpactordb", args.inpactordb),
            ("panteon", args.panteon),
        ]
        for src, path in sources:
            if path is None or not path.exists():
                print(f"[skip] {src}: no encontrado")
                continue
            print(f"[read] {src}: {path}")
            parser = PARSERS[src]
            for i, (h, s) in enumerate(iter_fasta(path)):
                L = len(s)
                if L < args.min_len or L > args.max_len:
                    dropped_len += 1
                    continue
                n_frac = (s.upper().count("N")) / max(L, 1)
                if n_frac > args.max_n_frac:
                    dropped_n += 1
                    continue
                sha = hashlib.sha1(s.upper().encode()).hexdigest()
                if sha in seen:
                    dropped_dup += 1
                    continue
                seen.add(sha)
                _, raw_lin, species = parser(h)
                canon = canon_lineage(raw_lin, lin_map)
                if canon == "UNKNOWN":
                    dropped_unknown += 1
                    continue
                seq_id = f"{src}_{i:07d}"
                out_fa.write(
                    f">{seq_id}|src={src}|lin={canon}|sp={species}|len={L}\n"
                )
                for k in range(0, L, 80):
                    out_fa.write(s[k : k + 80] + "\n")
                out_mf.write(
                    json.dumps(
                        {
                            "id": seq_id,
                            "source": src,
                            "orig_header": h,
                            "raw_lineage": raw_lin,
                            "canonical_lineage": canon,
                            "species": species,
                            "length": L,
                            "sha1": sha,
                        }
                    )
                    + "\n"
                )
                lineage_hist[canon] += 1
                src_hist[src] += 1
                kept += 1

    print(f"\n[ok  ] kept {kept:,}")
    print(f"[drop] len={dropped_len:,}  N-rich={dropped_n:,}  dup={dropped_dup:,}  lin_unknown={dropped_unknown:,}")
    print(f"[src ] {dict(src_hist)}")
    print(f"[top linajes]")
    for k, v in lineage_hist.most_common():
        print(f"  {k:15s} {v:>7,}  ({v/max(kept,1)*100:5.2f}%)")
    print(f"[out ] {args.out}")
    print(f"[out ] {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
