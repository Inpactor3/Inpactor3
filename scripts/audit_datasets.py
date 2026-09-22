"""
Descarga y auditoría de bases de datos candidatas para el corpus de Inpactor3.

URLs verificadas 2026-09-20 contra la API de Zenodo y las páginas oficiales:
  - InpactorDB v5 (Zenodo 6380332, DOI 10.5281/zenodo.6380332)  [CC-BY]
  - InpactorDB negatives (Zenodo 4543905)                        [CC-BY]
  - PanTEon v1.6.2 (Zenodo 21372179, DOI 10.5281/zenodo.21372179)[CC-BY]
  - Dfam 40.0 curated (dfam.org, release 2026-05-29)              [CC0]

Bases no descargadas por script:
  - REXDB (Bitbucket petrnovak/re_databases) → clon manual o vía TEsorter
  - GyDB (gydb.org/collection HMM) → navegación web requerida
  - RepetDB (INRAE URGI) → protegida por anti-scraping Anubis
  - RepBase (GIRI) → licencia de pago

Uso:
    python scripts/audit_datasets.py --data-dir data/raw --which inpactordb
    python scripts/audit_datasets.py --data-dir data/raw --which all
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

DATASETS = {
    "inpactordb_nr": {
        "url": "https://zenodo.org/api/records/6380332/files/InpactorDB_non_redundant_final_V5.zip/content",
        "filename": "InpactorDB_non_redundant_final_V5.zip",
        "size_mb": 154.9,
        "unzip_to": "inpactordb_nr",
        "license": "CC-BY-4.0",
    },
    "inpactordb_full": {
        "url": "https://zenodo.org/api/records/6380332/files/InpactorDB_redundant_final_V5.zip/content",
        "filename": "InpactorDB_redundant_final_V5.zip",
        "size_mb": 292.5,
        "unzip_to": "inpactordb_full",
        "license": "CC-BY-4.0",
    },
    "inpactordb_negatives": {
        "url": "https://zenodo.org/api/records/4543905/files/negative_instances_raw.zip/content",
        "filename": "negative_instances_raw.zip",
        "size_mb": 865.0,
        "unzip_to": "inpactordb_negatives",
        "license": "CC-BY-4.0",
    },
    "panteon_fasta": {
        "url": "https://zenodo.org/api/records/21372179/files/PanTEon_Database_v1.6.2.fasta/content",
        "filename": "PanTEon_Database_v1.6.2.fasta",
        "size_mb": 992.2,
        "license": "CC-BY-4.0",
    },
    "panteon_metadata": {
        "url": "https://zenodo.org/api/records/21372179/files/PanTEon_Database_metadata_v1.6.2.csv/content",
        "filename": "PanTEon_Database_metadata_v1.6.2.csv",
        "size_mb": 0.4,
        "license": "CC-BY-4.0",
    },
    "dfam_curated_consensus": {
        "url": "https://www.dfam.org/releases/current/families/FamDB/dfam40.curated.consensus.0.h5.gz",
        "filename": "dfam40.curated.consensus.0.h5.gz",
        "size_mb": 28.3,
        "license": "CC0",
    },
}


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} ya existe ({human(dest.stat().st_size)})")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[get ] {url}")
    req = Request(url, headers={"User-Agent": "Inpactor3-audit/0.1"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urlopen(req) as r, open(tmp, "wb") as f:
        total = int(r.headers.get("Content-Length", 0))
        got = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            got += len(chunk)
            if total:
                print(f"\r       {human(got)} / {human(total)}", end="", flush=True)
    print()
    tmp.rename(dest)
    print(f"[done] {dest.name} ({human(dest.stat().st_size)})")


def maybe_unzip(zip_path: Path, out_dir: Path) -> None:
    if not zip_path.exists() or zip_path.suffix != ".zip":
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    if any(out_dir.iterdir()):
        print(f"[skip] {out_dir} ya poblado")
        return
    print(f"[unzip] {zip_path.name} -> {out_dir}")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)


def iter_fasta_headers(path: Path):
    with open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


def parse_inpactordb(h: str) -> tuple[str, str, str]:
    """
    Header InpactorDB V5 real (separado por '-'):
      SUPERFAM-LINEAGE-Family-Species-Source-LENGTHbp-ID
    Ej.: >RLC-TORK/TAR-Brassicaceae-Arabidopsis_thaliana-Repbase-5186bp-1

    LINEAGE puede ser compuesto (TORK/TAR, ORYCO/IVANA, ALE/RETROFIT):
    lo dejamos tal cual; la normalización canónica se hace en merge_corpus.
    """
    parts = h.split("-")
    superfam = parts[0] if len(parts) > 0 else "unknown"
    lineage = parts[1] if len(parts) > 1 else "unknown"
    species = parts[3] if len(parts) > 3 else "unknown"
    return superfam, lineage, species


def parse_panteon(h: str) -> tuple[str, str, str]:
    """
    Header PanTEon real:
        SEQID#CLASSI/LTR/COPIA @Species with spaces
        SEQID#CLASSII/HELITRON/HELITRON @Species
    Devuelve (order, superfamily, species). Ejemplo:
        parse_panteon("PDB01#CLASSI/LTR/GYPSY @Puma concolor")
          → ("LTR", "GYPSY", "Puma_concolor")
    """
    if "#" not in h:
        return "unknown", "unknown", "unknown"
    _, rest = h.split("#", 1)
    if "@" in rest:
        classif, species_raw = rest.split("@", 1)
        species = species_raw.strip().replace(" ", "_")
    else:
        classif, species = rest.strip(), "unknown"
    parts = classif.strip().split("/")
    order = parts[1] if len(parts) > 1 else "unknown"
    superfam = parts[2] if len(parts) > 2 else "unknown"
    return order, superfam, species


def audit_fasta(fasta: Path, parser) -> dict:
    superfams: Counter = Counter()
    lineages: Counter = Counter()
    species: Counter = Counter()
    combo: defaultdict[tuple[str, str], int] = defaultdict(int)
    total = 0
    for h in iter_fasta_headers(fasta):
        sf, lin, sp = parser(h)
        superfams[sf] += 1
        lineages[lin] += 1
        species[sp] += 1
        combo[(sp, lin)] += 1
        total += 1
    return {
        "file": fasta.name,
        "total": total,
        "superfams": superfams,
        "lineages": lineages,
        "species": species,
        "combo": combo,
    }


def print_report(res: dict, top: int = 15) -> None:
    print(f"\n=== {res['file']} — {res['total']:,} secuencias ===")
    print("\n[Superfamilias / reinos]")
    for k, v in res["superfams"].most_common(top):
        print(f"  {k:25s} {v:>8,}  ({v/max(res['total'],1)*100:5.2f}%)")
    print(f"\n[Top {top} linajes]")
    for k, v in res["lineages"].most_common(top):
        print(f"  {k:25s} {v:>8,}  ({v/max(res['total'],1)*100:5.2f}%)")
    print(f"\n[Top {top} especies]")
    for k, v in res["species"].most_common(top):
        print(f"  {k:35s} {v:>8,}  ({v/max(res['total'],1)*100:5.2f}%)")
    n_sp = len(res["species"])
    n_lin = len(res["lineages"])
    counts = sorted(res["lineages"].values(), reverse=True)
    top10 = sum(counts[:10]) / max(res["total"], 1)
    print(f"\n[Cobertura] {n_sp} especies · {n_lin} linajes distintos")
    print(f"[Desbalance] top-10 linajes concentran {top10*100:.1f}% del corpus")


def write_tsv(res: dict, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("species\tlineage\tcount\n")
        for (sp, lin), n in sorted(res["combo"].items(), key=lambda x: -x[1]):
            f.write(f"{sp}\t{lin}\t{n}\n")
    print(f"[tsv ] {out}")


def process_dataset(name: str, data_dir: Path) -> None:
    if name not in DATASETS:
        print(f"[warn] dataset desconocido: {name}")
        return
    d = DATASETS[name]
    dest = data_dir / d["filename"]
    try:
        download(d["url"], dest)
    except Exception as e:
        print(f"[error] {name}: {e}", file=sys.stderr)
        return
    if "unzip_to" in d:
        maybe_unzip(dest, data_dir / d["unzip_to"])
        # localiza FASTAs dentro
        fastas = list((data_dir / d["unzip_to"]).rglob("*.fasta")) + list(
            (data_dir / d["unzip_to"]).rglob("*.fa")
        )
    elif dest.suffix in (".fasta", ".fa"):
        fastas = [dest]
    else:
        fastas = []

    parser = parse_panteon if name.startswith("panteon") else parse_inpactordb
    for fa in fastas:
        if fa.stat().st_size < 1000:
            continue
        res = audit_fasta(fa, parser)
        print_report(res)
        write_tsv(res, data_dir.parent / f"audit_{fa.stem}.tsv")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/raw")
    ap.add_argument(
        "--which",
        default="inpactordb_nr",
        help="Nombre de dataset o 'all' o coma-separado (inpactordb_nr,panteon_fasta,...)",
    )
    args = ap.parse_args()
    data_dir = Path(args.data_dir)

    if args.which == "all":
        names = list(DATASETS.keys())
    elif args.which == "list":
        for k, v in DATASETS.items():
            print(f"  {k:25s}  ~{v['size_mb']:>7.1f} MB  {v['license']}")
        return 0
    else:
        names = [x.strip() for x in args.which.split(",")]

    for n in names:
        process_dataset(n, data_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
