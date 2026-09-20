# Bases de datos candidatas para Inpactor3

Revisión de fuentes de LTR-retrotransposones (y TEs en general) evaluadas como
insumo del corpus de entrenamiento. Última verificación de URLs: **2026-09-20**.

## Fichas

### InpactorDB v5 (2022)
- Paper base: Orozco-Arias et al., *Genes* 12(2):190, 2021. DOI [10.3390/genes12020190](https://doi.org/10.3390/genes12020190).
- Zenodo: [record 6380332](https://zenodo.org/records/6380332), DOI 10.5281/zenodo.6380332.
- Archivos:
  - `InpactorDB_non_redundant_final_V5.zip` — 155 MB (FASTA no redundante)
  - `InpactorDB_redundant_final_V5.zip` — 292 MB (FASTA completo)
- Cobertura: Viridiplantae, ~195 especies.
- Header FASTA: `>SUPERFAM#LINEAGE#Family#Species#source#length#id`
- Clasificación: linajes tipo Wicker/REXDB fusionados (9 Copia + 7 Gypsy usados en Inpactor v1/v2).
- Licencia: **CC-BY 4.0**.
- **Rol**: núcleo del corpus.

### InpactorDB negatives
- Zenodo: [record 4543905](https://zenodo.org/records/4543905).
- `negative_instances_raw.zip` — 865 MB (secuencias genómicas no-LTR de 50 kb).
- Licencia CC-BY 4.0.
- **Rol**: ventanas negativas para la cabeza objectness del detector.

### PanTEon DB v1.6.2 (2025)
- Repo: [github.com/simonorozcoarias/PanTEon](https://github.com/simonorozcoarias/PanTEon).
- Zenodo: [record 21372179](https://zenodo.org/records/21372179), DOI 10.5281/zenodo.21372179.
- Archivos clave:
  - `PanTEon_Database_v1.6.2.fasta` — 992 MB
  - `PanTEon_Database_metadata_v1.6.2.csv` — 371 kB
  - `PanTEon_Database_v1.6.2_benchmark_edition.fasta` — 974 MB
  - Contenedores Singularity (CPU/GPU) y modelos entrenados por reino.
- ~240 000 TEs curados automáticamente (Dfam 3.9 + APTEdb + Ensembl 2025 + RepeatModeler2).
- Cobertura: Animalia + Plantae + Fungi.
- Licencia: **CC-BY 4.0**.
- **Rol**: ampliar taxonomía más allá de plantas; requiere dedup contra InpactorDB.

### GyDB 2.0
- Web: [gydb.org](https://gydb.org). Paper: Llorens et al., *NAR* 39:D70, 2011.
- Formato: HMMs (dominios GAG, PR, RT, RH, INT, ENV) + FASTA de proteínas + alineamientos.
- Cobertura: LTR-RTs y retrovirus, todos los reinos.
- Licencia: académica, cita obligatoria.
- **Rol**: **verificación de dominios**, no fuente de secuencia LTR-RT completa.

### REXDB (PGSB)
- Paper: Neumann et al., *Mobile DNA* 10:1, 2019. DOI [10.1186/s13100-018-0144-1](https://doi.org/10.1186/s13100-018-0144-1).
- Fuente vigente: [bitbucket.org/petrnovak/re_databases](https://bitbucket.org/petrnovak/re_databases) (Viridiplantae v4.0 + Metazoa v3.1).
- Formato: FASTA proteico de dominios + jerarquía taxonómica.
- Header: `>Cluster_Species__Class_I__LTR__Ty1_copia__Ale__RT` (jerarquía por `__`).
- Licencia: académica, cita.
- **Rol**: **fuente canónica de nomenclatura de linajes** (ver [`data/lineage_map.tsv`](../data/lineage_map.tsv)). Reclasificar el corpus con **TEsorter -db rexdb-plant**.

### Dfam 40.0 (release 2026-05-29)
- Paper base: Storer et al., *Mobile DNA* 12:2, 2021. DOI [10.1186/s13100-020-00230-y](https://doi.org/10.1186/s13100-020-00230-y).
- Descargas: [dfam.org/releases/current/families/FamDB](https://www.dfam.org/releases/current/families/FamDB/).
- Archivos relevantes:
  - `dfam40.0.h5.gz` (60 MB, índice global)
  - `dfam40.curated.consensus.0.h5.gz` (28 MB, consensos curados)
  - `dfam40.curated.hmm.*.h5.gz` (HMMs curados)
- Cobertura: pan-eucariota. Absorbió las familias derivadas de RepBase bajo CC0.
- Licencia: **CC0** (ideal para redistribución).
- **Rol**: fuente de LTR-RTs no vegetales (humano, insecto, etc.); filtrar por
  `classification LIKE '%LTR%'` y validar que trae 5'LTR–interno–3'LTR completo.

### RepetDB (INRAE URGI)
- Paper: Amselem et al., *Mobile DNA* 10:6, 2019. DOI [10.1186/s13100-019-0150-y](https://doi.org/10.1186/s13100-019-0150-y).
- Web: [urgi.versailles.inrae.fr/repetdb](https://urgi.versailles.inrae.fr/repetdb).
- Formato: FASTA + GFF; código Wicker de tres letras (RLC/RLG/RLX...).
- Cobertura: ~15 genomas vegetales/fúngicos.
- Licencia: INRAE open (Etalab).
- Descarga: sin dump masivo estable; el sitio está tras protección anti-scraping
  (Anubis). Requiere export interactivo o cliente autenticado. **Pendiente**.
- **Rol**: consensos plantas/hongos con Wicker estricto.

### TREP (Wicker)
- Web: [trep-db.uzh.ch](https://trep-db.uzh.ch/).
- Formato: FASTA con headers Wicker canónicos.
- ~4 000 secuencias, foco en trigo/gramíneas.
- Licencia: académica.
- **Rol**: validación cruzada de códigos Wicker.

### RepBase (GIRI)
- Paper: Bao, Kojima, Kohany, *Mobile DNA* 6:11, 2015.
- Web: [girinst.org/repbase](https://www.girinst.org/repbase).
- Licencia: **propietaria de pago desde 2019**; prohibida la redistribución.
- **Rol**: **descartada** para el corpus público de Inpactor3. Riesgo legal:
  entrenar y distribuir un modelo sobre secuencias RepBase puede constituir
  obra derivada. Uso interno solo con licencia institucional documentada.

### Ensembl Plants / Ensembl 2025
- FTP: [ftp.ensemblgenomes.org/pub/plants/current](https://ftp.ensemblgenomes.org/pub/plants/current/).
- Formato: genomas soft-masked + GFF3 de repeats por especie.
- Licencia: abierta.
- **Rol**: **sustrato** — reejecutar EDTA/LTR_retriever aquí y generar cajas
  reales para reemplazar los flancos aleatorios de la demo.

## Matriz de solapamiento

```
                Inpact PanTEon GyDB REXDB Dfam RepBase RepetDB Ensembl TREP
InpactorDB      -      HIGH    ref  ref   low  low     low     SRC     ref
PanTEon         HIGH   -       ref  ref   INT  low     low     SRC     ref
GyDB            ref    ref     -    FEEDS low  low     ref     -       ref
REXDB           ref    ref     FEEDS -    ref  ref     ref     -       ref
Dfam            low    INT     low  ref   -    ABS     med     med     ref
RepBase         low    low     low  ref   ABS  -       med     med     ref
RepetDB         low    low     ref  ref   med  med     -       med     ref
Ensembl         SRC    SRC     -    -     med  med     med     -       -
TREP            ref    ref     ref  ref   ref  ref     ref     -       -

IN=incluida-en, HIGH=solape alto, INT=integrada dentro, ABS=absorbe,
FEEDS=alimenta taxonomía, SRC=sustrato/genoma, ref=referencia, med/low=parcial
```

## Recomendación de fusión

**Corpus primario** (fusionar en el FASTA canónico):
1. InpactorDB V5 non-redundant.
2. PanTEon v1.6.2 (dedup contra 1).
3. Subset LTR de Dfam 40 curated (filtrar `LTR;*`, validar estructura).
4. RepetDB LTR (cuando se pueda descargar).

**Solo referencia** (no van al `X_train`):
- REXDB → mapa canónico de linajes; TEsorter para reclasificar homogéneamente.
- GyDB → HMMs para verificar dominios internos.
- TREP → validación Wicker.

**Descartar**:
- RepBase (licencia).
- APTEdb (redundante con Dfam y ya integrada en PanTEon).

## Pipeline de armonización

```
raw FASTAs → merge_corpus.py (dedup exacto + normalización taxonómica)
           → CD-HIT-est 95% / 80% cov (dedup difusa)   [siguiente PR]
           → TEsorter -db rexdb-plant (reclasificar)   [siguiente PR]
           → hmmscan GyDB (validar dominios)           [siguiente PR]
           → split estratificado por especie           [siguiente PR]
           → corpus final Inpactor3 v1 (Zenodo)
```

## Riesgo legal por fuente

| Base | Redistribuir seqs | Entrenar | Redistribuir pesos |
|---|---|---|---|
| InpactorDB | Sí (cita) | Sí | Sí |
| PanTEon | Sí (cita) | Sí | Sí |
| Dfam | Sí (CC0) | Sí | Sí |
| GyDB | Con permiso | Sí (académico) | Ambiguo |
| REXDB | Con cita | Sí | Sí |
| RepetDB | Sí (cita) | Sí | Sí |
| TREP | Con cita | Sí | Sí |
| **RepBase** | **No** | Solo con licencia | **Zona gris** |

Recomendación práctica: publicar Inpactor3 con pesos abiertos usando solo
CC0/CC-BY/Etalab → InpactorDB + PanTEon + Dfam + RepetDB + TREP.
