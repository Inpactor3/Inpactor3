# Comparación Inpactor2 vs Inpactor3 sobre PanTEon LTR

- Secuencias evaluadas: **30**
- Accuracy Inpactor2 (superfamilia): **0.267**
- Accuracy Inpactor3 (superfamilia): **0.133**
- Concordancia entre ambos: **0.467**

## Tabla completa

| seqid | especie | verdad | Inpactor2 lin | Inpactor2 SF | Inpactor3 lin | Inpactor3 SF |
|---|---|---|---|---|---|---|
| PDB00002101#CLASSI/L | Esox_lucius | **LARD** | RLG/REINA | GYPSY ✗ | REINA | GYPSY ✗ |
| PDB00000324#CLASSI/L | Anabarilius_grahami | **GYPSY** | RLC/ALE/RETROFIT | COPIA ✗ | REINA | GYPSY ✓ |
| PDB00001853#CLASSI/L | Ahaetulla_prasina | **LARD** | RLC/ALE/RETROFIT | COPIA ✗ | REINA | GYPSY ✗ |
| BEL-52_DRe#CLASSI/LT | Danio_rerio | **BELPAO** | RLG/GALADRIEL | GYPSY ✗ | REINA | GYPSY ✗ |
| ERV1-N6_DR#CLASSI/LT | Danio_rerio | **ERV** | RLC/ALE/RETROFIT | COPIA ✗ | ORYCO | COPIA ✗ |
| PDB00004205#CLASSI/L | Electrophorus_electr | **LARD** | RLC/TAR/TORK | COPIA ✗ | ORYCO | COPIA ✗ |
| PDB00002103#CLASSI/L | Esox_lucius | **TRIM** | RLC/ALE/RETROFIT | COPIA ✗ | ORYCO | COPIA ✗ |
| MER34B-int#CLASSI/LT | Balaenoptera_acutoro | **ERV** | RLC/TAR/TORK | COPIA ✗ | ORYCO | COPIA ✗ |
| PDB00003964#CLASSI/L | Anastrepha_obliqua | **GYPSY** | RLG/CRM | GYPSY ✓ | ORYCO | COPIA ✗ |
| PDB00001223#CLASSI/L | Samia_ricini | **COPIA** | RLG/TEKAY/DEL | GYPSY ✗ | ORYCO | COPIA ✓ |
| PDB00000764#CLASSI/L | Triplophysa_siluroid | **TRIM** | RLC/ALE/RETROFIT | COPIA ✗ | REINA | GYPSY ✗ |
| PDB00001859#CLASSI/L | Ahaetulla_prasina | **COPIA** | RLC/ALE/RETROFIT | COPIA ✓ | REINA | GYPSY ✗ |
| BEL1_DR#CLASSI/LTR/B | Danio_rerio | **BELPAO** | RLC/IVANA/ORYCO | COPIA ✗ | REINA | GYPSY ✗ |
| BEL-56_DRe#CLASSI/LT | Danio_rerio | **BELPAO** | RLG/CRM | GYPSY ✗ | TORK | COPIA ✗ |
| PDB00002250#CLASSI/L | Tinamus_guttatus | **TRIM** | RLC/ALE/RETROFIT | COPIA ✗ | ALE | COPIA ✗ |
| PDB00001297#CLASSI/L | Talpa_occidentalis | **LARD** | RLG/TAT | GYPSY ✗ | REINA | GYPSY ✗ |
| PDB00004011#CLASSI/L | Anastrepha_obliqua | **COPIA** | RLC/IVANA/ORYCO | COPIA ✓ | ORYCO | COPIA ✓ |
| MicOch-1.327#CLASSI/ | Microtus_ochrogaster | **ERV** | RLC/ALE/RETROFIT | COPIA ✗ | ORYCO | COPIA ✗ |
| PDB00000657#CLASSI/L | Anarsia_innoxiella | **LARD** | RLC/ALE/RETROFIT | COPIA ✗ | REINA | GYPSY ✗ |
| PDB00001257#CLASSI/L | Samia_ricini | **COPIA** | RLC/ALE/RETROFIT | COPIA ✓ | REINA | GYPSY ✗ |
| PDB00000236#CLASSI/L | Anabarilius_grahami | **TRIM** | RLG/REINA | GYPSY ✗ | TORK | COPIA ✗ |
| PDB00004758#CLASSI/L | Heterocephalus_glabe | **TRIM** | RLC/ALE/RETROFIT | COPIA ✗ | REINA | GYPSY ✗ |
| TolMat-5.333#CLASSI/ | Tolypeutes_matacus | **ERV** | RLG/REINA | GYPSY ✗ | ORYCO | COPIA ✗ |
| OndZib-1.1#CLASSI/LT | Ondatra_zibethicus | **ERV** | RLC/ALE/RETROFIT | COPIA ✗ | ORYCO | COPIA ✗ |
| PDB00004048#CLASSI/L | Anastrepha_obliqua | **GYPSY** | RLG/REINA | GYPSY ✓ | ORYCO | COPIA ✗ |
| BEL13_DR#CLASSI/LTR/ | Danio_rerio | **BELPAO** | RLC/ALE/RETROFIT | COPIA ✗ | ORYCO | COPIA ✗ |
| BEL-48_DRe#CLASSI/LT | Danio_rerio | **BELPAO** | RLG/REINA | GYPSY ✗ | REINA | GYPSY ✗ |
| Gypsy36_DR#CLASSI/LT | Danio_rerio | **GYPSY** | RLG/ATHILA | GYPSY ✓ | ORYCO | COPIA ✗ |
| PDB00002540#CLASSI/L | Conomelus_anceps | **GYPSY** | RLG/ATHILA | GYPSY ✓ | ORYCO | COPIA ✗ |
| PDB00003017#CLASSI/L | Ricordea_florida | **COPIA** | RLC/ALE/RETROFIT | COPIA ✓ | ORYCO | COPIA ✓ |

## Confusión Inpactor2 (filas=verdad, columnas=predicho)

| ↓verdad \ pred→ | BELPAO | COPIA | ERV | GYPSY | LARD | TRIM |
|---|---|---|---|---|---|---|
| **BELPAO** | 0 | 2 | 0 | 3 | 0 | 0 |
| **COPIA** | 0 | 4 | 0 | 1 | 0 | 0 |
| **ERV** | 0 | 4 | 0 | 1 | 0 | 0 |
| **GYPSY** | 0 | 1 | 0 | 4 | 0 | 0 |
| **LARD** | 0 | 3 | 0 | 2 | 0 | 0 |
| **TRIM** | 0 | 4 | 0 | 1 | 0 | 0 |

## Confusión Inpactor3 (filas=verdad, columnas=predicho)

| ↓verdad \ pred→ | BELPAO | COPIA | ERV | GYPSY | LARD | TRIM |
|---|---|---|---|---|---|---|
| **BELPAO** | 0 | 2 | 0 | 3 | 0 | 0 |
| **COPIA** | 0 | 3 | 0 | 2 | 0 | 0 |
| **ERV** | 0 | 5 | 0 | 0 | 0 | 0 |
| **GYPSY** | 0 | 4 | 0 | 1 | 0 | 0 |
| **LARD** | 0 | 1 | 0 | 4 | 0 | 0 |
| **TRIM** | 0 | 3 | 0 | 2 | 0 | 0 |
