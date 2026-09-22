# Comparación Inpactor2 vs Inpactor3

- Referencia: **inpactor2** (14 cajas)
- Evaluado : **inpactor3** (385 cajas)
- Umbral IoU 1D: **0.3**
- Cajas que coinciden posicionalmente: **12**
- Sólo en referencia (FN): **2**
- Sólo en inpactor3 (FP): **373**

## Resumen por herramienta

### Inpactor2
- Predicciones totales: **14**
- Longitud media: 6359 bp
- Linajes distintos: 8
- Top linajes:
    - ALE               4
    - TORK              2
    - TAT               2
    - ORYCO             2
    - SIRE              1
    - REINA             1
    - CRM               1
    - ATHILA            1

### Inpactor3
- Predicciones totales: **385**
- Longitud media: 2683 bp
- Linajes distintos: 4
- Top linajes:
    - REINA           232
    - ORYCO           106
    - ALE              39
    - TORK              8

## Métricas por linaje (TP requiere IoU y linaje correctos)

| Linaje | TP | FP | FN | Precisión | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| ALE | 0 | 39 | 4 | 0.000 | 0.000 | 0.000 |
| ATHILA | 0 | 0 | 1 | 0.000 | 0.000 | 0.000 |
| CRM | 0 | 0 | 1 | 0.000 | 0.000 | 0.000 |
| ORYCO | 0 | 106 | 2 | 0.000 | 0.000 | 0.000 |
| REINA | 1 | 231 | 0 | 0.004 | 1.000 | 0.009 |
| SIRE | 0 | 0 | 1 | 0.000 | 0.000 | 0.000 |
| TAT | 0 | 0 | 2 | 0.000 | 0.000 | 0.000 |
| TORK | 0 | 8 | 2 | 0.000 | 0.000 | 0.000 |

**Macro promedio**: P=0.001  R=0.125  F1=0.001

## Confusión de linaje (filas = referencia, columnas = predicho)

| ref \ pred | ALE | CRM | ORYCO | REINA | SIRE | TAT | TORK |
|---|---|---|---|---|---|---|---|
| **ALE** | 0 | 0 | 0 | 4 | 0 | 0 | 0 |
| **CRM** | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| **ORYCO** | 1 | 0 | 0 | 1 | 0 | 0 | 0 |
| **REINA** | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| **SIRE** | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| **TAT** | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| **TORK** | 1 | 0 | 0 | 1 | 0 | 0 | 0 |
