#!/usr/bin/env bash
# Demo end-to-end de Inpactor3 para presentacion (~2 min en CPU).
#
# Muestra en orden:
#   1) Auditoria de InpactorDB V5 (ya descargada)
#   2) Corpus canonico fusionado (65605 secuencias)
#   3) Arquitectura Yoro1D entrenando 2 epocas
#   4) Inferencia con NMS sobre 2 ventanas de validacion

set -e
cd "$(dirname "$0")/.."

echo "============================================================"
echo "  Inpactor3 - demo de arquitectura YORO 1D"
echo "============================================================"
echo

echo "[1/3] Auditoria InpactorDB V5 (67k LTR-RTs, 181 especies)"
echo "------------------------------------------------------------"
python3 scripts/audit_datasets.py --which inpactordb_nr 2>&1 | tail -25

echo
echo "[2/3] Corpus canonico fusionado"
echo "------------------------------------------------------------"
echo "  Total secuencias: $(grep -c '^>' data/corpus/inpactor3_v0.fasta)"
echo "  Primer header:    $(head -1 data/corpus/inpactor3_v0.fasta)"
echo "  Linajes canonicos:"
grep '^>' data/corpus/inpactor3_v0.fasta \
    | sed 's/.*|lin=\([^|]*\).*/\1/' \
    | sort | uniq -c | sort -rn | head -13 \
    | awk '{printf "    %-12s %6d\n", $2, $1}'

echo
echo "[3/3] Entrenamiento Yoro1D (2 epocas, 128 ventanas, CPU)"
echo "------------------------------------------------------------"
python3 scripts/train_demo.py \
    --fasta data/corpus/inpactor3_v0.fasta \
    --epochs 2 --batch 2 --limit-ltrs 500 --windows 64 --val-windows 8

echo
echo "============================================================"
echo "  OK - checkpoint guardado en models/inpactor3_demo.pt"
echo "============================================================"
