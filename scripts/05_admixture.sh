#!/usr/bin/env bash
# 05_admixture.sh — ADMIXTURE para K = 2..7 con validación cruzada (bloque B).
# Tardo como 3.5 horas. 
#
# Decisiones
#  - Submuestra aleatoria de ADMIX_SNPS SNPs (semilla fija): ADMIXTURE escala con #SNPs x #personas;
#    ~60 mil SNPs independientes son suficientes para ancestría continental y la me da para correrlo en poquitas horas en mi compu portatil. 
#  - --cv (validación cruzada 5 veces): el K con menor error CV es el más apoyado por los datos.
#  - -s SEED: la misma semilla da el mismo resultado (reproducibilidad) ps: esta en config.
set -euo pipefail
source config/config.sh
ANC=${DATA}/ancestria
ADMIX_SNPS=${ADMIX_SNPS:-60000}
K_MIN=${K_MIN:-2}
K_MAX=${K_MAX:-7}

plink2 --bfile "${ANC}/podado_chr13_17" --thin-count "$ADMIX_SNPS" --seed "$SEED" \
       --make-bed --out "${ANC}/admix_input"

cd "$ANC"
for K in $(seq "$K_MIN" "$K_MAX"); do
  if [ -f "admix_input.${K}.Q" ]; then echo ">> K=${K} ya existe, se salta"; continue; fi
  echo ">> ADMIXTURE K=${K}  ($(date +%H:%M))"
  admixture --cv -s "$SEED" -j"$THREADS" admix_input.bed "$K" > "admixture_K${K}.log"
  grep "CV error" "admixture_K${K}.log"
done
cd - > /dev/null

grep -h "CV error" "${ANC}"/admixture_K*.log | sed -E 's/.*K=([0-9]+)\): (.*)/\1\t\2/' \
  | sort -n | { echo -e "K\tCV_error"; cat; } > results/tablas/B2_admixture_cv.tsv
cat results/tablas/B2_admixture_cv.tsv
python scripts/utils/graficos_ancestria.py admixture
echo ">> Listo: results/figuras/B3_admixture_*.png"
