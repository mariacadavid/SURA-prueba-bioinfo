#!/usr/bin/env bash
# 06_farmacogenomica.sh — Diplotipos y fenotipos PGx con PharmCAT (módulo Farmacogenómica).
# Tarda: 30-90 min. La primera vez PharmCAT descarga el genoma de referencia GRCh38 (~1 GB).
#
# Decisiones:
#  - PharmCAT 3.4.0: implementa las definiciones de alelos de PharmVar/CPIC y asigna fenotipos según CPIC.
#  - Se corre en las 2.504 personas no emparentadas (todas las poblaciones), para comparar CLM con
#    EUR, AFR y EAS con el mismo método. Los resultados se reportan con foco en AMR (CLM, MXL, PEL, PUR).
#  - --absent-to-ref: el VCF 30x solo lista posiciones donde alguien tiene una variante; una posición
#    ausente significa "igual a la referencia". Sin esta opción PharmCAT la trataría como dato faltante
#    y casi ningún gen tendría llamado. Limitación: una posición de baja cobertura también estaría ausente.
#  - CYP2D6 no se tipifica desde VCF (deleciones/duplicaciones del gen); requiere BAM y herramientas
#    como Cyrius o StellarPGx. No está entre los 9 genes pedidos.
set -euo pipefail
source config/config.sh
PGX=${DATA}/pgx
mkdir -p "$PGX" results/tablas results/figuras

cp "${DATA}/qc/unrelated_2504.txt" "${PGX}/muestras.txt"

python tools/pharmcat/pharmcat_pipeline "${REG}/pgx_regions.vcf.gz" \
  -S "${PGX}/muestras.txt" \
  --absent-to-ref \
  -g CYP2C19,CYP2C9,CYP3A5,DPYD,TPMT,NUDT15,SLCO1B1,UGT1A1,VKORC1 \
  -reporterCallsOnlyTsv \
  -cp "$THREADS" \
  -o "${PGX}/pharmcat"

python scripts/utils/resumen_pgx.py "${PGX}/pharmcat" "${META}/samples_3202_ped_population.txt"
echo ">> Listo: results/tablas/D*.tsv y results/figuras/D*.png"