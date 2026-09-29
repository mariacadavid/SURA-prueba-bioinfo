#!/usr/bin/env bash
# 01_descarga.sh — Descarga metadatos y los VCF completos de chr13 y chr17 (30x, GRCh38).
# 30x GRCh38 lo usa ClinVar, PharmCAT y gnomAD.
# Si se corta internet, volver a correrlo: curl -C - continúa donde quedó.
set -euo pipefail
source config/config.sh

baja () {  # baja URL DESTINO (reanudable, con reintentos)
  echo ">> $(basename "$2")"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$2" "$1"
}

# 1) Metadatos
#    - 3202 muestras 30x con familia, sexo, población y superpoblación
baja "${BASE_URL}/20130606_g1k_3202_samples_ped_population.txt" "${META}/samples_3202_ped_population.txt"
#    - 2504 muestras de la Fase 3 (no emparentadas): se usan para análisis poblacionales
baja "${PHASE3_URL}/integrated_call_samples_v3.20130502.ALL.panel" "${META}/phase3_2504.panel"

# 2) VCF completos de chr13 y chr17 + índice (.tbi)
for c in $CHROMS_MAIN; do
  baja "$(vcf_url "$c")"      "${RAW}/${c}.vcf.gz"
  baja "$(vcf_url "$c").tbi"  "${RAW}/${c}.vcf.gz.tbi"
done

# 3) Verificación rápida
echo "---- Verificación"
head -3 "${META}/samples_3202_ped_population.txt"
echo "Muestras en panel Fase 3: $(tail -n +2 "${META}/phase3_2504.panel" | wc -l)"
for c in $CHROMS_MAIN; do
  n=$(bcftools query -l "${RAW}/${c}.vcf.gz" | wc -l)
  echo "${c}: ${n} muestras, $(du -h "${RAW}/${c}.vcf.gz" | cut -f1)"
done
echo ">> Descarga completa"
