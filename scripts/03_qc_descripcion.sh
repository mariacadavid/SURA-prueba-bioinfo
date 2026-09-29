#!/usr/bin/env bash
# 03_qc_descripcion.sh — Control de calidad y descripción del dataset (bloque A).
# Uso: conda activate bioinfo-prueba && bash scripts/03_qc_descripcion.sh
#
# Decisiones (y por qué):
#  - Solo las 2.504 personas no emparentadas de la Fase 3: los 602 tríos agregados en la versión 30x
#    incluyen hijos, y los parientes inflan frecuencias y distorsionan PCA/ADMIXTURE.
#  - Solo SNPs bialélicos A/C/G/T: PCA, ADMIXTURE y MAF asumen dos alelos; indels y SVs se describen aparte.
#  - FILTER=PASS: el proyecto ya marcó las variantes de baja calidad.
#  - Categorías de frecuencia (MAF): común >= 5%, baja frecuencia 1-5%, rara < 1%. Son los cortes estándar.
set -euo pipefail
source config/config.sh
QC=${DATA}/qc
mkdir -p "$QC" results/tablas

# 1) Lista de no emparentados = IDs del panel de la Fase 3 (2.504)
tail -n +2 "${META}/phase3_2504.panel" | cut -f1 > "${QC}/unrelated_2504.txt"

for c in $CHROMS_MAIN; do
  echo ">> ${c}"
  # 2) Conteo de TODOS los tipos de variante antes de filtrar (para la descripción)
  bcftools view -S "${QC}/unrelated_2504.txt" --force-samples -f PASS,. -Ou "${RAW}/${c}.vcf.gz" \
    | bcftools stats - > "${QC}/${c}.bcftools_stats.txt"

  # 3) PLINK2: SNPs bialélicos, no emparentados, IDs únicos chr:pos:ref:alt
  plink2 --vcf "${RAW}/${c}.vcf.gz" \
         --keep "${QC}/unrelated_2504.txt" \
         --var-filter \
         --snps-only just-acgt --max-alleles 2 \
         --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 50 truncate \
         --threads "$THREADS" \
         --make-pgen --out "${QC}/${c}"

  # 4) Frecuencias, faltantes por persona y por variante, Hardy-Weinberg
  plink2 --pfile "${QC}/${c}" --freq --missing --hardy --threads "$THREADS" --out "${QC}/${c}"
done

# 5) Resumen en tablas
python scripts/utils/resumen_qc.py "$QC" "$META" results/tablas
echo ">> Tablas en results/tablas/"
