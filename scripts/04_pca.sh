#!/usr/bin/env bash
# 04_pca.sh — Prepara SNPs independientes y corre el PCA (bloque B).
#
# Decisiones:
#  - MAF >= 5%: las variantes comunes informan sobre ancestría; las raras agregan ruido y costo.
#  - Sin filtro de Hardy-Weinberg sobre todas las poblaciones juntas: mezclar poblaciones
#    rompe HWE por la propia estructura (efecto Wahlund) y quitaríamos justo los SNPs informativos.
#  - Poda por LD (ventana 50 SNPs, paso 5, r² > 0.2): PCA y ADMIXTURE asumen SNPs independientes;
#    bloques en LD pesarían de más.
#  - Se excluye la inversión 17q21.31 (chr17:45.4-46.8 Mb, GRCh38): región de LD de largo alcance
#    que puede dominar un componente principal sin reflejar ancestría global.
#  - Salida con cromosomas numéricos (--output-chr 26): ADMIXTURE no acepta "chr13".
set -euo pipefail
source config/config.sh
QC=${DATA}/qc
ANC=${DATA}/ancestria
mkdir -p "$ANC" results/tablas

printf "chr17\t45400000\t46800000\tinversion_17q21.31\n" > config/ld_largo_alcance.bed

for c in $CHROMS_MAIN; do
  echo ">> Poda LD ${c}"
  plink2 --pfile "${QC}/${c}" --maf 0.05 --geno 0.02 \
         --exclude bed1 config/ld_largo_alcance.bed \
         --indep-pairwise 50 5 0.2 --threads "$THREADS" --out "${ANC}/${c}_poda"
  plink2 --pfile "${QC}/${c}" --extract "${ANC}/${c}_poda.prune.in" \
         --make-pgen --out "${ANC}/${c}_podado"
done

# Unir los dos cromosomas en un solo archivo PLINK 1 (.bed), que usan PCA y ADMIXTURE
printf "%s\n" ${CHROMS_MAIN} | sed "s|^|${ANC}/|; s|$|_podado|" > "${ANC}/lista_merge.txt"
plink2 --pmerge-list "${ANC}/lista_merge.txt" pfile --make-bed --output-chr 26 \
       --out "${ANC}/podado_chr13_17"
echo ">> SNPs para ancestría: $(wc -l < "${ANC}/podado_chr13_17.bim" | tr -d ' ')"

# PCA: 10 componentes 
plink2 --bfile "${ANC}/podado_chr13_17" --pca 10 --threads "$THREADS" --out "${ANC}/pca"

python scripts/utils/graficos_ancestria.py pca
echo ">> Listo: results/figuras/B1_pca_*.png"
