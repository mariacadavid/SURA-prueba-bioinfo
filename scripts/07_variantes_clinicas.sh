#!/usr/bin/env bash
# 07_variantes_clinicas.sh — Variantes en genes de enfermedad accionables (bloque C).
# note: lo que pesa es descargar ClinVar, ~100 MB
#
# Decisiones (y por qué):
#  - Genes: los de la lista ACMG SF (hallazgos secundarios accionables) que están en chr13 y chr17.
#    Criterio publicado y verificable, en lugar de una selección a gusto.
#  - Clasificación clínica, enfermedad, nivel de revisión (estrellas) e IMPACTO FUNCIONAL: todo sale
#    del VCF de ClinVar (campos CLNSIG, CLNDN, CLNREVSTAT y MC). No se necesita VEP ni SnpEff.
#  - Solo las 2.504 personas no emparentadas, SNPs e indels bialélicos con FILTER=PASS.
set -euo pipefail
source config/config.sh
CLI=${DATA}/clinicas
mkdir -p "$CLI" results/tablas results/figuras

# 1) Regiones de los 6 genes (GRCh38, coordenadas de Ensembl + 2 kb de margen a cada lado)
cat > config/genes_clinicos.bed <<'EOF'
chr13	32313086	32402268	BRCA2
chr13	48301710	48601436	RB1
chr13	51928436	52014198	ATP7B
chr17	7659779	7689546	TP53
chr17	43042292	43172245	BRCA1
chr17	80099526	80121881	GAA
EOF
sed 's/^chr//' config/genes_clinicos.bed > "${CLI}/genes_sin_chr.bed"   # ClinVar usa "13", no "chr13"
printf "13\tchr13\n17\tchr17\n" > "${CLI}/rename_chrs.txt"

# 2) Descargar ClinVar (GRCh38) si no está
CV=${CLI}/clinvar.vcf.gz
if [ ! -s "$CV" ]; then
  echo ">> Descargando ClinVar (GRCh38)..."
  curl -fL --retry 5 --retry-delay 10 -C - -o "$CV" \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  curl -fL --retry 5 --retry-delay 10 -C - -o "${CV}.tbi" \
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi
# Guardar la fecha de la versión usada (reproducibilidad)
bcftools view -h "$CV" | grep -i "^##fileDate\|^##source" > "${CLI}/clinvar_version.txt" || true
cat "${CLI}/clinvar_version.txt"

# 3) Recortar ClinVar a los 6 genes y renombrar cromosomas a "chrN"
bcftools view -R "${CLI}/genes_sin_chr.bed" "$CV" -Ou \
  | bcftools annotate --rename-chrs "${CLI}/rename_chrs.txt" -Oz -o "${CLI}/clinvar_genes.vcf.gz"
bcftools index -f -t "${CLI}/clinvar_genes.vcf.gz"
echo ">> Variantes de ClinVar en los 6 genes: $(bcftools view -H "${CLI}/clinvar_genes.vcf.gz" | wc -l | tr -d ' ')"

# 4) Recortar 1000 Genomas a los mismos genes
for c in chr13 chr17; do
  awk -v c="$c" '$1==c' config/genes_clinicos.bed > "${CLI}/tmp_${c}.bed"
  bcftools view -R "${CLI}/tmp_${c}.bed" -S "${DATA}/qc/unrelated_2504.txt" --force-samples \
    -f PASS,. -m2 -M2 -Oz -o "${CLI}/g1k_${c}.vcf.gz" "${RAW}/${c}.vcf.gz"
  bcftools index -f -t "${CLI}/g1k_${c}.vcf.gz"
done
bcftools concat "${CLI}/g1k_chr13.vcf.gz" "${CLI}/g1k_chr17.vcf.gz" -Oz -o "${CLI}/g1k_genes.vcf.gz"
bcftools index -f -t "${CLI}/g1k_genes.vcf.gz"
echo ">> Variantes de 1000G en los 6 genes: $(bcftools view -H "${CLI}/g1k_genes.vcf.gz" | wc -l | tr -d ' ')"

# 5) Anotar con ClinVar (bcftools empareja por cromosoma, posición, REF y ALT)
bcftools annotate -a "${CLI}/clinvar_genes.vcf.gz" \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/MC,INFO/GENEINFO,INFO/ALLELEID \
  "${CLI}/g1k_genes.vcf.gz" -Oz -o "${CLI}/anotado.vcf.gz"
bcftools index -f -t "${CLI}/anotado.vcf.gz"

# 6) Tabla de sitios (sin genotipos, para no cargar 2.504 columnas de todo)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/MC\n' \
  "${CLI}/anotado.vcf.gz" > "${CLI}/sitios.tsv"

# 7) Python elige los sitios con clasificación en ClinVar y escribe sus posiciones
python scripts/utils/variantes_clinicas.py sitios "$CLI"

# 8) Genotipos solo de esos sitios (pocos, así que es liviano)
bcftools view -R "${CLI}/posiciones.txt" "${CLI}/anotado.vcf.gz" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > "${CLI}/genotipos.tsv"
bcftools query -l "${CLI}/anotado.vcf.gz" > "${CLI}/muestras.txt"

# 9) Reporte: frecuencias por población, filtros y cruce con ancestría
python scripts/utils/variantes_clinicas.py reporte "$CLI"
echo ">> Listo: results/tablas/C*.tsv y results/figuras/C1_*.png"
