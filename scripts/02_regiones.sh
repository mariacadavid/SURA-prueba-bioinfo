#!/usr/bin/env bash
# 02_regiones.sh — Extraer solo las regiones necesarias de otros cromosomas, sin descargarlos completos:
#   (a) los 9 genes de farmacogenómica (coordenadas oficiales de PharmCAT, GRCh38)
#   (b) los 12 SNPs de nutrigenómica (posiciones GRCh38 consultadas en Ensembl)
# bcftools lee el VCF remoto usando su índice y baja únicamente esas regiones.
set -euo pipefail
source config/config.sh

# (a) Regiones PGx: filtramos el BED de PharmCAT a los 9 genes de interes
#     POI = posición adicional que PharmCAT usa para CYP2C (rs12777823), se conserva.
# Si PharmCAT se instaló como .jar suelto (versión anterior del setup), bajamos el paquete completo,
# que trae el BED de regiones y el preprocesador de VCF.
if [ ! -f tools/pharmcat/pharmcat_regions.bed ]; then
  echo ">> Descargando paquete PharmCAT 3.4.0 (regiones + preprocesador)..."
  mkdir -p tools/pharmcat
  curl -fL -o tools/pharmcat-pipeline.tar.gz \
    "https://github.com/PharmGKB/PharmCAT/releases/download/v3.4.0/pharmcat-pipeline-3.4.0.tar.gz"
  tar -xzf tools/pharmcat-pipeline.tar.gz -C tools/pharmcat && rm tools/pharmcat-pipeline.tar.gz
  pip install -q -r tools/pharmcat/requirements.txt
fi
GENES="DPYD|UGT1A1|TPMT|CYP3A5|CYP2C19|CYP2C9|SLCO1B1|NUDT15|VKORC1"
grep -E "PX=(${GENES})$|POI" tools/pharmcat/pharmcat_regions.bed > config/pgx_regions.bed
echo ">> Regiones PGx:"; cat config/pgx_regions.bed

# (b) SNPs de nutrigenómica: posiciones GRCh38 desde la API de Ensembl -> config/nutri_snps.tsv
python scripts/utils/ensembl_positions.py config/nutri_snps.tsv
awk 'NR>1{print $2"\t"$3-1"\t"$3"\t"$1}' config/nutri_snps.tsv > config/nutri_regions.bed
echo ">> SNPs de nutrigenómica:"; cat config/nutri_snps.tsv

# (c) Extraer por cromosoma y unir
extrae () {  # extrae BED SALIDA
  local bed=$1 out=$2 parts=()
  for c in $(cut -f1 "$bed" | sort -u); do
    local tmp="${REG}/tmp_$(basename "$out" .vcf.gz)_${c}.vcf.gz"
    echo "   ${c}..."
    awk -v c="$c" '$1==c' "$bed" > "${REG}/tmp_${c}.bed"   # (awk: el grep de Mac no tiene -P)
    bcftools view -R "${REG}/tmp_${c}.bed" -Oz -o "$tmp" "$(vcf_url "$c")"
    parts+=("$tmp")
  done
  bcftools concat "${parts[@]}" -Ou | bcftools sort -Oz -o "$out"
  bcftools index -t "$out"
  rm -f "${parts[@]}" "${REG}"/tmp_*.bed
  echo "   -> $out: $(bcftools view -H "$out" | wc -l | tr -d ' ') variantes"
}

echo ">> Extrayendo PGx"
extrae config/pgx_regions.bed   "${REG}/pgx_regions.vcf.gz"
echo ">> Extrayendo nutrigenómica"
extrae config/nutri_regions.bed "${REG}/nutri_snps.vcf.gz"
# Los índices remotos (.tbi) que bcftools baja a la carpeta actual se pueden borrar
rm -f ./*.vcf.gz.tbi
echo ">> Listoo"
