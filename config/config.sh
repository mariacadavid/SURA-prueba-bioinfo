# config.sh — rutas y parámetros comunes. Todos los scripts lo cargan con: source config/config.sh
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"
VCF_DIR="${BASE_URL}/working/20220422_3202_phased_SNV_INDEL_SV"
PHASE3_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# Nombre del VCF de un cromosoma (uso: vcf_url chr17)
vcf_url () { echo "${VCF_DIR}/1kGP_high_coverage_Illumina.$1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"; }

CHROMS_MAIN="chr13 chr17"            # cromosomas completos que pide la prueba
SEED=20260926                        # semilla fija para todo lo aleatorio (ADMIXTURE, simulaciones)
THREADS=6                            # núcleos a usar (el MacBook Air M2 tiene 8)

DATA=data
RAW=${DATA}/raw
META=${DATA}/meta
REG=${DATA}/regiones
LOGS=logs
mkdir -p "$RAW" "$META" "$REG" "$LOGS" results
