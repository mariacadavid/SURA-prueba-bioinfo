#!/usr/bin/env bash
# 00_check.sh — Verifica que todas las herramientas respondan.
ok=0; fail=0
check () { if eval "$2" >/dev/null 2>&1; then echo "OK    $1"; ok=$((ok+1)); else echo "FALLA $1"; fail=$((fail+1)); fi; }
check "bcftools"   "bcftools --version"
check "tabix"      "tabix --version"
check "plink2"     "plink2 --version"
check "admixture"  "admixture --help | grep -qi admixture"
check "java"       "java -version"
check "snpEff"     "snpEff -version"
check "python"     "python -c 'import pandas, numpy, scipy, matplotlib, seaborn, statsmodels, gseapy'"
check "PharmCAT"   "java -jar tools/pharmcat/pharmcat.jar -version"
echo "---- $ok OK, $fail con falla"
# Prueba de red: leer la cabecera del VCF de chr17 directo del servidor de 1000 Genomas
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
if bcftools view -h "$URL" 2>/dev/null | tail -1 | grep -q "#CHROM"; then
  echo "OK    acceso remoto a 1000 Genomas"
else
  echo "FALLA acceso remoto a 1000 Genomas (revisar internet o el nombre del archivo)"
fi
