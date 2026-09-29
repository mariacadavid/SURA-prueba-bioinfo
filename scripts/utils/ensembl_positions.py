#!/usr/bin/env python
"""Consulta en Ensembl (GRCh38) la posición de los SNPs de nutrigenómica y la guarda en un TSV.
Por qué: el VCF 30x no trae rsIDs, así que hay que buscar cada SNP por cromosoma y posición.
Uso: python scripts/utils/ensembl_positions.py config/nutri_snps.tsv
"""
import json, sys, urllib.request

SNPS = {  # rsID: gen (lista de la prueba)
    "rs9939609": "FTO", "rs17782313": "MC4R", "rs7903146": "TCF7L2",
    "rs738409": "PNPLA3", "rs58542926": "TM6SF2", "rs1260326": "GCKR",
    "rs662799": "APOA5", "rs174547": "FADS1", "rs1801133": "MTHFR",
    "rs4988235": "LCT", "rs429358": "APOE", "rs7412": "APOE",
}

req = urllib.request.Request(
    "https://rest.ensembl.org/variation/homo_sapiens",
    data=json.dumps({"ids": list(SNPS)}).encode(),
    headers={"Content-Type": "application/json", "Accept": "application/json"},
)
data = json.load(urllib.request.urlopen(req, timeout=60))

rows = []
for rs, gene in SNPS.items():
    maps = [m for m in data[rs]["mappings"]
            if m["assembly_name"] == "GRCh38" and m["seq_region_name"].isalnum()
            and "_" not in m["seq_region_name"]]
    m = maps[0]
    rows.append((rs, "chr" + m["seq_region_name"], m["start"], m["allele_string"], gene))

rows.sort(key=lambda r: (int(r[1][3:]) if r[1][3:].isdigit() else 99, r[2]))
with open(sys.argv[1], "w") as f:
    f.write("rsid\tchrom\tpos\talleles\tgene\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
