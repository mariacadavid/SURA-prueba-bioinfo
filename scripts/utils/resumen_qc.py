#!/usr/bin/env python
"""Resume el QC en tablas para el informe (bloque A).
Uso: python scripts/utils/resumen_qc.py data/qc data/meta results/tablas
"""
import sys, glob, os, re
import pandas as pd

qc, meta, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)

# --- Variantes por cromosoma y clases de frecuencia
filas = []
for f in sorted(glob.glob(f"{qc}/chr*.afreq")):
    c = os.path.basename(f).split(".")[0]
    fr = pd.read_csv(f, sep="\t")
    maf = fr["ALT_FREQS"].where(fr["ALT_FREQS"] <= 0.5, 1 - fr["ALT_FREQS"])
    # tipos de variante antes de filtrar (bcftools stats, sección SN)
    sn = {}
    with open(f"{qc}/{c}.bcftools_stats.txt") as h:
        for line in h:
            if line.startswith("SN\t"):
                _, _, k, v = line.rstrip("\n").split("\t")
                sn[k.rstrip(":")] = int(v)
    filas.append({
        "cromosoma": c,
        "variantes_totales": sn.get("number of records"),
        "SNPs_todos": sn.get("number of SNPs"),
        "indels": sn.get("number of indels"),
        "multialelicos": sn.get("number of multiallelic sites"),
        "SNPs_bialelicos_QC": len(fr),
        "monomorficos_en_2504": int((maf == 0).sum()),
        "comunes_MAF>=5%": int((maf >= 0.05).sum()),
        "baja_frec_1-5%": int(((maf >= 0.01) & (maf < 0.05)).sum()),
        "raras_<1%": int(((maf > 0) & (maf < 0.01)).sum()),
        "singletons": int(((ac := (fr["ALT_FREQS"] * fr["OBS_CT"]).round()) == 1).sum()
                          + ((fr["OBS_CT"] - ac) == 1).sum()),
    })
tv = pd.DataFrame(filas)
tv.to_csv(f"{out}/A1_variantes_por_cromosoma.tsv", sep="\t", index=False)
print(tv.to_string(index=False))

# --- Datos faltantes
sm = pd.concat([pd.read_csv(f, sep="\t").assign(cromosoma=os.path.basename(f).split(".")[0])
                for f in sorted(glob.glob(f"{qc}/chr*.smiss"))])
col_id = "#IID" if "#IID" in sm.columns else "IID"
miss = sm.groupby(col_id)["F_MISS"].max()
faltantes = pd.DataFrame([{
    "individuos_analizados": miss.size,
    "F_MISS_maximo": miss.max(),
    "individuos_con_F_MISS>2%": int((miss > 0.02).sum()),
}])
faltantes.to_csv(f"{out}/A2_datos_faltantes.tsv", sep="\t", index=False)
print(faltantes.to_string(index=False))

# --- Individuos por población y continente
ped = pd.read_csv(f"{meta}/samples_3202_ped_population.txt", sep=r"\s+")
ped.columns = [re.sub(r"[^a-z]", "", c.lower()) for c in ped.columns]  # tolera variaciones del encabezado
id_col = next(c for c in ped.columns if "sample" in c)
pop_col = next(c for c in ped.columns if c == "population")
sup_col = next(c for c in ped.columns if "superpop" in c)
unrel = set(open(f"{qc}/unrelated_2504.txt").read().split())
ped["no_emparentado"] = ped[id_col].isin(unrel)

cont = (ped.groupby(sup_col)
          .agg(muestras_3202=(id_col, "size"), no_emparentados=("no_emparentado", "sum"))
          .assign(proporcion=lambda d: (d.no_emparentados / d.no_emparentados.sum()).round(3))
          .reset_index().rename(columns={sup_col: "superpoblacion"}))
cont.to_csv(f"{out}/A3_individuos_por_continente.tsv", sep="\t", index=False)
print(cont.to_string(index=False))

pops = (ped[ped.no_emparentado].groupby([sup_col, pop_col]).size()
          .rename("n").reset_index().rename(columns={sup_col: "superpoblacion", pop_col: "poblacion"}))
pops.to_csv(f"{out}/A4_individuos_por_poblacion.tsv", sep="\t", index=False)

# Individuos del panel de la Fase 3 que faltan en el VCF 30x
faltan = unrel - set(ped[id_col])
print(f"Individuos de la Fase 3 ausentes en la metadata 30x: {len(faltan)}")
