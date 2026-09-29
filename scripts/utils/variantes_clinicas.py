#!/usr/bin/env python
"""Variantes en genes de enfermedad accionables (bloque C).

Dos modos:
  sitios  <dir>  -> lee sitios.tsv, se queda con las variantes que ClinVar clasifica y
                    escribe posiciones.txt para que bcftools extraiga solo esos genotipos.
  reporte <dir>  -> lee genotipos.tsv, calcula frecuencias por población, aplica los filtros
                    clínicos y cruza los portadores con la ancestría.

Criterios de filtrado (explícitos y justificados en el informe):
  1. Clasificación de ClinVar = Pathogenic o Likely_pathogenic, sin clasificaciones en conflicto.
  2. Nivel de revisión >= 2 estrellas (criterios aportados, varios laboratorios, sin conflicto).
  3. Frecuencia GLOBAL (las 2.504 personas) compatible con la prevalencia de la enfermedad:
       - dominante (BRCA1, BRCA2, TP53, RB1): AF global <= 0,1 %  (es decir, <= ~5 portadores)
       - recesivo  (ATP7B, GAA): AF global <= 2 % (frecuencia de portadores)
     Se usa la frecuencia global y no la maxima por superpoblacion porque con n = 504-661 por
     superpoblacion un solo portador ya da AF ~ 0,08 %, y el filtro descartaria variantes raras
     legitimas. La AF maxima por superpoblacion se reporta igual, como senal de alerta adicional.
     Una variante "patogenica" mas frecuente que eso es casi con certeza una clasificacion
     erronea, no un hallazgo: contradice la prevalencia conocida de la enfermedad.
"""
import sys, os, re
import numpy as np
import pandas as pd

MODO, CLI = sys.argv[1], sys.argv[2]
TAB, FIG = "results/tablas", "results/figuras"
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)

GENES = {  # gen: (herencia, enfermedad, umbral de frecuencia alelica)
    "BRCA1":  ("dominante", "Cancer de mama y ovario hereditario", 0.001),
    "BRCA2":  ("dominante", "Cancer de mama y ovario hereditario", 0.001),
    "TP53":   ("dominante", "Sindrome de Li-Fraumeni", 0.001),
    "RB1":    ("dominante", "Retinoblastoma", 0.001),
    "ATP7B":  ("recesivo",  "Enfermedad de Wilson", 0.02),
    "GAA":    ("recesivo",  "Enfermedad de Pompe", 0.02),
}
# Nivel de revision de ClinVar -> estrellas
ESTRELLAS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_multiple_submitters": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
}
COLS = ["CHROM", "POS", "REF", "ALT", "GENEINFO", "CLNSIG", "CLNREVSTAT", "CLNDN", "MC"]


def gen_de(geneinfo, chrom, pos):
    """Gen al que pertenece la variante: primero por GENEINFO de ClinVar, si no por la region."""
    for g in GENES:
        if g in str(geneinfo):
            return g
    REG = {("chr13", 32313086, 32402268): "BRCA2", ("chr13", 48301710, 48601436): "RB1",
           ("chr13", 51928436, 52014198): "ATP7B", ("chr17", 7659779, 7689546): "TP53",
           ("chr17", 43042292, 43172245): "BRCA1", ("chr17", 80099526, 80121881): "GAA"}
    for (c, a, b), g in REG.items():
        if chrom == c and a <= int(pos) <= b:
            return g
    return "otro"


def estrellas(rev):
    return ESTRELLAS.get(str(rev).strip(), np.nan)


def impacto(mc):
    """Impacto funcional legible a partir del campo MC de ClinVar (terminos Sequence Ontology)."""
    s = str(mc)
    if s in (".", "nan", ""):
        return "no anotado"
    orden = ["nonsense", "frameshift", "splice_donor", "splice_acceptor", "start_lost",
             "stop_lost", "missense", "inframe", "splice_site", "synonymous",
             "5_prime_UTR", "3_prime_UTR", "intron", "non-coding"]
    for t in orden:
        if t in s:
            return {"nonsense": "stop ganado (nonsense)", "frameshift": "cambio de marco (frameshift)",
                    "splice_donor": "sitio donador de splicing", "splice_acceptor": "sitio aceptor de splicing",
                    "start_lost": "perdida del codon de inicio", "stop_lost": "perdida del codon de parada",
                    "missense": "missense", "inframe": "indel en marco", "splice_site": "sitio de splicing",
                    "synonymous": "sinonima", "5_prime_UTR": "UTR 5'", "3_prime_UTR": "UTR 3'",
                    "intron": "intronica", "non-coding": "no codificante"}[t]
    return s.split("|")[-1] if "|" in s else s


def carga_sitios():
    d = pd.read_csv(f"{CLI}/sitios.tsv", sep="\t", names=COLS, dtype=str, na_values=["."])
    d["POS"] = d.POS.astype(int)
    return d


# ---------------- MODO 1: elegir los sitios que ClinVar clasifica ----------------
if MODO == "sitios":
    d = carga_sitios()
    con_clinvar = d[d.CLNSIG.notna()].copy()
    con_clinvar[["CHROM", "POS"]].to_csv(f"{CLI}/posiciones.txt", sep="\t", header=False, index=False)
    print(f"Variantes de 1000G en los 6 genes: {len(d)}")
    print(f"  de ellas, con clasificacion en ClinVar: {len(con_clinvar)}")
    sys.exit(0)


# ---------------- MODO 2: reporte completo ----------------
sitios = carga_sitios()
sitios = sitios[sitios.CLNSIG.notna()].copy()

muestras = [l.strip() for l in open(f"{CLI}/muestras.txt") if l.strip()]
gt = pd.read_csv(f"{CLI}/genotipos.tsv", sep="\t", header=None,
                 names=["CHROM", "POS", "REF", "ALT"] + muestras, dtype=str)
gt["POS"] = gt.POS.astype(int)

# Metadatos de poblacion: se reutiliza la tabla del PCA (seccion B)
ped = pd.read_csv(f"{TAB}/B1_pca_coordenadas.tsv", sep="\t")[["IID", "pop", "superpop"]]
ped = ped.rename(columns={"pop": "poblacion", "superpop": "superpoblacion"})
mapa_sp = dict(zip(ped.IID, ped.superpoblacion))
mapa_pop = dict(zip(ped.IID, ped.poblacion))
SUPER = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# Conteo de alelos alternativos por individuo: 0, 1 o 2 (el VCF viene faseado: 0|1, 1|1, ...)
G = gt[muestras].apply(lambda col: col.str.count("1").fillna(0)).astype(int)
G.index = gt.index

filas, portadores = [], []
for i, r in gt.iterrows():
    dosis = G.loc[i]
    fila = {"CHROM": r.CHROM, "POS": r.POS, "REF": r.REF, "ALT": r.ALT}
    for sp in SUPER:
        ids = [m for m in muestras if mapa_sp.get(m) == sp]
        n = len(ids)
        ac = int(dosis[ids].sum())
        fila[f"AF_{sp}"] = round(ac / (2 * n), 6) if n else np.nan
        fila[f"port_{sp}"] = int((dosis[ids] > 0).sum())
    ids_clm = [m for m in muestras if mapa_pop.get(m) == "CLM"]
    fila["AF_CLM"] = round(int(dosis[ids_clm].sum()) / (2 * len(ids_clm)), 6) if ids_clm else np.nan
    fila["port_CLM"] = int((dosis[ids_clm] > 0).sum())
    fila["AF_max_superpob"] = max(fila[f"AF_{sp}"] for sp in SUPER)
    fila["AF_global"] = round(int(dosis.sum()) / (2 * len(muestras)), 6)
    fila["portadores_total"] = int((dosis > 0).sum())
    fila["homocigotos"] = int((dosis == 2).sum())
    filas.append(fila)
    for m in muestras:
        if dosis[m] > 0:
            portadores.append({"IID": m, "CHROM": r.CHROM, "POS": r.POS, "REF": r.REF, "ALT": r.ALT,
                               "copias": int(dosis[m]), "poblacion": mapa_pop.get(m),
                               "superpoblacion": mapa_sp.get(m)})
frec = pd.DataFrame(filas)

# Unir anotacion de ClinVar con las frecuencias
d = sitios.merge(frec, on=["CHROM", "POS", "REF", "ALT"], how="inner")
d["gen"] = [gen_de(gi, c, p) for gi, c, p in zip(d.GENEINFO, d.CHROM, d.POS)]
d["herencia"] = d.gen.map(lambda g: GENES.get(g, ("?",))[0])
d["enfermedad_gen"] = d.gen.map(lambda g: GENES.get(g, ("", ""))[1] if g in GENES else "")
d["estrellas"] = d.CLNREVSTAT.map(estrellas)
d["impacto_funcional"] = d.MC.map(impacto)
d["umbral_AF"] = d.gen.map(lambda g: GENES[g][2] if g in GENES else np.nan)

# --- Filtros, en cascada
d["f1_PLP"] = d.CLNSIG.fillna("").str.contains("Pathogenic", case=False) & \
              ~d.CLNSIG.fillna("").str.contains("Conflicting", case=False)
d["f2_estrellas>=2"] = d.estrellas >= 2
d["f3_AF_plausible"] = d.AF_global <= d.umbral_AF
d["pasa_todos"] = d.f1_PLP & d["f2_estrellas>=2"] & d.f3_AF_plausible

orden = ["gen", "herencia", "enfermedad_gen", "CHROM", "POS", "REF", "ALT", "impacto_funcional",
         "CLNSIG", "estrellas", "CLNREVSTAT", "CLNDN", "portadores_total", "homocigotos",
         "AF_global", "umbral_AF", "AF_max_superpob"] + [f"AF_{s}" for s in SUPER] + ["AF_CLM"] + \
        [f"port_{s}" for s in SUPER] + ["port_CLM", "f1_PLP", "f2_estrellas>=2", "f3_AF_plausible", "pasa_todos"]
d[orden].sort_values(["pasa_todos", "gen", "POS"], ascending=[False, True, True]) \
        .to_csv(f"{TAB}/C1_variantes_clinvar.tsv", sep="\t", index=False)

plp = d[d.f1_PLP].copy()
plp[orden].to_csv(f"{TAB}/C2_patogenicas.tsv", sep="\t", index=False)

# --- Embudo de filtrado (la tabla mas informativa del bloque)
embudo = pd.DataFrame([
    {"paso": "1. Variantes de 1000G con clasificacion en ClinVar", "variantes": len(d),
     "portadores": int(d.portadores_total.sum())},
    {"paso": "2. Clasificadas Pathogenic / Likely pathogenic (sin conflicto)", "variantes": int(d.f1_PLP.sum()),
     "portadores": int(d[d.f1_PLP].portadores_total.sum())},
    {"paso": "3. + nivel de revision >= 2 estrellas", "variantes": int((d.f1_PLP & d["f2_estrellas>=2"]).sum()),
     "portadores": int(d[d.f1_PLP & d["f2_estrellas>=2"]].portadores_total.sum())},
    {"paso": "4. + frecuencia compatible con la enfermedad", "variantes": int(d.pasa_todos.sum()),
     "portadores": int(d[d.pasa_todos].portadores_total.sum())},
])
embudo.to_csv(f"{TAB}/C3_embudo_filtrado.tsv", sep="\t", index=False)

# --- Portadores de variantes que pasan todos los filtros, por ancestria
port = pd.DataFrame(portadores)
finales = d[d.pasa_todos][["CHROM", "POS", "REF", "ALT", "gen", "CLNSIG", "estrellas", "impacto_funcional"]]
port_f = port.merge(finales, on=["CHROM", "POS", "REF", "ALT"], how="inner") if len(finales) else port.iloc[0:0]
port_f.to_csv(f"{TAB}/C4_portadores_por_ancestria.tsv", sep="\t", index=False)

# --- Impacto funcional: cuantas variantes de cada tipo
imp = (d.groupby(["gen", "impacto_funcional"]).size().rename("n").reset_index()
         .sort_values(["gen", "n"], ascending=[True, False]))
imp.to_csv(f"{TAB}/C5_impacto_funcional.tsv", sep="\t", index=False)

# --- Salida por pantalla
pd.set_option("display.width", 220)
print("\n=== Embudo de filtrado ===")
print(embudo.to_string(index=False))
print("\n=== Impacto funcional por gen (top 5) ===")
print(imp.groupby("gen").head(5).to_string(index=False))
print("\n=== Variantes P/LP descartadas por frecuencia (senal de clasificacion dudosa) ===")
desc = d[d.f1_PLP & d["f2_estrellas>=2"] & ~d.f3_AF_plausible]
cols_d = ["gen", "CHROM", "POS", "REF", "ALT", "impacto_funcional", "CLNSIG", "estrellas",
          "AF_global", "umbral_AF", "AF_max_superpob", "portadores_total"]
print(desc[cols_d].to_string(index=False) if len(desc) else "  (ninguna)")
print("\n=== Variantes que pasan TODOS los filtros ===")
cols_f = ["gen", "herencia", "CHROM", "POS", "REF", "ALT", "impacto_funcional", "CLNSIG",
          "estrellas", "portadores_total", "homocigotos"] + [f"port_{s}" for s in SUPER] + ["port_CLM"]
print(d[d.pasa_todos][cols_f].to_string(index=False) if d.pasa_todos.any() else "  (ninguna)")
if len(port_f):
    print("\n=== Portadores por superpoblacion ===")
    print(port_f.groupby(["superpoblacion", "gen"]).size().rename("portadores").to_string())

# --- Figura: embudo + portadores por superpoblacion
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
AZUL, NARANJA, TINTA, TINTA2, GRILLA = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "xtick.color": TINTA2,
                     "ytick.color": TINTA2, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
y = np.arange(len(embudo))[::-1]
axes[0].barh(y, embudo.variantes, color=AZUL, height=0.6)
for yi, v in zip(y, embudo.variantes):
    axes[0].text(v, yi, f" {v}", va="center", fontsize=9, color=TINTA)
axes[0].set_yticks(y)
axes[0].set_yticklabels([p.split(". ", 1)[1][:42] for p in embudo.paso], fontsize=8)
axes[0].set_xlabel("Variantes"); axes[0].set_title("Embudo de filtrado clinico", loc="left", color=TINTA)
axes[0].grid(axis="x", color=GRILLA, lw=0.6); axes[0].set_axisbelow(True)

if len(port_f):
    c = port_f.groupby("superpoblacion").size().reindex(SUPER).fillna(0)
    axes[1].bar(np.arange(len(SUPER)), c.values, color=NARANJA, width=0.6)
    for i, v in enumerate(c.values):
        axes[1].text(i, v, f"{int(v)}", ha="center", va="bottom", fontsize=9, color=TINTA)
    axes[1].set_xticks(np.arange(len(SUPER))); axes[1].set_xticklabels(SUPER)
    axes[1].set_ylabel("Portadores"); axes[1].grid(axis="y", color=GRILLA, lw=0.6); axes[1].set_axisbelow(True)
    axes[1].set_title("Portadores de variantes P/LP por superpoblacion", loc="left", color=TINTA)
else:
    axes[1].text(0.5, 0.5, "Ninguna variante paso todos los filtros", ha="center", va="center",
                 fontsize=10, color=TINTA2); axes[1].axis("off")
fig.tight_layout()
fig.savefig(f"{FIG}/C1_variantes_clinicas.png", bbox_inches="tight", dpi=150)
plt.close(fig)
