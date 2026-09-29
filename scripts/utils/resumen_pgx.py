#!/usr/bin/env python
"""Resume los resultados de PharmCAT (módulo Farmacogenómica).
Uso: python scripts/utils/resumen_pgx.py <carpeta_salida_pharmcat> <metadata_3202>

Genera:
  D1_diplotipos_AMR.tsv            diplotipo y fenotipo por individuo AMR y gen (entregable pedido)
  D2_fenotipos_por_grupo.tsv       frecuencia de cada fenotipo por población/superpoblación
  D3_diplotipos_frecuentes_AMR.tsv diplotipos más frecuentes por gen en CLM, MXL, PEL, PUR
  D4_pruebas_CLM_vs_otros.tsv      CLM vs EUR, AFR, EAS: Fisher (fenotipo no normal vs normal) + FDR
  D1_fenotipos_CLM_vs_otros.png    figura comparativa
"""
import sys, glob, os, re
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR, META = sys.argv[1], sys.argv[2]
TAB, FIG = "results/tablas", "results/figuras"
os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)
GENES = ["CYP2C19", "CYP2C9", "CYP3A5", "DPYD", "TPMT", "NUDT15", "SLCO1B1", "UGT1A1", "VKORC1"]
CATS = ["Pobre", "Intermedio", "Normal", "Rápido", "Ultrarrápido", "Indeterminado"]
GRUPOS = ["CLM", "MXL", "PEL", "PUR", "EUR", "AFR", "EAS"]

def bh_fdr(p):
    """Corrección de Benjamini-Hochberg (tasa de falsos descubrimientos)."""
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def categoria(gen, fenotipo, diplotipo):
    """Lleva el fenotipo CPIC de PharmCAT a una escala común ordenada.
    SLCO1B1 usa 'función' (pobre/disminuida/normal/aumentada) -> misma escala.
    VKORC1 no tiene fenotipo en PharmCAT: se deriva de rs9923231 (-1639G>A; C>T en GRCh38):
    0 copias T = Normal, 1 = Intermedio (sensibilidad aumentada a warfarina), 2 = Pobre (alta sensibilidad)."""
    d = str(diplotipo or "")
    if gen == "VKORC1":
        if "rs9923231" not in d: return "Indeterminado"
        return {0: "Normal", 1: "Intermedio", 2: "Pobre"}[d.count("variant")]
    p = str(fenotipo or "").lower().strip()
    if not p or p in ("nan", "n/a", "no result") or "indeterminate" in p: return "Indeterminado"
    if "ultrarapid" in p: return "Ultrarrápido"
    if "rapid" in p or "increased" in p: return "Rápido"
    if "poor" in p: return "Pobre"
    if "intermediate" in p or "decreased" in p: return "Intermedio"
    if "normal" in p: return "Normal"
    return "Indeterminado"

# --- 1) Leer un report.tsv por persona
filas = []
for f in glob.glob(f"{DIR}/*.report.tsv"):
    muestra = os.path.basename(f).split(".")[-3]
    t = pd.read_csv(f, sep="\t", skiprows=1, dtype=str)
    t = t[t["Gene"].isin(GENES)]
    for _, r in t.iterrows():
        filas.append({"IID": muestra, "gen": r["Gene"], "diplotipo": r["Source Diplotype"],
                      "fenotipo_CPIC": r.get("Phenotype"), "puntaje_actividad": r.get("Activity Score"),
                      "categoria": categoria(r["Gene"], r.get("Phenotype"), r["Source Diplotype"])})
if not filas: sys.exit(f"No se encontraron *.report.tsv en {DIR}")
pgx = pd.DataFrame(filas)

ped = pd.read_csv(META, sep=r"\s+")
ped.columns = [re.sub(r"[^a-z]", "", c.lower()) for c in ped.columns]
idc = next(c for c in ped.columns if "sample" in c); supc = next(c for c in ped.columns if "superpop" in c)
ped = ped.rename(columns={idc: "IID", "population": "poblacion", supc: "superpoblacion"})
pgx = pgx.merge(ped[["IID", "poblacion", "superpoblacion"]], on="IID", how="left")
pgx["grupo"] = np.where(pgx.superpoblacion == "AMR", pgx.poblacion, pgx.superpoblacion)
print(f"Individuos procesados: {pgx.IID.nunique()}")

# --- 2) D1: diplotipos por individuo AMR
amr = pgx[pgx.superpoblacion == "AMR"].sort_values(["poblacion", "IID", "gen"])
amr[["IID", "poblacion", "gen", "diplotipo", "fenotipo_CPIC", "puntaje_actividad", "categoria"]] \
    .to_csv(f"{TAB}/D1_diplotipos_AMR.tsv", sep="\t", index=False)

# --- 3) D2: frecuencia de fenotipos por grupo
cnt = (pgx.groupby(["gen", "grupo", "categoria"]).size().rename("n").reset_index())
cnt["total"] = cnt.groupby(["gen", "grupo"]).n.transform("sum")
cnt["proporcion"] = (cnt.n / cnt.total).round(4)
cnt.to_csv(f"{TAB}/D2_fenotipos_por_grupo.tsv", sep="\t", index=False)

# --- 4) D3: diplotipos más frecuentes en AMR
d3 = (amr.groupby(["gen", "poblacion", "diplotipo"]).size().rename("n").reset_index())
d3["proporcion"] = (d3.n / d3.groupby(["gen", "poblacion"]).n.transform("sum")).round(3)
d3 = d3.sort_values(["gen", "poblacion", "n"], ascending=[True, True, False]).groupby(["gen", "poblacion"]).head(5)
d3.to_csv(f"{TAB}/D3_diplotipos_frecuentes_AMR.tsv", sep="\t", index=False)

# --- 5) D4: CLM vs EUR, AFR, EAS por fenotipo accionable (dirección clínica de cada gen) + IC95 + FDR
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pruebas_pgx import pruebas
pr = pruebas(cnt)
pr.to_csv(f"{TAB}/D4_pruebas_CLM_vs_otros.tsv", sep="\t", index=False)
print(pr.drop(columns=["farmaco_ejemplo", "n_CLM", "n_otro"]).to_string(index=False))

# --- 6) Figura: fenotipos por grupo, un panel por gen (escala divergente centrada en Normal)
COL = {"Pobre": "#a52a2a", "Intermedio": "#e8927c", "Normal": "#d9d8d3",
       "Rápido": "#86b6ef", "Ultrarrápido": "#1c5cab", "Indeterminado": "#ffffff"}
TINTA, TINTA2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c9c8c2", "xtick.color": TINTA2, "ytick.color": TINTA2})
fig, axes = plt.subplots(3, 3, figsize=(12, 8.5), sharex=True)
for ax, g in zip(axes.flat, GENES):
    tab = (cnt[cnt.gen == g].pivot(index="grupo", columns="categoria", values="proporcion")
             .reindex(index=GRUPOS, columns=CATS).fillna(0))
    y = np.arange(len(GRUPOS))[::-1]
    izq = np.zeros(len(GRUPOS))
    for c in CATS:
        v = tab[c].values
        indet = c == "Indeterminado"
        ax.barh(y, v, left=izq, color=COL[c], edgecolor="#a8a7a1" if indet else "white",
                linewidth=0.8 if indet else 1.5, height=0.72, hatch="////" if indet else None)
        for yi, xi, vi in zip(y, izq, v):
            if vi >= 0.12 and c != "Normal":
                ax.text(xi + vi / 2, yi, f"{vi*100:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if c in ("Pobre", "Ultrarrápido") else TINTA)
        izq += v
    ax.set_yticks(y); ax.set_yticklabels(GRUPOS)
    for lbl in ax.get_yticklabels():
        if lbl.get_text() == "CLM": lbl.set_fontweight("bold"); lbl.set_color(TINTA)
    ax.set_title(g, loc="left", fontsize=10, color=TINTA, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xticklabels(["0", "25", "50", "75", "100%"])
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, fc=COL[c], ec="#a8a7a1", hatch="////" if c == "Indeterminado" else None) for c in CATS]
fig.legend(handles, ["Pobre", "Intermedio", "Normal", "Rápido", "Ultrarrápido", "Indeterminado / sin llamado"],
           ncol=6, loc="upper center", bbox_to_anchor=(0.5, 0.985), frameon=False)
fig.suptitle("Fenotipos farmacogenéticos (CPIC) por población — números = % de individuos",
             x=0.01, y=1.02, ha="left", color=TINTA, fontsize=11)
fig.text(0.01, -0.035, "Porcentajes sobre el total de individuos, incluidos los indeterminados. SLCO1B1 en escala de función "
         "(pobre/disminuida/normal/aumentada). VKORC1 según rs9923231: Pobre = TT (alta sensibilidad a warfarina), Intermedio = CT.\n"
         "CYP3A5: 'Pobre' (*3/*3) es el genotipo más común en EUR; el fenotipo clínicamente accionable es ser expresor "
         "(Normal/Intermedio), que requiere más tacrolimus.", fontsize=7.5, color=TINTA2, va="top")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(f"{FIG}/D1_fenotipos_CLM_vs_otros.png", bbox_inches="tight", dpi=150); plt.close(fig)
