#!/usr/bin/env python
"""Figuras de ancestría (bloque B).
Uso: python scripts/utils/graficos_ancestria.py pca        (tras 04_pca.sh)
     python scripts/utils/graficos_ancestria.py admixture  (tras 05_admixture.sh)
"""
import sys, glob, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

ANC, META, FIG, TAB = "data/ancestria", "data/meta", "results/figuras", "results/tablas"
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)

# Paleta categórica (orden fijo) + marcadores como segunda codificación (no depender solo del color)
PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
MARCAS = ["o", "s", "^", "D", "v", "P", "X", "*"]
SUPER = ["AFR", "AMR", "EAS", "EUR", "SAS"]
AMR = ["CLM", "MXL", "PEL", "PUR"]
TINTA, TINTA2, GRILLA, GRIS = "#0b0b0b", "#52514e", "#e6e5e0", "#c9c8c2"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRIS, "axes.labelcolor": TINTA2,
                     "xtick.color": TINTA2, "ytick.color": TINTA2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})

def metadata():
    ped = pd.read_csv(f"{META}/samples_3202_ped_population.txt", sep=r"\s+")
    ped.columns = [re.sub(r"[^a-z]", "", c.lower()) for c in ped.columns]
    idc = next(c for c in ped.columns if "sample" in c)
    supc = next(c for c in ped.columns if "superpop" in c)
    return ped.rename(columns={idc: "IID", "population": "pop", supc: "superpop"})[["IID", "pop", "superpop"]]

def etiqueta_centroide(ax, x, y, texto):
    ax.annotate(texto, (np.median(x), np.median(y)), fontsize=9, fontweight="bold", color=TINTA,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

def pca():
    ev = pd.read_csv(f"{ANC}/pca.eigenvec", sep=r"\s+")
    ev = ev.rename(columns={ev.columns[0]: "IID"}) if "IID" not in ev.columns else ev
    if "#FID" in ev.columns: ev = ev.drop(columns="#FID")
    ev["IID"] = ev["IID"].astype(str)
    ev = ev.merge(metadata(), on="IID", how="left")
    val = np.loadtxt(f"{ANC}/pca.eigenval")
    var = val / val.sum() * 100  # % de la varianza entre los 10 PCs calculados
    ev.to_csv(f"{TAB}/B1_pca_coordenadas.tsv", sep="\t", index=False)
    pd.DataFrame({"PC": [f"PC{i+1}" for i in range(len(val))], "autovalor": val,
                  "pct_varianza_entre_10PC": var.round(2)}).to_csv(f"{TAB}/B1_pca_autovalores.tsv", sep="\t", index=False)

    # Figura 1: todas las superpoblaciones, PC1-PC2 y PC3-PC4
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (a, b) in zip(axes, [(1, 2), (3, 4)]):
        for i, sp in enumerate(SUPER):
            d = ev[ev.superpop == sp]
            ax.scatter(d[f"PC{a}"], d[f"PC{b}"], s=10, c=PALETA[i], marker=MARCAS[i],
                       alpha=0.75, linewidths=0, label=f"{sp} (n={len(d)})")
            etiqueta_centroide(ax, d[f"PC{a}"], d[f"PC{b}"], sp)
        ax.set_xlabel(f"PC{a} ({var[a-1]:.1f}%)"); ax.set_ylabel(f"PC{b} ({var[b-1]:.1f}%)")
        ax.grid(color=GRILLA, lw=0.6); ax.set_axisbelow(True)
    axes[0].legend(frameon=False, fontsize=8, markerscale=1.8, loc="best")
    fig.suptitle("PCA de 2.504 individuos no emparentados, chr13 + chr17", color=TINTA, x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/B1_pca_superpoblaciones.png", bbox_inches="tight"); plt.close(fig)

    # Figura 2: poblaciones AMR sobre el resto en gris
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    otros = ev[ev.superpop != "AMR"]
    ax.scatter(otros.PC1, otros.PC2, s=6, c=GRIS, linewidths=0, alpha=0.5, label="Otras superpoblaciones")
    for sp in ["AFR", "EUR", "EAS", "SAS"]:
        d = ev[ev.superpop == sp]; etiqueta_centroide(ax, d.PC1, d.PC2, sp)
    for i, p in enumerate(AMR):
        d = ev[ev["pop"] == p]
        ax.scatter(d.PC1, d.PC2, s=16, c=PALETA[i], marker=MARCAS[i], linewidths=0.4,
                   edgecolors="white", label=f"{p} (n={len(d)})")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.grid(color=GRILLA, lw=0.6); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, markerscale=1.5)
    ax.set_title("Poblaciones latinoamericanas (AMR) en el espacio del PCA", color=TINTA, loc="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/B1_pca_amr.png", bbox_inches="tight"); plt.close(fig)
    print("PCA: % varianza (entre 10 PCs):", dict(zip([f"PC{i+1}" for i in range(4)], var[:4].round(1))))

def admixture():
    fam = pd.read_csv(f"{ANC}/admix_input.fam", sep=r"\s+", header=None)
    ids = fam[1].astype(str)
    meta = metadata().set_index("IID").loc[ids].reset_index()
    orden_sp = {s: i for i, s in enumerate(["AFR", "EUR", "EAS", "SAS", "AMR"])}
    NOMBRES = {"AFR": "Africana", "EUR": "Europea", "EAS": "Asia oriental", "SAS": "Asia del sur",
               "AMR": "Indígena americana"}

    def etiquetar(q):
        """Nombre de cada componente: ancestría de la superpoblación donde es máximo + población con el valor más alto."""
        qdf = pd.DataFrame(q)
        por_sp = qdf.groupby(meta.superpop.values).mean(); por_pop = qdf.groupby(meta["pop"].values).mean()
        return [f"{NOMBRES.get(por_sp[j].idxmax(), por_sp[j].idxmax())} (máx. en {por_pop[j].idxmax()})"
                for j in range(q.shape[1])]
    qs = sorted(glob.glob(f"{ANC}/admix_input.*.Q"), key=lambda f: int(f.split(".")[-2]))
    if not qs: sys.exit("No hay archivos .Q todavía")

    # Orden de individuos: superpoblación, población y componente principal (con el K mayor)
    qmax = pd.read_csv(qs[-1], sep=r"\s+", header=None)
    meta["dom"] = qmax.values.argmax(1); meta["domv"] = qmax.values.max(1)
    meta["o_sp"] = meta.superpop.map(orden_sp)
    orden = meta.sort_values(["o_sp", "pop", "dom", "domv"], ascending=[True, True, True, False]).index

    fig, axes = plt.subplots(len(qs), 1, figsize=(13, 1.25 * len(qs) + 1), sharex=True)
    axes = np.atleast_1d(axes)
    resumen = []
    prev = None  # matriz Q (ya alineada) del K anterior
    for ax, f in zip(axes, qs):
        K = int(f.split(".")[-2])
        q = pd.read_csv(f, sep=r"\s+", header=None).values
        if prev is None:
            # primer K: componentes ordenados por la superpoblación donde son máximos
            medias = pd.DataFrame(q).groupby(meta.superpop.values).mean()
            comp = sorted(range(K), key=lambda k: (orden_sp.get(medias[k].idxmax(), 9), -medias[k].max()))
        else:
            # K siguientes: cada componente hereda el color del componente del K anterior con el que
            # más se correlaciona (asignación húngara); el componente nuevo toma el siguiente color.
            # Así un mismo color = misma ancestría en todas las filas.
            corr = np.corrcoef(prev.T, q.T)[:prev.shape[1], prev.shape[1]:]
            fila, col = linear_sum_assignment(-np.nan_to_num(corr))
            comp = list(col[np.argsort(fila)]) + [k for k in range(K) if k not in col]
        q = q[:, comp]
        prev = q
        qo = q[orden]
        base = np.zeros(len(qo))
        for j in range(K):
            ax.bar(np.arange(len(qo)), qo[:, j], bottom=base, width=1.0, color=PALETA[j % 8], linewidth=0)
            base += qo[:, j]
        ax.set_ylim(0, 1); ax.set_yticks([]); ax.set_ylabel(f"K={K}", rotation=0, ha="right", va="center", color=TINTA)
        ax.tick_params(axis="x", length=0)
        for sp in ax.spines.values(): sp.set_visible(False)
        # proporciones medias por población, en formato largo, con el nombre del componente EN ESTE K
        et_k = etiquetar(q)
        m = pd.DataFrame(q, columns=et_k); m["pop"] = meta["pop"].values
        for p_, fila_ in m.groupby("pop").mean().iterrows():
            for j, (c, v) in enumerate(fila_.items()):
                resumen.append({"K": K, "poblacion": p_, "color": j + 1, "componente": c, "proporcion": round(v, 3)})

    etiquetas = etiquetar(prev)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=PALETA[j % 8]) for j in range(len(etiquetas))]
    fig.legend(handles, etiquetas, ncol=min(4, len(etiquetas)), loc="upper center", bbox_to_anchor=(0.5, 0.0),
               frameon=False, fontsize=9,
               title=f"Componente de ancestría (nombres según K={prev.shape[1]}; en K menores un color puede agrupar "
                     "ancestrías que luego se separan)", title_fontsize=9)

    # separadores y nombres de población
    pops = meta.loc[orden, ["pop", "superpop"]].reset_index(drop=True)
    cortes = pops.index[pops["pop"] != pops["pop"].shift()].tolist() + [len(pops)]
    for ax in axes:
        for c in cortes[1:-1]: ax.axvline(c, color="white", lw=0.8)
    centros = [(a + b) / 2 for a, b in zip(cortes[:-1], cortes[1:])]
    axes[-1].set_xticks(centros)
    axes[-1].set_xticklabels([pops.loc[c, "pop"] for c in cortes[:-1]], rotation=90, fontsize=8)
    for lbl, c in zip(axes[-1].get_xticklabels(), cortes[:-1]):
        if pops.loc[c, "pop"] == "CLM": lbl.set_fontweight("bold"); lbl.set_color(TINTA)
    fig.suptitle("ADMIXTURE: proporciones de ancestría por individuo (chr13 + chr17)", color=TINTA, x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/B3_admixture_barras.png", bbox_inches="tight"); plt.close(fig)
    res = pd.DataFrame(resumen)
    res.to_csv(f"{TAB}/B3_admixture_medias_por_poblacion.tsv", sep="\t", index=False)
    ancho = res[res.poblacion.isin(AMR)].pivot_table(index=["K", "componente"], columns="poblacion",
                                                     values="proporcion").reindex(columns=AMR)
    ancho.to_csv(f"{TAB}/B3_admixture_AMR_resumen.tsv", sep="\t")
    print(ancho.round(3).to_string())

    cv = f"{TAB}/B2_admixture_cv.tsv"
    if os.path.exists(cv):
        d = pd.read_csv(cv, sep="\t")
        fig, ax = plt.subplots(figsize=(5, 3.2))
        ax.plot(d.K, d.CV_error, color=PALETA[0], lw=2, marker="o", ms=8, mec="white", mew=1.5)
        best = d.loc[d.CV_error.idxmin()]
        ax.annotate(f"mínimo: K={int(best.K)}", (best.K, best.CV_error), xytext=(0, 12),
                    textcoords="offset points", ha="center", fontsize=9, color=TINTA)
        ax.set_xlabel("K (poblaciones ancestrales)"); ax.set_ylabel("Error de validación cruzada")
        ax.grid(color=GRILLA, lw=0.6); ax.set_axisbelow(True)
        fig.tight_layout(); fig.savefig(f"{FIG}/B2_admixture_cv.png", bbox_inches="tight"); plt.close(fig)
        print(d.to_string(index=False))

if __name__ == "__main__":
    {"pca": pca, "admixture": admixture}[sys.argv[1]]()
