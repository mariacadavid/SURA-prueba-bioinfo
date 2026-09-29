#!/usr/bin/env python
"""Pruebas CLM vs EUR, AFR, EAS por FENOTIPO ACCIONABLE, respetando la dirección clínica de cada gen.
Por qué: agrupar todo lo "no normal" mezcla efectos opuestos (p. ej. CYP2C19 pobre vs ultrarrápido)
y en CYP3A5 lo accionable es ser EXPRESOR (normal/intermedio), no lo contrario.

Uso: python scripts/utils/pruebas_pgx.py results/tablas/D2_fenotipos_por_grupo.tsv
Salida: results/tablas/D4_pruebas_CLM_vs_otros.tsv
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

RED = {"Pobre", "Intermedio"}
# gen -> [(fenotipo accionable, categorías que lo componen, fármaco de ejemplo)]
ACCIONABLE = {
    "CYP2C19": [("función disminuida (PM+IM)", RED, "clopidogrel"),
                ("función aumentada (RM+UM)", {"Rápido", "Ultrarrápido"}, "IBP, citalopram/escitalopram")],
    "CYP2C9":  [("función disminuida (PM+IM)", RED, "warfarina, fenitoína, AINE")],
    "CYP3A5":  [("expresor (NM+IM)", {"Normal", "Intermedio"}, "tacrolimus (requiere más dosis)")],
    "DPYD":    [("función disminuida (PM+IM)", RED, "5-FU, capecitabina")],
    "TPMT":    [("función disminuida (PM+IM)", RED, "tiopurinas")],
    "NUDT15":  [("función disminuida (PM+IM)", RED, "tiopurinas")],
    "SLCO1B1": [("función disminuida o pobre", RED, "simvastatina y otras estatinas")],
    "UGT1A1":  [("función disminuida (PM+IM)", RED, "irinotecán, atazanavir")],
    "VKORC1":  [("portador rs9923231-T (sensibilidad)", RED, "warfarina (menor dosis)")],
}

def wilson(k, n, z=1.96):
    """Intervalo de confianza del 95 % de Wilson para una proporción (mejor que el normal con n pequeño)."""
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (max(0, c - h), min(1, c + h))

def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out

def pruebas(cnt):
    """cnt: tabla con columnas gen, grupo, categoria, n (conteos por fenotipo). Excluye indeterminados."""
    filas = []
    for gen, lista in ACCIONABLE.items():
        sub = cnt[(cnt.gen == gen) & (cnt.categoria != "Indeterminado")]
        def cuenta(g, cats):
            s = sub[sub.grupo == g]; k = int(s[s.categoria.isin(cats)].n.sum()); return k, int(s.n.sum())
        for fen, cats, farmaco in lista:
            k1, n1 = cuenta("CLM", cats)
            for otro in ["EUR", "AFR", "EAS"]:
                k2, n2 = cuenta(otro, cats)
                OR, p = fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])
                l1, u1 = wilson(k1, n1); l2, u2 = wilson(k2, n2)
                filas.append({"gen": gen, "fenotipo_accionable": fen, "farmaco_ejemplo": farmaco,
                              "comparacion": f"CLM vs {otro}",
                              "CLM_%": round(100 * k1 / n1, 1), "CLM_IC95": f"{100*l1:.1f}-{100*u1:.1f}",
                              "otro_%": round(100 * k2 / n2, 1), "otro_IC95": f"{100*l2:.1f}-{100*u2:.1f}",
                              "n_CLM": n1, "n_otro": n2,
                              "odds_ratio": round(OR, 2) if np.isfinite(OR) else OR, "p": p})
    pr = pd.DataFrame(filas)
    pr["q_FDR"] = bh_fdr(pr.p.values)
    pr["significativo_q<0.05"] = pr.q_FDR < 0.05
    return pr

if __name__ == "__main__":
    cnt = pd.read_csv(sys.argv[1], sep="\t")
    pr = pruebas(cnt)
    salida = sys.argv[2] if len(sys.argv) > 2 else "results/tablas/D4_pruebas_CLM_vs_otros.tsv"
    pr.to_csv(salida, sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(pr.drop(columns=["farmaco_ejemplo", "n_CLM", "n_otro"]).to_string(index=False))
