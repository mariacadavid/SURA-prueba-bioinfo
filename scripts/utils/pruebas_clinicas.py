#!/usr/bin/env python
"""Pruebas de agrupamiento por ancestria de los portadores de variantes P/LP (bloque C).

Pregunta: la tasa de portadores de variantes patogenicas o probablemente patogenicas,
?es distinta entre superpoblaciones? El enunciado pide decir si las variantes con
potencial patogenico se agrupan por tipo de ancestria.

Que hace, y por que cada prueba:
 1. Chi-cuadrado de independencia sobre la tabla 5 x 2 (superpoblacion x portador/no portador).
    Es la prueba global: una sola hipotesis nula, "la tasa es la misma en las cinco".
 2. Fisher exacta para las comparaciones por pares contra EUR. Con 3-4 portadores en una
    celda la aproximacion del chi-cuadrado deja de ser confiable; la exacta no la necesita.
 3. Correccion Benjamini-Hochberg sobre esas cuatro pruebas por pares (q_FDR), porque mirar
    cuatro comparaciones y quedarse con la menor infla el error tipo I.
 4. Prueba de permutaciones sobre la tabla completa, como respaldo del chi-cuadrado: baraja
    la etiqueta de superpoblacion y cuenta cuantas veces el estadistico barajado iguala o
    supera el observado. No asume frecuencias esperadas grandes, que es justo el supuesto
    que aqui queda al limite.

Una persona que lleve dos variantes P/LP cuenta como UN portador, no dos.

Uso:
  python scripts/utils/pruebas_clinicas.py
  python scripts/utils/pruebas_clinicas.py results/tablas/C4_portadores_por_ancestria.tsv \
                                           results/tablas/B1_pca_coordenadas.tsv
Salida: results/tablas/C6_pruebas_ancestria.tsv (y el resumen por pantalla)
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

SUPER = ["AFR", "AMR", "EAS", "EUR", "SAS"]
REFERENCIA = "EUR"          # grupo contra el que se comparan los demas
SEMILLA = 42                # las permutaciones son aleatorias: se fija para que sea reproducible
N_PERM = 100_000

TAB = "results/tablas"
IN_PORT = f"{TAB}/C4_portadores_por_ancestria.tsv"
IN_META = f"{TAB}/B1_pca_coordenadas.tsv"
OUT = f"{TAB}/C6_pruebas_ancestria.tsv"


def wilson(k, n, z=1.96):
    """IC 95 % de Wilson para una proporcion. Con k = 3 y n = 489 el IC normal se sale de [0, 1]."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def bh_fdr(p):
    """Benjamini-Hochberg. Misma implementacion que pruebas_pgx.py, para que ambos modulos coincidan."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def chi2_2col(k, n, p):
    """Chi-cuadrado de Pearson para una tabla C x 2, en forma cerrada.

    Con solo dos columnas el estadistico se reduce a sum_i (k_i - n_i p)^2 / (n_i p (1-p)),
    donde p es la tasa global de portadores. Vale para una tabla o para muchas a la vez,
    que es lo que permite hacer las permutaciones sin un bucle.
    """
    return (((k - n * p) ** 2) / (n * p * (1 - p))).sum(axis=-1)


def permutacion(k, n, n_perm=N_PERM, semilla=SEMILLA):
    """p por permutaciones (prueba exacta aproximada) para la tabla completa.

    Bajo la hipotesis nula los portadores se reparten al azar entre las personas,
    manteniendo fijos el total de portadores y el tamano de cada superpoblacion. Eso es
    exactamente una hipergeometrica multivariante, asi que en vez de barajar un vector de
    2.504 individuos se muestrea directamente el reparto de portadores. Devuelve la
    proporcion de repartos cuyo estadistico iguala o supera el observado, con la correccion
    (x + 1) / (m + 1), que evita reportar p = 0 cuando ninguno lo alcanza.

    No depende del supuesto de frecuencias esperadas grandes, que es el que aqui queda justo
    en el limite, de modo que sirve de respaldo al chi-cuadrado.
    """
    rng = np.random.default_rng(semilla)
    n = np.asarray(n, float)
    p = k.sum() / n.sum()
    obs = chi2_2col(np.asarray(k, float), n, p)
    sim = rng.multivariate_hypergeometric(n.astype(int), int(k.sum()), size=n_perm).astype(float)
    extremos = int((chi2_2col(sim, n, p) >= obs - 1e-12).sum())
    return (extremos + 1) / (n_perm + 1)


def tabla_conteos(portadores, meta):
    """Portadores unicos y tamano de muestra por superpoblacion."""
    n = (meta.drop_duplicates("IID")
             .superpoblacion.value_counts()
             .reindex(SUPER).fillna(0).astype(int))
    k = (portadores.drop_duplicates("IID")      # una persona con dos variantes es UN portador
                   .superpoblacion.value_counts()
                   .reindex(SUPER).fillna(0).astype(int))
    return k.values, n.values


def pruebas(k, n):
    tab = np.array([[ki, ni - ki] for ki, ni in zip(k, n)])
    chi2, p_glob, dof, esperadas = chi2_contingency(tab, correction=False)

    i_ref = SUPER.index(REFERENCIA)
    filas = []
    for i, sp in enumerate(SUPER):
        if sp == REFERENCIA:
            continue
        OR, p = fisher_exact([[k[i_ref], n[i_ref] - k[i_ref]], [k[i], n[i] - k[i]]])
        l1, u1 = wilson(k[i_ref], n[i_ref])
        l2, u2 = wilson(k[i], n[i])
        filas.append({"comparacion": f"{REFERENCIA} vs {sp}",
                      f"{REFERENCIA}_%": round(100 * k[i_ref] / n[i_ref], 2),
                      f"{REFERENCIA}_IC95": f"{100*l1:.2f}-{100*u1:.2f}",
                      "otro_%": round(100 * k[i] / n[i], 2),
                      "otro_IC95": f"{100*l2:.2f}-{100*u2:.2f}",
                      "portadores_otro": int(k[i]), "n_otro": int(n[i]),
                      "odds_ratio": round(OR, 2) if np.isfinite(OR) else OR,
                      "p": p})
    pares = pd.DataFrame(filas)
    pares["q_FDR"] = bh_fdr(pares.p.values)
    pares["significativo_q<0.05"] = pares.q_FDR < 0.05
    return chi2, p_glob, dof, esperadas, pares


if __name__ == "__main__":
    f_port = sys.argv[1] if len(sys.argv) > 1 else IN_PORT
    f_meta = sys.argv[2] if len(sys.argv) > 2 else IN_META
    salida = sys.argv[3] if len(sys.argv) > 3 else OUT

    portadores = pd.read_csv(f_port, sep="\t")
    meta = pd.read_csv(f_meta, sep="\t")
    if "superpoblacion" not in meta.columns:        # B1 trae la columna como "superpop"
        meta = meta.rename(columns={"superpop": "superpoblacion"})

    k, n = tabla_conteos(portadores, meta)
    chi2, p_glob, dof, esperadas, pares = pruebas(k, n)
    p_perm = permutacion(k, n)

    pd.set_option("display.width", 250)
    print("=== Portadores de variantes P/LP por superpoblacion ===")
    print(pd.DataFrame({"portadores": k, "n": n,
                        "tasa_%": (100 * k / n).round(2),
                        "esperados_si_H0": esperadas[:, 0].round(1)},
                       index=SUPER).to_string())
    filas_dobles = len(portadores) - portadores.IID.nunique()
    if filas_dobles:
        print(f"\nNota: {filas_dobles} fila(s) corresponden a personas con mas de una variante P/LP; "
              "se cuentan una sola vez.")

    print(f"\n=== Prueba global (chi-cuadrado de independencia, tabla {len(SUPER)} x 2) ===")
    print(f"chi2 = {chi2:.3f}   gl = {dof}   p = {p_glob:.4f}")
    minima = esperadas.min()
    if minima < 5:
        print(f"AVISO: la frecuencia esperada mas baja es {minima:.1f} (< 5). El chi-cuadrado queda "
              "al limite de su supuesto; la prueba de permutaciones de abajo es la que manda.")
    print(f"p por permutaciones ({N_PERM:,} repartos simulados, semilla {SEMILLA}) = {p_perm:.4f}")

    print(f"\n=== Comparaciones por pares contra {REFERENCIA} (Fisher exacta + BH) ===")
    print(pares.to_string(index=False))

    pares.to_csv(salida, sep="\t", index=False)
    with open(salida.replace(".tsv", "_global.txt"), "w") as fh:
        fh.write(f"chi2\t{chi2:.6f}\ngl\t{dof}\np_chi2\t{p_glob:.6f}\n"
                 f"p_permutaciones\t{p_perm:.6f}\nn_permutaciones\t{N_PERM}\nsemilla\t{SEMILLA}\n"
                 f"esperada_minima\t{minima:.3f}\n")
    print(f"\n>> Escrito: {salida} y {salida.replace('.tsv', '_global.txt')}")
