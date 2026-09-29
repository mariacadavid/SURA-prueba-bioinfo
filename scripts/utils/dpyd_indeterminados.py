#!/usr/bin/env python
"""¿Qué hay detrás de los llamados 'Indeterminate' de DPYD? (control de calidad del módulo PGx)

PharmCAT marca Indeterminate cuando un haplotipo lleva una COMBINACIÓN de variantes que no
corresponde a ningún alelo con nombre definido, aunque las variantes sí estén clasificadas.
Este script abre cada diplotipo indeterminado, busca variantes con función alterada según CPIC
y separa dos grupos:
  - "sin impacto": todas las variantes del haplotipo son de función normal -> equivalen a NM.
  - "con variante accionable": hay al menos una variante de función disminuida o nula
    -> clínicamente serían metabolizadores intermedios o pobres y PharmCAT no los reporta como tales.

Uso: python scripts/utils/dpyd_indeterminados.py data/pgx/pharmcat data/meta/samples_3202_ped_population.txt
"""
import sys, glob, os, re
import pandas as pd

# Variantes DPYD con función alterada (CPIC, consultado 27-09-2026). Valor de actividad entre paréntesis.
SIN_FUNCION = ["c.1905+1G>A", "c.1679T>G", "c.1898delC", "c.295_298delTCAT", "c.703C>T", "c.1156G>T",
               "c.2983G>T", "c.61C>T", "c.601A>C", "c.632A>G", "c.1024G>A", "c.1057C>T", "c.1475C>T",
               "c.1484A>G", "c.1774C>T", "c.1775G>A", "c.1777G>A", "c.2021G>A", "c.2639G>T",
               "c.2872A>G", "c.2933A>G"]
DISMINUIDA = ["HapB3", "c.1129-5923C>G", "c.2846A>T", "c.1236G>A", "c.557A>G", "c.868A>G",
              "c.1314T>G", "c.2279C>T"]

def clasifica_haplotipo(h):
    if any(v in h for v in SIN_FUNCION): return 0.0
    if any(v in h for v in DISMINUIDA):  return 0.5
    return 1.0

def parte_diplotipo(d):
    """Separa el diplotipo en sus dos haplotipos respetando los corchetes de las combinaciones."""
    nivel, corte = 0, None
    for i, c in enumerate(d):
        if c == "[": nivel += 1
        elif c == "]": nivel -= 1
        elif c == "/" and nivel == 0: corte = i; break
    return (d[:corte], d[corte+1:]) if corte else (d, d)

def main(DIR, META):
    TAB = "results/tablas"; os.makedirs(TAB, exist_ok=True)
    filas = []
    for f in glob.glob(f"{DIR}/*.report.tsv"):
        muestra = os.path.basename(f).split(".")[-3]
        for linea in open(f):
            if linea.startswith("DPYD\t"):
                c = linea.rstrip("\n").split("\t")
                dip, fen = c[1], c[2]
                h1, h2 = parte_diplotipo(dip)
                a1, a2 = clasifica_haplotipo(h1), clasifica_haplotipo(h2)
                filas.append({"IID": muestra, "diplotipo": dip, "fenotipo_PharmCAT": fen,
                              "AS_reconstruido": a1 + a2})
                break
    d = pd.DataFrame(filas)

    ped = pd.read_csv(META, sep=r"\s+")
    ped.columns = [re.sub(r"[^a-z]", "", x.lower()) for x in ped.columns]
    idc = next(x for x in ped.columns if "sample" in x); supc = next(x for x in ped.columns if "superpop" in x)
    ped = ped.rename(columns={idc: "IID", "population": "poblacion", supc: "superpoblacion"})
    d = d.merge(ped[["IID", "poblacion", "superpoblacion"]], on="IID", how="left")
    d["grupo"] = d.superpoblacion.where(d.superpoblacion != "AMR", d.poblacion)

    def fen(a):  # rangos CPIC para DPYD
        return "Pobre" if a <= 0.5 else ("Intermedio" if a <= 1.5 else "Normal")
    d["fenotipo_rescatado"] = d.AS_reconstruido.map(fen)

    ind = d[d.fenotipo_PharmCAT.str.contains("Indeterminate", na=False)].copy()
    print(f"Individuos con DPYD: {len(d)} · indeterminados: {len(ind)} ({100*len(ind)/len(d):.1f} %)")
    print("\n-- Fenotipo reconstruido de los indeterminados:")
    print(ind.fenotipo_rescatado.value_counts().to_string())
    print("\n-- Indeterminados que SÍ llevan una variante de función alterada:")
    oculta = ind[ind.fenotipo_rescatado != "Normal"]
    print(f"   {len(oculta)} individuos")
    if len(oculta):
        print(oculta.groupby(["grupo", "fenotipo_rescatado"]).size().to_string())
        print("\n-- Diplotipos involucrados:")
        print(oculta.diplotipo.value_counts().to_string())

    # Tabla final: % con función disminuida antes y después del rescate
    G = ["CLM", "MXL", "PEL", "PUR", "EUR", "AFR", "EAS", "SAS"]
    res = []
    for g in G:
        s = d[d.grupo == g]
        val = s[~s.fenotipo_PharmCAT.str.contains("Indeterminate", na=False)]
        res.append({"grupo": g, "n": len(s),
                    "indeterminados_%": round(100 * (len(s) - len(val)) / len(s), 1),
                    "PharmCAT_disminuida_%": round(100 * val.fenotipo_PharmCAT.str.contains("Intermediate|Poor").sum() / len(val), 1),
                    "rescatado_disminuida_%": round(100 * (s.fenotipo_rescatado != "Normal").sum() / len(s), 1)})
    r = pd.DataFrame(res)
    r.to_csv(f"{TAB}/D5_dpyd_indeterminados.tsv", sep="\t", index=False)
    d.to_csv(f"{TAB}/D5_dpyd_por_individuo.tsv", sep="\t", index=False)
    print("\n-- Tabla D5 (results/tablas/D5_dpyd_indeterminados.tsv):")
    print(r.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
