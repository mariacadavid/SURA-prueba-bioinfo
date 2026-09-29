# A. Descripción del dataset

## Fuente y decisiones

- **Datos:** 1000 Genomas, secuenciación de alta cobertura (30x, Illumina), ensamblaje GRCh38, panel faseado de SNV, indels y variantes estructurales (NYGC, 20220422; Byrska-Bishop et al., *Cell* 2022). Cromosomas 13 y 17.
- **Nota sobre la fuente:** la tabla de la prueba describe datos 30x GRCh38, pero el enlace apunta a `release/20130502` (Fase 3, baja cobertura, GRCh37). Se eligió la versión 30x GRCh38 porque es la que describe el enunciado y el ensamblaje de ClinVar, gnomAD y PharmCAT. De `release/20130502` se tomó el panel de las 2.504 muestras de la Fase 3 para definir el conjunto no emparentado.
- **Individuos:** de 3.202 muestras se analizaron las **2.504 no emparentadas**. Las 698 restantes son familiares agregados en la versión 30x (602 tríos); incluir parientes infla frecuencias y sesga PCA y ADMIXTURE.
- **Filtros:** FILTER = PASS; SNPs bialélicos A/C/G/T para los análisis de frecuencia y ancestría. Indels y variantes estructurales se reportan pero no se usan en PCA/ADMIXTURE.
- **Categorías de frecuencia (MAF, alelo menor en las 2.504):** común ≥ 5 %, baja frecuencia 1–5 %, rara < 1 %.

## Tabla A1. Variantes por cromosoma

| | chr13 | chr17 |
|---|---:|---:|
| Variantes totales en el VCF | 2.509.179 | 2.073.624 |
| SNPs | 2.186.334 | 1.767.777 |
| Indels | 319.428 | 302.526 |
| Otras (variantes estructurales) | 3.417 | 3.321 |
| SNPs por Mb | ~19.100 | ~21.200 |
| SNPs monomórficos en las 2.504 | 28.664 | 24.143 |
| **SNPs polimórficos en las 2.504** | **2.157.670** | **1.743.634** |
| Comunes (MAF ≥ 5 %) | 253.510 (11,7 %) | 190.239 (10,9 %) |
| Baja frecuencia (1–5 %) | 200.679 (9,3 %) | 153.959 (8,8 %) |
| Raras (< 1 %) | 1.703.481 (79,0 %) | 1.399.436 (80,3 %) |
| Singletons (una sola copia) | 385.448 (17,9 %) | 323.501 (18,6 %) |

Porcentajes sobre los SNPs polimórficos.

**Lectura:**
- **~80 % de los SNPs son raros.** Es el patrón esperado en poblaciones humanas que crecieron rápido en tamaño. La mayoría de variantes raras son privadas de una población o de un continente, y solo la secuenciación (no los arrays) las captura.
- **Solo ~11 % son comunes**, y son las que se usan para ancestría.
- **chr17 es más denso en SNPs por Mb que chr13.** chr17 es rico en genes; chr13 es acrocéntrico y su brazo corto casi no se secuencia, lo que baja su densidad aparente.
- **El VCF no tiene sitios multialélicos** porque el panel viene "normalizado": cada alelo alternativo está en su propia fila. Los SNPs monomórficos (~25 mil por cromosoma) son variantes que solo aparecen en los 698 familiares excluidos.

## Tabla A2. Datos faltantes

| Individuos analizados | Proporción máxima de genotipos faltantes | Individuos con > 2 % faltantes |
|---:|---:|---:|
| 2.504 | 0 | 0 |

**Lectura:** no hay genotipos faltantes, pero esto **no** significa que la secuenciación sea perfecta. El panel fue **faseado estadísticamente**, y ese proceso imputa los genotipos faltantes. En datos propios de SURA (sin fasear) se esperaría missingness real, y se aplicarían filtros por individuo (> 2–5 %) y por variante. Tampoco hay individuos ausentes: las 2.504 muestras de la Fase 3 están en la versión 30x.

## Tabla A3. Individuos por continente (superpoblación)

| Superpoblación | Muestras 30x | No emparentados | Proporción |
|---|---:|---:|---:|
| AFR – África | 893 | 661 | 26,4 % |
| EAS – Asia oriental | 585 | 504 | 20,1 % |
| EUR – Europa | 633 | 503 | 20,1 % |
| SAS – Asia del sur | 601 | 489 | 19,5 % |
| AMR – Américas (mezcladas) | 490 | 347 | 13,9 % |
| **Total** | **3.202** | **2.504** | 100 % |

Detalle AMR: CLM (Medellín) 94 · MXL 64 · PEL 85 · PUR 104. Detalle completo por población en `results/tablas/A4_individuos_por_poblacion.tsv`.

**Lectura:** AMR es la superpoblación más pequeña (14 %) y CLM tiene solo 94 personas de un área metropolitana. Esto limita extrapolar frecuencias a Colombia, un país con alta heterogeneidad regional en ancestría (andina, caribe, pacífica, amazónica).
