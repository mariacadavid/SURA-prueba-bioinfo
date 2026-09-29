# C. Análisis de variantes en genes asociados a enfermedad

## C1. Métodos

**Selección de genes.** En vez de elegir genes a criterio propio, se aplicó un criterio publicado: la lista **ACMG SF** (*Secondary Findings*), que reúne los genes en los que un hallazgo tiene **utilidad clínica accionable**. Restringida a los cromosomas del ejercicio, quedan 6 genes:

| Cromosoma | Gen | Herencia | Enfermedad |
|---|---|---|---|
| 13 | BRCA2 | Dominante | Cáncer de mama y ovario hereditario (HBOC) |
| 13 | RB1 | Dominante | Retinoblastoma |
| 13 | ATP7B | Recesiva | Enfermedad de Wilson |
| 17 | BRCA1 | Dominante | HBOC |
| 17 | TP53 | Dominante | Síndrome de Li-Fraumeni |
| 17 | GAA | Recesiva | Enfermedad de Pompe |

**Datos y anotación.** Regiones de los 6 genes (coordenadas de Ensembl, GRCh38, con 2 kb de margen), en las 2.504 personas no emparentadas, variantes bialélicas con FILTER = PASS. La **clasificación clínica**, la **enfermedad**, el **nivel de revisión** y el **impacto funcional** se tomaron del VCF de ClinVar (campos `CLNSIG`, `CLNDN`, `CLNREVSTAT` y `MC`); la versión descargada queda registrada en `data/clinicas/clinvar_version.txt`.

**Criterios de filtrado**, aplicados en cascada y declarados de antemano:

1. Clasificación **Pathogenic** o **Likely pathogenic**, excluyendo las clasificaciones en conflicto.
2. **Nivel de revisión ≥ 2 estrellas** (criterios aportados por varios laboratorios sin conflicto, panel de expertos o guía de práctica clínica).
3. **Frecuencia global compatible con la prevalencia de la enfermedad**: ≤ 0,1 % para genes dominantes y ≤ 2 % para recesivos (frecuencia de portadores). Una variante clasificada como patogénica que sea más frecuente que eso contradice la epidemiología de la enfermedad y es sospechosa de estar mal clasificada.

Se usó la frecuencia **global** y no la máxima por superpoblación porque, con 347 a 661 personas por superpoblación, un solo portador ya produce una frecuencia alélica de ~0,08 %, y el filtro habría descartado variantes raras legítimas. La frecuencia máxima por superpoblación se reporta igual, como señal de alerta adicional.

## C2. Resultados

**Tabla C3. Embudo de filtrado**

| Paso | Variantes |
|---|---:|
| 1. Variantes en los 6 genes con clasificación en ClinVar | **2.104** |
| 2. Clasificadas Pathogenic / Likely pathogenic, sin conflicto | **12** |
| 3. + nivel de revisión ≥ 2 estrellas | 12 |
| 4. + frecuencia compatible con la enfermedad | **12** |

**Resultado destacable: los filtros 3 y 4 no eliminaron ninguna variante.** Las 12 candidatas ya tenían respaldo sólido: **9 de 12 están revisadas por panel de expertos (3 estrellas)** y las 3 restantes tienen 2 estrellas. Ninguna resultó demasiado frecuente para su enfermedad. Esto no era el resultado esperado, y tiene una explicación: BRCA1, BRCA2 y GAA están entre los genes **mejor curados de ClinVar**, con paneles de expertos dedicados que ya depuraron las clasificaciones dudosas. En un panel de genes menos curados, el mismo filtro sí habría recortado.

El contraste del paso 1 al 2 es igualmente informativo: de 2.104 variantes con clasificación, **solo el 0,6 % es patogénica**. El resto son benignas, probablemente benignas o de significado incierto (VUS). Es el patrón esperado en una cohorte de voluntarios sanos.

**Tabla C2. Las 12 variantes patogénicas o probablemente patogénicas**

| Gen | Posición (GRCh38) | Impacto funcional | Clasificación | Estrellas | Portadores | AF global |
|---|---|---|---|---:|---:|---:|
| BRCA2 | chr13:32339596 C>CTA | Truncante | Pathogenic | 3 | 1 | 0,02 % |
| BRCA2 | chr13:32357820 G>GA | Frameshift | Pathogenic | 3 | 1 | 0,02 % |
| ATP7B | chr13:51950132 C>T | Missense | Pathogenic/Likely path. | 2 | 5 | 0,10 % |
| ATP7B | chr13:51974407 G>T | Nonsense | Pathogenic/Likely path. | 2 | 2 | 0,04 % |
| **GAA** | **chr17:80104542 T>G** | **Splicing (intrón 1)** | **Pathogenic** | **3** | **14** | **0,28 %** |
| GAA | chr17:80105873 G>A | Missense | Likely pathogenic | 3 | 1 | 0,02 % |
| GAA | chr17:80107705 C>T | Missense | Likely pathogenic | 3 | 2 | 0,04 % |
| GAA | chr17:80107866 G>A | Missense | Pathogenic | 3 | 1 | 0,02 % |
| GAA | chr17:80112666 G>A | Missense | Pathogenic | 3 | 1 | 0,02 % |
| GAA | chr17:80112929 G>A | Missense | Pathogenic | 3 | 1 | 0,02 % |
| GAA | chr17:80117016 G>C | Missense | Pathogenic | 3 | 2 | 0,04 % |
| GAA | chr17:80118271 C>T | Nonsense | Pathogenic | 3 | 5 | 0,10 % |

**Impacto funcional:** 4 variantes truncantes (nonsense o frameshift), 7 missense y 1 de splicing. **Distribución por gen:** GAA 8 variantes (27 portadores), ATP7B 2 (7 portadores), BRCA2 2 (2 portadores). **BRCA1, TP53 y RB1: ninguna variante patogénica.**

**Hallazgo identificado: `chr17:80104542 T>G` es rs386834236 = GAA c.-32-13T>G**, la variante más frecuente de la enfermedad de Pompe de inicio tardío. Con 14 portadores es la más común de la lista, y su distribución (EUR 0,70 %, AMR 0,58 %, AFR 0,15 %, SAS 0,10 %, EAS 0 %) concuerda con lo descrito: es una variante de origen europeo. Que el flujo de trabajo la recupere de forma independiente, con la frecuencia y el patrón poblacional esperados, es una **validación del método**.

## C3. ¿Se agrupan por ancestría?

**Tabla C4. Portadores de variantes P/LP por superpoblación**

| Superpoblación | Portadores | n | Tasa |
|---|---:|---:|---:|
| EUR | 13 | 503 | **2,58 %** |
| AMR | 6 | 347 | 1,73 % |
| AFR | 10 | 661 | 1,51 % |
| EAS | 4 | 504 | 0,79 % |
| SAS | 3 | 489 | 0,61 % |
| **Total** | **36** | **2.504** | **1,44 %** |

En CLM hubo **1 portador de 94** (1,06 %), de la variante GAA chr17:80117016.

**Prueba estadística.** La comparación global entre las cinco superpoblaciones **no es significativa** (chi-cuadrado = 8,72; 4 grados de libertad; p = 0,068). En las comparaciones de EUR contra cada una, EUR vs. EAS (p = 0,029) y EUR vs. SAS (p = 0,021) son nominalmente significativas, pero **ninguna sobrevive la corrección por comparaciones múltiples** (q = 0,059 en ambas, Benjamini-Hochberg).

**Conclusión: no hay evidencia estadística de agrupamiento por ancestría** con este tamaño muestral. Y aunque la hubiera, habría que interpretarla con cuidado por dos razones:

1. **El resultado está dominado por un solo gen.** 27 de los 36 portadores lo son de GAA, y 14 de una sola variante europea (c.-32-13T>G). El patrón observado refleja la distribución de esa variante, no una propiedad general de la carga de variantes patogénicas.
2. **Sesgo de catálogo.** Solo se pueden detectar variantes que **ya están en ClinVar**, y ClinVar está enriquecido en variantes descritas en poblaciones europeas. Una mayor tasa aparente en EUR es exactamente lo que produciría ese sesgo, aunque la carga real fuera igual en todas las poblaciones. **La ancestría con más hallazgos puede ser simplemente la mejor estudiada.**

## C4. Interpretación crítica

**Se detectaron menos portadores de BRCA1/BRCA2 de los esperados.** Se encontraron 2 portadores en 2.504 personas (1 de cada 1.252). La frecuencia de portadores de variantes patogénicas en BRCA1 o BRCA2 en población general se estima en torno a 1 de cada 300 a 500, lo que daría entre 5 y 8 portadores esperados. Tres explicaciones, en orden de peso:

1. **Solo se detecta lo que ya está catalogado.** La mayoría de las variantes patogénicas de BRCA1 y BRCA2 son truncantes **privadas de una familia**. Una variante nueva, que no esté en ClinVar, es invisible para este método aunque esté en los datos. Para encontrarlas haría falta clasificar de novo con criterios ACMG/AMP, lo que excede el alcance de este ejercicio.
2. **Los datos son un panel de referencia faseado**, no genotipos de grado clínico; los indels raros pueden estar subrepresentados.
3. **1000 Genomas son voluntarios adultos sanos**, con cierta depleción esperable de portadores de condiciones de alta penetrancia, aunque el efecto a estas frecuencias es pequeño.

**El mismo sesgo, en dos direcciones.** En el módulo de farmacogenómica el catálogo incompleto produjo **subdetección en ancestría africana** (sección D3). Aquí produce **subdetección de variantes privadas y una aparente concentración en EUR**. Es el mismo problema visto desde dos ángulos: *lo que no está en el catálogo no existe para el análisis*, y los catálogos se construyeron sobre poblaciones europeas.

**Implicación para un programa poblacional.** Un tamizaje basado únicamente en variantes ya clasificadas subestimará el riesgo, y lo hará de forma desigual entre ancestrías. Un programa de SURA necesitaría capacidad de **clasificación de novo con criterios ACMG/AMP** y participación en el depósito de variantes a ClinVar, no solo consulta pasiva.

## C5. Limitaciones

- **Solo variantes ya clasificadas en ClinVar.** No se clasificaron variantes de novo; las variantes nuevas o privadas quedan fuera (ver C4).
- **No se evaluaron variantes estructurales ni deleciones de exones**, que son un mecanismo relevante en BRCA1 y BRCA2.
- **El impacto funcional proviene del campo `MC` de ClinVar.** Cuando una variante tiene varias consecuencias anotadas, el script reporta la más severa según un orden fijo de prioridad.
- **Poder estadístico insuficiente** para la pregunta de agrupamiento por ancestría: con 36 portadores repartidos en 5 superpoblaciones y un solo gen dominando el conteo, solo se detectarían diferencias muy grandes.
- **Solo 6 genes en 2 cromosomas.** No es un tamizaje de hallazgos secundarios completo, que requeriría los 81 genes de la lista ACMG SF.
- Las clasificaciones de ClinVar **cambian con el tiempo**; el resultado corresponde a la versión descargada, que queda registrada en el repositorio.
