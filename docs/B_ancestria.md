# B. Análisis de ancestría genética: PCA y ADMIXTURE en 2.504 individuos de 1000 Genomas (chr13 + chr17)

## Métodos (paso a paso)

1. **Individuos:** 2.504 no emparentados (26 poblaciones, 5 superpoblaciones).
2. **SNPs:** bialélicos, FILTER = PASS, MAF ≥ 5 %, faltantes ≤ 2 %. Sin filtro de Hardy-Weinberg sobre el conjunto total, porque la estructura poblacional lo rompe por sí misma (efecto Wahlund) y el filtro eliminaría justo los SNPs informativos.
3. **Región excluida:** inversión 17q21.31 (chr17:45,4–46,8 Mb), una región de LD de largo alcance que puede dominar un componente principal sin reflejar ancestría global.
4. **Poda por LD** por cromosoma con PLINK 2 (`--indep-pairwise 50 5 0.2`) y unión de chr13 + chr17.
5. **PCA** con PLINK 2 (10 componentes).
6. **ADMIXTURE 1.3.0** con K = 2 a 7, validación cruzada de 5 particiones (`--cv`) y semilla fija, sobre una submuestra aleatoria de 60.000 SNPs podados (semilla fija). Este número basta para la estructura continental y hace viable la corrida en un portátil.
7. **Colores alineados entre valores de K:** cada componente hereda el color del componente del K anterior con el que más se correlaciona (asignación húngara), así un mismo color representa la misma ancestría en todas las filas. Cada componente se nombra según la superpoblación y la población donde es máximo.

## PCA

**Figura B1a – PCA de las cinco superpoblaciones** (`results/figuras/B1_pca_superpoblaciones.png`). **Figura B1b – Poblaciones AMR sobre el resto** (`B1_pca_amr.png`).

| Componente | % de la varianza (entre los 10 PCs calculados) | Qué separa |
|---|---:|---|
| PC1 | 60,3 % | África vs. resto del mundo |
| PC2 | 22,6 % | Asia oriental vs. Europa (Asia del sur y AMR quedan entre ambos) |
| PC3 | 6,7 % | Asia del sur |
| PC4 | 5,0 % | Eje indígena americano (AMR) |
| PC5–PC10 | ≤ 1,1 % cada uno | Estructura fina (ruido a esta escala) |

Nota: los porcentajes son relativos a los 10 PCs calculados, no a la varianza genética total. Sirven para comparar la importancia relativa de los ejes.

**Lectura:**
- Las cinco superpoblaciones forman grupos compactos, salvo **AMR**, que se dispersa entre Europa, África y el extremo indígena americano. Es el patrón esperado en poblaciones mezcladas recientemente.
- **ACB y ASW** (afrodescendientes en el Caribe y EE. UU.) forman un "puente" entre África y Europa sobre PC1, por su mezcla europea.
- **En AMR**, la posición de cada persona refleja su mezcla:
  - **PUR** y **CLM** se extienden desde Europa hacia África sobre PC1. Algunos individuos de PUR tienen ancestría africana alta.
  - **PEL** se extiende hacia el extremo negativo de PC4 y hacia Asia oriental en PC2, por la afinidad genética indígena americana–asiática.
  - **MXL** queda en posición intermedia.
- La mediana de PC4 ordena las poblaciones AMR igual que la ancestría indígena estimada por ADMIXTURE:

| | PEL | MXL | CLM | PUR |
|---|---:|---:|---:|---:|
| Mediana PC4 | −0,072 | −0,033 | −0,009 | −0,001 |
| Ancestría indígena (ADMIXTURE K=5) | 77,9 % | 47,1 % | 24,6 % | 14,9 % |

Los dos métodos coinciden, lo que da confianza en el resultado.

## ADMIXTURE

**Figura B2 – Error de validación cruzada por K** (`B2_admixture_cv.png`).

| K | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| Error CV | 0,4823 | 0,4624 | 0,4572 | **0,4535** | 0,4533 | 0,4529 |

El error cae fuerte hasta K=5 y luego se aplana. El mínimo formal está en K=7, pero la mejora de K=5 a K=7 es de 0,0006 (0,1 %), y K=7 es el borde del rango explorado.

**Decisión:** se reporta **K=5 como modelo principal**, porque es el "codo" de la curva y el modelo más simple que captura las cinco grandes ancestrías, incluida la indígena americana. **K=6–7 se usan como estructura fina.** Para afirmar un K óptimo habría que explorar K > 7 con varias semillas.

**Figura B3 – Proporciones de ancestría por individuo** (`B3_admixture_barras.png`). Los colores están alineados entre valores de K y la leyenda nombra cada componente según la superpoblación y la población donde es máximo.

Colores (nombres según K=7): **azul** = africana occidental (máx. en GWD) · **verde oscuro** = africana de Nigeria/África oriental (máx. en LWK) · **naranja** = europea del sur (máx. en TSI) · **violeta** = europea del norte (máx. en FIN) · **verde azulado** = Asia oriental · **amarillo** = Asia del sur · **rosado** = indígena americana (máx. en PEL). En K bajos, un mismo color agrupa ancestrías que se separan después. Por ejemplo, en K=2 el naranja representa a toda la población no africana.

| K | Nuevo componente | Observación (medias por población) |
|---|---|---|
| 2 | África vs. resto | Fracción no africana: ACB 11,9 %, ASW 22,5 %. Fracción africana: PUR 16,2 %, CLM 9,6 %. |
| 3 | Asia oriental | Asia del sur aparece como mezcla (≈ 60–70 % "europeo"). PEL (51 %) y MXL (32 %) muestran un componente "asiático" que en realidad es afinidad indígena americana, porque aún no hay un componente indígena propio. FIN: 7 % oriental. |
| 4 | Asia del sur | Los cinco grupos SAS se separan (74–90 %). |
| 5 | **Indígena americano** | PEL 77,9 % · MXL 47,1 % · CLM 24,6 % · PUR 14,9 %. |
| 6 | División dentro de África | África occidental (GWD 91,6 %, MSL 76,0 %) vs. Nigeria/África oriental (LWK 85,5 %, ESN 59,5 %, YRI 54,8 %). YRI y ESN quedan repartidas casi a la mitad: es un gradiente geográfico oeste–este, no dos poblaciones discretas. |
| 7 | División dentro de Europa | Sur (TSI 82,4 %, IBS 77,3 %) vs. norte (FIN 84,2 %). CEU y GBR quedan ~55/43. |

**Tabla B3. Ancestría media de las poblaciones AMR (K=5)**

| Población | Europea | Indígena americana | Africana | Asia sur + oriental |
|---|---:|---:|---:|---:|
| **CLM (Medellín)** | **64,5 %** | **24,6 %** | **8,4 %** | 2,5 % |
| MXL | 44,9 % | 47,1 % | 3,9 % | 4,1 % |
| PEL | 17,6 % | 77,9 % | 2,5 % | 2,0 % |
| PUR | 67,9 % | 14,9 % | 14,6 % | 2,6 % |

**Lectura para CLM:**
- En promedio, las personas de Medellín tienen **~65 % de ancestría europea, ~25 % indígena americana y ~8 % africana**. Es un patrón de mezcla de tres vías consistente con la historia de la región.
- La proporción **varía bastante entre individuos**: las barras de CLM no son uniformes y en el PCA los individuos se dispersan hacia Europa o hacia África.
- Con K=7, de ese 64 % europeo, **54 puntos corresponden al componente del sur de Europa** (IBS/TSI) y 10 al del norte. Es consistente con el origen ibérico de la colonización.
- Los residuos de "Asia" (< 3 %) en AMR probablemente reflejan la afinidad indígena–asiática y ruido del modelo, no ancestría asiática reciente.

## Cómo la ancestría puede afectar resultados clínicos y genómicos

1. **Frecuencias de referencia:** una variante rara en europeos puede ser común en población indígena o africana, y viceversa. Clasificar patogenicidad con frecuencias europeas puede producir **falsos positivos** (variantes benignas frecuentes en ancestrías poco representadas, que se reportan como VUS o patogénicas) y **falsos negativos**.
2. **Puntajes de riesgo poligénico (PRS):** la mayoría se entrenaron en europeos; su precisión cae ~2 veces en latinoamericanos y ~4,5 veces en africanos (Martin et al., 2019). En Colombia, además, la precisión puede variar **entre individuos** según su proporción de cada ancestría.
3. **Farmacogenómica:** varios alelos de función reducida tienen frecuencias muy distintas entre ancestrías, por ejemplo en CYP2C19, CYP3A5, NUDT15 y UGT1A1. La ancestría individual ayuda a anticipar fenotipos, pero **no reemplaza el genotipo**.
4. **Estudios de asociación y proteómica:** si no se ajusta por ancestría (con PCs o proporciones de ADMIXTURE como covariables), la estructura poblacional produce **asociaciones espurias**.
5. **Colombia no es CLM:** la ancestría varía mucho por región (andina, caribe, pacífica, amazónica). Un programa poblacional requiere referencias locales.

## Limitaciones

- Solo 2 cromosomas y 60.000 SNPs. Suficiente para estructura continental, no para estructura fina ni ancestría local.
- Una sola corrida por K. Lo ideal son varias semillas y alinearlas (por ejemplo, con pong o CLUMPAK).
- 1000 Genomas no tiene una población de referencia indígena americana no mezclada. El componente indígena se infiere de PEL y MXL.
- Los componentes de ADMIXTURE son abstracciones estadísticas, no poblaciones históricas reales.

## Referencias

- Alexander DH, Novembre J, Lange K. Fast model-based estimation of ancestry in unrelated individuals. *Genome Research* 2009;19:1655–1664.
- Chang CC et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* 2015;4:7.
- Byrska-Bishop M et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 2022;185:3426–3440.
- Martin AR et al. Clinical use of current polygenic risk scores may exacerbate health disparities. *Nature Genetics* 2019;51:584–591.
