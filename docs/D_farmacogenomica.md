# D. Módulo Farmacogenómica

## D1. Métodos

- **Datos:** regiones de los 9 genes (coordenadas oficiales de PharmCAT, GRCh38) extraídas del VCF 30x de 1000 Genomas, en las 2.504 personas no emparentadas.
- **Herramienta:** PharmCAT 3.4.0 (preprocesador + *named allele matcher* + fenotipador). Asigna alelos estrella según las definiciones de PharmVar/CPIC y fenotipos según las tablas vigentes de CPIC.
- **Opción clave:** `--absent-to-ref`. El VCF solo lista posiciones con alguna variante, así que una posición ausente se interpreta como referencia. Limitación: una posición con baja cobertura también estaría ausente.
- **VKORC1:** PharmCAT no asigna fenotipo. Se usa el genotipo de rs9923231 (−1639G>A; C>T en GRCh38): portador T = mayor sensibilidad a warfarina.
- **Estadística:** para cada gen se define el **fenotipo accionable según su dirección clínica** (tabla D4). Se compara CLM contra EUR, AFR y EAS con la prueba exacta de Fisher, IC 95 % de Wilson y corrección de Benjamini-Hochberg (FDR) sobre las 30 comparaciones. Los indeterminados se excluyen.
- **Entregables:** `D1_diplotipos_AMR.tsv` (diplotipo y fenotipo por individuo AMR y gen), `D2_fenotipos_por_grupo.tsv`, `D3_diplotipos_frecuentes_AMR.tsv`, `D4_pruebas_CLM_vs_otros.tsv` y la figura `D1_fenotipos_CLM_vs_otros.png`.

**Decisión metodológica.** Una primera versión comparaba "fenotipo no normal vs. normal". Eso mezcla efectos opuestos: en CYP2C19, "no normal" suma metabolizadores pobres y ultrarrápidos. En CYP3A5 lo clínicamente accionable es ser *expresor* (fenotipo normal o intermedio), porque requiere más tacrolimus. Con la definición direccional cambia la conclusión para CYP2C19: CLM **no** difiere de EUR en función disminuida, pero sí tiene **menos** metabolizadores rápidos/ultrarrápidos.

## D2. Resultados

**Tabla D2. % de individuos con fenotipo accionable** (excluye indeterminados; n: CLM 94, MXL 64, PEL 85, PUR 104, EUR 503, AFR 661, EAS 504, SAS 489)

| Gen | Fenotipo accionable | CLM | MXL | PEL | PUR | EUR | AFR | EAS |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CYP2C19 | función disminuida (PM+IM) | 21,5 | 25,0 | 11,8 | 26,0 | 28,6 | 38,0 | 62,1 |
| CYP2C19 | función aumentada (RM+UM) | 19,4 | 17,2 | 7,1 | 29,8 | 33,6 | 32,1 | 1,4 |
| CYP2C9 | función disminuida (PM+IM) | 35,1 | 25,4 | 7,1 | 36,3 | 36,6 | 21,4 | 7,9 |
| CYP3A5 | expresor (NM+IM) | 30,9 | 39,1 | 18,8 | 36,5 | 10,3 | 78,5 | 49,4 |
| DPYD | función disminuida (PM+IM) | 1,6 | 2,3 | 0,0 | 0,0 | 2,9 | 0,2 | 0,0 |
| TPMT | función disminuida (PM+IM) | 6,4 | 9,5 | 14,1 | 17,3 | 7,6 | 13,6 | 4,4 |
| NUDT15 | función disminuida (PM+IM) | 8,5 | 12,5 | 22,4 | 2,9 | 1,0 | 0,6 | 20,4 |
| SLCO1B1 | función disminuida o pobre | 33,7 | 15,6 | 26,5 | 25,7 | 30,2 | 19,6 | 24,3 |
| UGT1A1 | función disminuida (PM+IM) | 58,5 | 65,6 | 72,9 | 61,5 | 52,1 | 75,6 | 43,1 |
| VKORC1 | portador rs9923231-T | 64,9 | 73,4 | 62,4 | 61,5 | 61,8 | 10,4 | 98,2 |

**Validación con la literatura.** Las frecuencias de las poblaciones de referencia coinciden con lo esperado: CYP3A5 no expresor (\*3/\*3) en 90 % de EUR; CYP2C19 metabolizador pobre en 12 % de EAS; UGT1A1 pobre (\*28/\*28) en 23 % de AFR; VKORC1-T en 98 % de EAS; NUDT15 disminuido en 20 % de EAS. Esto respalda la calidad del flujo de trabajo.

**Tabla D4. CLM frente a otras superpoblaciones (q < 0,05 tras FDR)**

| Gen | Fenotipo | CLM % (IC 95 %) | vs EUR | vs AFR | vs EAS |
|---|---|---|---|---|---|
| CYP2C19 | disminuida | 21,5 (14,4–30,9) | = | **menor** (38,0) | **menor** (62,1) |
| CYP2C19 | aumentada | 19,4 (12,6–28,5) | **menor** (33,6) | **menor** (32,1) | **mayor** (1,4) |
| CYP2C9 | disminuida | 35,1 (26,2–45,2) | = | **mayor** (21,4) | **mayor** (7,9) |
| CYP3A5 | expresor | 30,9 (22,4–40,8) | **mayor** (10,3) | **menor** (78,5) | **menor** (49,4) |
| NUDT15 | disminuida | 8,5 (4,4–15,9) | **mayor** (1,0) | **mayor** (0,6) | **menor** (20,4) |
| SLCO1B1 | disminuida | 33,7 (24,9–43,8) | = | **mayor** (19,6) | = |
| UGT1A1 | disminuida | 58,5 (48,4–67,9) | = | **menor** (75,6) | **mayor** (43,1) |
| VKORC1 | portador T | 64,9 (54,8–73,8) | = | **mayor** (10,4) | **menor** (98,2) |
| DPYD, TPMT | disminuida | — | = | = | = |

"=" significa sin diferencia significativa. Entre paréntesis, el % del grupo de comparación. Tabla completa con OR, p y q en `D4_pruebas_CLM_vs_otros.tsv`.

## D3. Discusión e implicaciones clínicas

1. **CLM no es "europeo" en farmacogenómica.** Se parece a EUR en CYP2C9, SLCO1B1, UGT1A1 y VKORC1, pero difiere en genes donde pesa la ancestría indígena o africana:
   - **NUDT15:** 8,5 % de CLM tiene función disminuida, contra 1 % en EUR (OR 9,3). En PEL llega a 22 %. Un programa que solo genotipe TPMT, como es habitual en guías basadas en europeos, dejaría sin detectar a la mayoría de los colombianos con riesgo de mielotoxicidad por tiopurinas (leucemia linfoblástica aguda pediátrica, enfermedad inflamatoria intestinal).
   - **CYP3A5:** 31 % de CLM son expresores, contra 10 % en EUR. En trasplante, esos pacientes necesitan dosis iniciales más altas de tacrolimus para llegar al rango terapéutico.
   - **CYP2C19:** CLM tiene menos metabolizadores rápidos/ultrarrápidos (\*17) que EUR, con igual proporción de función disminuida. Eso afecta a clopidogrel (~1 de cada 5 con activación reducida), inhibidores de bomba de protones y antidepresivos ISRS.
2. **Alto impacto poblacional:**
   - **SLCO1B1:** 1 de cada 3 en CLM tiene función disminuida o pobre → riesgo de miopatía con simvastatina, uno de los fármacos más usados.
   - **Warfarina:** 65 % de CLM porta VKORC1-T y 35 % tiene CYP2C9 disminuido. Un algoritmo de dosis que integre ambos genes es relevante.
   - **UGT1A1:** 15 % es metabolizador pobre → toxicidad por irinotecán.
3. **La heterogeneidad dentro de AMR es grande.** Por ejemplo, NUDT15 va de 2,9 % en PUR a 22,4 % en PEL, y CYP3A5 expresor de 18,8 % a 39,1 %. Sigue el mismo orden que la ancestría indígena y africana de la sección B, así que las frecuencias de una población latinoamericana no se pueden extrapolar a otra.
4. **DPYD:** los portadores de variantes de función disminuida son raros (CLM 1,6 %, 1 de 94), coherente con lo reportado. Con n = 94 el poder estadístico no alcanza para compararlos.

**Control de calidad: los indeterminados de DPYD.** El 26,4 % de los individuos (661/2.504) queda sin fenotipo. La causa no es falta de datos: PharmCAT marca *Indeterminate* cuando un haplotipo lleva una **combinación de variantes que no corresponde a ningún alelo con nombre definido** (p. ej. `[c.85T>C (*9A) + c.496A>G]`), aunque cada variante por separado sí esté clasificada por CPIC.

Se reconstruyó el puntaje de actividad de cada haplotipo indeterminado a partir de la función CPIC de sus variantes (`scripts/utils/dpyd_indeterminados.py`). De los 661 indeterminados, **574 solo contienen variantes de función normal** (equivalen a metabolizadores normales) y **87 contienen al menos una variante de función disminuida**, casi siempre **HapB3** (`c.1129-5923C>G, c.1236G>A`, valor 0,5) o **c.557A>G** (0,5). Esas 87 personas son clínicamente **metabolizadoras intermedias** y deberían recibir la mitad de la dosis inicial de fluoropirimidinas, pero el reporte automático no las señala.

**Tabla D5. Portadores de DPYD de función disminuida: reportados por PharmCAT vs. reconstruidos**

| Grupo | n | Indeterminados % | PharmCAT n (%) | Reconstruido n (%) | Ocultos | Factor |
|---|---:|---:|---:|---:|---:|---:|
| CLM | 94 | 33,0 | 1 (1,1) | 3 (3,2) | 2 | 3,0× |
| MXL | 64 | 31,2 | 1 (1,6) | 3 (4,7) | 2 | 3,0× |
| PEL | 85 | 25,9 | 0 (0,0) | 2 (2,4) | 2 | — |
| PUR | 104 | 26,9 | 0 (0,0) | 2 (1,9) | 2 | — |
| EUR | 503 | 30,6 | 10 (2,0) | 36 (7,2) | 26 | 3,6× |
| AFR | 661 | 35,6 | 1 (0,2) | 34 (5,1) | 33 | **34×** |
| EAS | 504 | 13,7 | 0 (0,0) | 0 (0,0) | 0 | — |
| SAS | 489 | 20,9 | 12 (2,5) | 32 (6,5) | 20 | 2,7× |

Ambas columnas usan el mismo denominador (todos los individuos), porque un indeterminado es, en la práctica clínica, un paciente a quien el reporte **no le señala riesgo**.

**Lectura.** La subdetección afecta a todas las poblaciones, pero es **desigual**: en EUR el reporte capta 10 de 36 portadores reales (3,6×), mientras que en **AFR capta 1 de 34 (34×)**. La razón es `c.557A>G`, una variante de función disminuida frecuente en poblaciones de ancestría africana que aquí aparece casi siempre acompañada de `c.85T>C (*9A)` en el mismo haplotipo, una combinación sin alelo nombrado. Es decir, **la herramienta rinde peor justamente en las poblaciones menos representadas en los catálogos de alelos**. En CLM el efecto es menor en números absolutos (1 → 3 portadores) pero mantiene la dirección.

**Advertencia metodológica.** Esta reconstrucción es una **hipótesis de trabajo, no un llamado clínico**. PharmCAT es conservador a propósito: un haplotipo que lleva HapB3 *más* otra variante no es un alelo caracterizado, y no puede descartarse que la combinación se comporte distinto. Lo correcto no es reclasificar automáticamente, sino **marcar estos casos para revisión manual**.

**Implicación para SURA.** Un reporte automático de farmacogenómica **no se entrega sin auditar los indeterminados**: aquí esconden hasta el 97 % de los portadores accionables en una superpoblación. Un programa clínico necesita (a) una regla explícita de manejo de indeterminados, (b) revisión por un profesional antes de emitir el reporte, y (c) medir la tasa de indeterminados **por ancestría** como indicador de equidad del servicio.

## D4. Limitaciones

- n = 94 en CLM: IC amplios y poder limitado para fenotipos raros (DPYD, TPMT pobre, NUDT15 pobre).
- CLM (Medellín) no representa a Colombia (ver sección B).
- `--absent-to-ref` asume referencia en posiciones sin datos.
- Solo se evalúan las variantes que definen alelos conocidos. Las variantes raras o nuevas, más frecuentes en poblaciones poco estudiadas, no se interpretan.
- CYP2D6, uno de los genes farmacogenéticos más relevantes, no se puede tipificar desde este VCF.
- Indeterminados: DPYD 26,4 % en total (13,7–35,6 % según población; ver control de calidad arriba) y SLCO1B1 18,3 % en AFR, probablemente por la misma causa: combinaciones de variantes sin alelo nombrado. Convendría auditar SLCO1B1 con el mismo método.
