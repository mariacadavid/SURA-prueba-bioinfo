# D. Panel de farmacogenómica preventiva para afiliados de SURA

## 1. Objetivo y alcance

**Objetivo:** genotipar **una sola vez en la vida** a cada afiliado en un conjunto acotado de genes con evidencia clínica accionable, guardar el resultado en la historia clínica y **hacerlo visible al médico en el momento de prescribir**.

**Decisión de diseño: panel mínimo, no panel máximo.** Se propone un MVP de **7 genes** en lugar de los 12–14 de los paneles europeos. El criterio no es técnico sino de implementación: el cuello de botella de un programa preventivo no es genotipar, sino **integrarlo al flujo de prescripción y demostrar impacto**. Un panel pequeño permite validar esa cadena completa antes de escalar. La ruta de expansión está definida en la sección 6.

**Criterios de inclusión de un gen (los cuatro deben cumplirse):**

1. **Evidencia:** par gen–fármaco de **nivel A de CPIC** (existe recomendación clínica accionable), respaldado por etiqueta de FDA o EMA cuando aplica.
2. **Volumen:** el fármaco se prescribe de forma relevante en la población afiliada de SURA (dato interno, paso obligatorio antes de cerrar el panel).
3. **Impacto:** el fenotipo accionable tiene frecuencia apreciable en población colombiana.
4. **Viabilidad:** el genotipo se puede determinar de forma confiable con la tecnología elegida.

## 2. Genes del panel mínimo

| Gen | Fármacos (CPIC nivel A) | Acción clínica | % accionable en CLM* |
|---|---|---|---:|
| **CYP2C19** | Clopidogrel; omeprazol, pantoprazol, lansoprazol; citalopram, escitalopram, sertralina, amitriptilina; voriconazol | Cambiar antiagregante en metabolizador pobre o intermedio; ajustar dosis de IBP e ISRS | 21,5 % función disminuida |
| **CYP2C9** | Warfarina; fenitoína; celecoxib, ibuprofeno, meloxicam, piroxicam | Dosis inicial de warfarina; reducir dosis de AINE y fenitoína | 35,1 % |
| **VKORC1** | Warfarina | Dosis inicial (junto con CYP2C9) | 64,9 % portador rs9923231-T |
| **SLCO1B1** | **Todas las estatinas** (simvastatina, atorvastatina, rosuvastatina, pravastatina, lovastatina, fluvastatina, pitavastatina) | Limitar dosis o cambiar de estatina para evitar miopatía | 33,7 % función disminuida |
| **DPYD** | Capecitabina, 5-fluorouracilo | **Reducir dosis inicial 50 %**; evita toxicidad grave o fatal | 1,6 % (3,2 % al reconstruir indeterminados) |
| **TPMT** | Azatioprina, mercaptopurina, tioguanina | Reducir dosis para evitar mielosupresión | 6,4 % |
| **NUDT15** | Azatioprina, mercaptopurina, tioguanina | Igual que TPMT | **8,5 %** |

\* Calculado en este ejercicio sobre 94 individuos CLM de 1000 Genomas (sección D2). Los intervalos de confianza son amplios; ver limitaciones.

**Por qué estos siete y no otros.** El panel cubre cuatro escenarios clínicos de alto volumen o alto riesgo en una aseguradora con EPS e IPS:

- **Cardiovascular:** clopidogrel (CYP2C19), warfarina (CYP2C9 + VKORC1) y estatinas (SLCO1B1). Es el grupo de mayor volumen de prescripción.
- **Oncología:** DPYD. Bajo número de pacientes, pero **la consecuencia es muerte evitable**: la recomendación CPIC es de fuerza *fuerte* y la EMA recomienda evaluar el déficit de DPD antes de iniciar fluoropirimidinas.
- **Inmunosupresión y hematooncología pediátrica:** TPMT + NUDT15.
- **Salud mental y gastroenterología:** CYP2C19 cubre ISRS e inhibidores de bomba de protones sin costo adicional, porque ya está en el panel.

**El argumento decisivo de este análisis: NUDT15 no es opcional en Colombia.** En los datos de CLM, el 8,5 % tiene función disminuida (IC 95 % 4,4–15,9) frente al **1,0 % en europeos**. Un panel copiado de guías europeas que incluyera solo TPMT dejaría sin detectar a la mayoría de los colombianos en riesgo de mielotoxicidad por tiopurinas. **Incluir TPMT sin NUDT15 sería replicar en Colombia un sesgo de la evidencia europea.**

**Nota de costo marginal.** CYP4F2 (warfarina, nivel A) es un solo SNP (rs2108622) y puede añadirse sin costo relevante si la tecnología ya cubre la región. Su efecto sobre la dosis es menor que el de CYP2C9 y VKORC1.

## 3. Tecnología: NGS dirigido

**Recomendación: panel de NGS dirigido (captura o amplicones) sobre las regiones completas de los 7 genes**, con profundidad suficiente para detectar variantes de número de copias.

| Opción | A favor | En contra | Veredicto |
|---|---|---|---|
| **Array de genotipado** | Costo más bajo, alto rendimiento, interpretación sencilla | **Solo detecta las variantes diseñadas en el chip**, catalogadas mayoritariamente en poblaciones europeas. No detecta variantes raras ni nuevas. Difícil con repeticiones (UGT1A1\*28) y CNV | **No** para Colombia |
| **NGS dirigido** | Cubre el gen completo: detecta variantes raras, nuevas y combinaciones. Permite reinterpretar desde los datos guardados. Escalable | Costo intermedio, requiere infraestructura bioinformática | **Sí** |
| **WGS** | Cubre todo, reutilizable para tamizaje monogénico y PRS | Costo y almacenamiento altos; exige gobernanza de hallazgos incidentales; innecesario para un MVP | Fase posterior |

**La justificación viene de este mismo trabajo.** El análisis de DPYD (sección D3) mostró que **87 individuos portan variantes de función disminuida que el reporte automático no señala**, porque están en combinaciones de variantes sin alelo nombrado. En ancestría africana la subdetección fue de **34 veces** (1 de 34 portadores detectados), impulsada por `c.557A>G`, una variante frecuente en esa ancestría. Un array diseñado sobre catálogos europeos **no llevaría esa variante**, así que el problema sería aún mayor y además invisible. El NGS dirigido, al secuenciar el gen completo, sí la captura.

**Requisitos técnicos mínimos:**
- Profundidad ≥ 100× en las regiones de interés, incluidas las intrónicas que definen alelos (p. ej. HapB3 en `c.1129-5923C>G`).
- Detección de CNV y de la repetición TA de UGT1A1 cuando se incluya en fase 2.
- Interpretación con **PharmCAT**, con versión de guías CPIC registrada en cada informe.
- **Regla explícita de indeterminados:** todo resultado indeterminado se marca para revisión manual antes de emitirse, con la lista de variantes detectadas y su función según CPIC.

## 4. Evidencia que respalda el panel

1. **Guías CPIC de nivel A** para los 7 genes (verificado en la base de datos de CPIC). Nivel A significa que existe una recomendación clínica concreta y accionable.
2. **Ensayo PREPARE** (consorcio U-PGx, *The Lancet* 2023): panel preventivo de 50 variantes en 12 genes, 6.944 pacientes, diseño aleatorizado por conglomerados con cruce. Entre quienes tenían resultado accionable, las reacciones adversas clínicamente relevantes bajaron de **27,7 % a 21,0 % (OR 0,70; IC 95 % 0,54–0,91)**. Es la mejor evidencia disponible de que un panel preventivo reduce daño.
3. **Reguladores:** la EMA recomienda evaluar el déficit de DPD antes de fluoropirimidinas; la FDA incluye información farmacogenética en las etiquetas de clopidogrel, warfarina, azatioprina y capecitabina, entre otras.
4. **Guías DPWG** (consorcio holandés), que cubren pares adicionales y coinciden con CPIC en los genes de este panel.

**Matiz honesto sobre PREPARE, que debe declararse:** el **97,7 % de los participantes era de ancestría europea, mediterránea o de Medio Oriente**. La mejor evidencia disponible **no se generó en una población como la colombiana**. Esto no invalida el panel, porque el mecanismo farmacológico es el mismo, pero sí obliga a **medir el desempeño localmente** en lugar de asumir que el efecto se traslada igual. Es la razón principal para estructurar el programa como un piloto con medición de desenlaces.

## 5. Limitaciones de usar frecuencias de 1000 Genomas para la población colombiana

Las frecuencias de la sección D2 sirven para **dimensionar y priorizar**, no para parametrizar un servicio clínico. Siete limitaciones, de la más a la menos evidente:

1. **CLM no es Colombia.** Son 94 personas de **un área metropolitana** (Medellín). Colombia tiene estructura poblacional regional marcada: andina, caribe, pacífica (alta ancestría africana) y amazónica (alta ancestría indígena). La sección B mostró que la proporción de ancestría varía mucho **incluso entre individuos de CLM**.
2. **Precisión insuficiente para planear.** Con n = 94, los intervalos de confianza son amplios: NUDT15 disminuido es 8,5 % con IC 95 % de **4,4 % a 15,9 %**. Estimar el costo de un programa sobre ese rango es poco defendible.
3. **El muestreo no fue poblacional.** 1000 Genomas reclutó voluntarios sanos con criterios de ancestría declarada, no una muestra representativa de la población.
4. **No hay datos clínicos.** 1000 Genomas no tiene diagnósticos, prescripciones ni desenlaces. No permite estimar cuántos afiliados **realmente reciben** cada fármaco, que es el paso 2 de los criterios de inclusión. Ese dato solo lo tiene SURA.
5. **Los catálogos de alelos también son eurocéntricos.** Las variantes están mejor caracterizadas en europeos, así que las variantes de función alterada propias de ancestrías indígena o africana están subrepresentadas. **Este trabajo lo demostró empíricamente** con la subdetección de 34 veces en DPYD en ancestría africana.
6. **Los genotipos no son de grado clínico.** El panel 30x está faseado e imputado estadísticamente; por eso la sección A mostró 0 % de datos faltantes, lo cual no refleja lo que ocurriría con muestras clínicas reales.
7. **Genes no cubiertos.** CYP2D6 no se puede tipificar de forma confiable desde un VCF, así que 1000 Genomas no permite estimar sus frecuencias para el panel.

**Qué hacer en su lugar.** Usar 1000 Genomas para el diseño inicial y, en el piloto, **generar frecuencias alélicas propias de la población afiliada, estratificadas por región**. Esa base de datos local es, por sí sola, un activo: hoy no existe para Colombia con este nivel de detalle, y es la entrada necesaria para la fase de expansión.

## 6. Ruta de expansión

| Fase | Añadir | Justificación | Requisito nuevo |
|---|---|---|---|
| **2** | **CYP2D6** | Codeína, tramadol, tamoxifeno, ondansetrón, antidepresivos tricíclicos y paroxetina. Alto volumen | Método adicional para deleciones, duplicaciones e híbridos (lectura larga o caller especializado) |
| **2** | **UGT1A1** | Irinotecán, atazanavir. 14,9 % metabolizador pobre en CLM | Detección de la repetición TA (\*28) |
| **2** | **CYP3A5** | Tacrolimus. **31 % de expresores en CLM vs. 10 % en EUR** | Ninguno. Por bajo volumen poblacional, puede ofrecerse como prueba dirigida en programas de trasplante en lugar de tamizaje |
| **3** | **HLA-B** | Abacavir, alopurinol, carbamazepina | Tipificación HLA; los SNP indicadores usados en arrays se validaron en europeos y su desempeño varía por ancestría |
| **3** | **G6PD** | **Primaquina y tafenoquina** (malaria), dapsona, rasburicasa, nitrofurantoína | Gen ligado al X; conviene complementar con actividad enzimática. **Prioridad regional**: focalizar en zonas endémicas de malaria (Pacífico, Amazonía) más que a nivel nacional |
| **3** | CYP2B6, RYR1/CACNA1S, MT-RNR1 | Efavirenz; hipertermia maligna; ototoxicidad por aminoglucósidos | Ninguno relevante |

## 7. Implementación y métricas

**Sin estos tres componentes el panel no genera valor:**

1. **Integración con la prescripción.** Alerta automática en el momento de formular, no un PDF archivado. Es el punto donde fracasan la mayoría de los programas preventivos.
2. **Política de reinterpretación.** Las guías CPIC y los valores de actividad cambian con el tiempo: en este trabajo, el valor de CYP2D6\*41 pasó de 0,5 en la guía de 2018 a 0,25 en la tabla vigente. Debe definirse cada cuánto se reinterpretan los datos guardados y cómo se notifica un cambio de resultado.
3. **Consentimiento y privacidad.** Un dato genético guardado de por vida exige reglas de acceso, uso secundario y gobernanza, bajo la normativa colombiana de datos sensibles.

**Indicadores del piloto:**

| Indicador | Por qué |
|---|---|
| % de afiliados del piloto genotipados | Cobertura |
| % de prescripciones de fármacos del panel con resultado disponible | Mide si el modelo preventivo funciona |
| % de alertas accionables atendidas por el médico | Mide adopción real, no solo entrega |
| **Tasa de indeterminados por ancestría** | **Indicador de equidad**, derivado del hallazgo de DPYD |
| Eventos adversos evitados y hospitalizaciones | Desenlace y base del modelo económico |
| Costo por evento adverso evitado | Sustento para escalar |

**Modelo económico.** El ahorro proviene de hospitalizaciones y eventos adversos evitados, de fallas terapéuticas evitadas (clopidogrel sin efecto en metabolizadores pobres) y de menos pruebas repetidas, porque el genotipo se hace una sola vez. El costo por afiliado se amortiza entre todos los fármacos del panel a lo largo de la vida, que es la ventaja estructural del modelo preventivo frente al reactivo. Las cifras concretas requieren los datos de prescripción y de costo de SURA.
