# G. Caso de negocio: Programa de Genómica Poblacional para SURA

## G1. Propuesta de proyecto

### El problema que resuelve

SURA gasta hoy en tres cosas que la genómica puede reducir:

1. **Eventos adversos a medicamentos.** Hospitalizaciones, urgencias y días de incapacidad por reacciones prevenibles.
2. **Fallas terapéuticas silenciosas.** Pacientes que reciben un fármaco que en ellos no funciona, por ejemplo clopidogrel en metabolizadores pobres de CYP2C19, y que vuelven al sistema con el evento que se quería prevenir.
3. **Diagnóstico tardío de condiciones prevenibles.** Cáncer de mama detectado en estadio avanzado en una portadora de BRCA que nunca supo que lo era; un infarto en un paciente con hipercolesterolemia familiar no diagnosticada.

**La ventaja estructural de SURA.** Un laboratorio independiente que vende una prueba genética paga el costo y **otro actor captura el ahorro**. SURA integra aseguradora, EPS, IPS y laboratorio: **paga la prueba y captura el ahorro que esa prueba genera**. Es la razón por la que este programa tiene sentido económico en SURA y no en la mayoría de organizaciones. Ese argumento debe encabezar la presentación a la Gerencia.

### Objetivos

| | Objetivo | Indicador |
|---|---|---|
| Primario | Reducir eventos adversos a medicamentos prevenibles en la población afiliada | Tasa de eventos adversos en fármacos del panel |
| Secundario | Detectar portadores de condiciones accionables antes del evento clínico | Portadores identificados y adherencia al tamizaje indicado |
| Habilitador | Generar la **primera base de frecuencias farmacogenéticas de población colombiana** con datos propios | Cobertura por región y número de variantes caracterizadas |

El tercer objetivo suele omitirse y es el más valioso a largo plazo: este trabajo mostró que **las frecuencias de referencia disponibles no representan a Colombia**, y esa base de datos no existe hoy.

### A quién va dirigido

- **Fase piloto:** afiliados que ya están en contacto con el sistema por una de las condiciones de alto impacto, es decir pacientes cardiovasculares, oncológicos y de enfermedad inflamatoria intestinal. Tienen la mayor probabilidad de recibir un fármaco del panel en el corto plazo, lo que permite medir impacto rápido.
- **Fase de expansión:** afiliados adultos en general, priorizando por edad y por carga de prescripción.

### Métodos, herramientas y alcance

Todo el componente técnico está descrito y ejecutado en las secciones A a D de este informe:

| Componente | Método | Estado en este ejercicio |
|---|---|---|
| Genotipado | NGS dirigido, panel de 7 genes (ver `D_panel_pgx_sura.md`) | Diseñado |
| Interpretación PGx | PharmCAT, guías CPIC con versión registrada | Ejecutado sobre 2.504 individuos |
| Variantes clínicas | ClinVar con filtros explícitos de clasificación, nivel de revisión y frecuencia | Ejecutado sobre 6 genes ACMG SF |
| Ancestría | PCA y ADMIXTURE | Ejecutado |
| Reproducibilidad | conda con versiones fijadas, scripts numerados, semillas fijas | Implementado |

**Fuera de alcance del programa** (declarado de entrada, no descubierto sobre la marcha): secuenciación de genoma completo, puntajes poligénicos con uso clínico, y nutrigenómica. Las razones están en G2.

### Modelo de ingreso y de ahorro

**Ahorro** (aseguradora y EPS), que es la fuente principal:
- Hospitalizaciones y urgencias evitadas por eventos adversos.
- Fallas terapéuticas evitadas: menos reintervenciones y menos progresión de enfermedad.
- Detección temprana: el tratamiento de un cáncer en estadio temprano cuesta una fracción del de uno avanzado.
- Prueba única de por vida: el genotipo no se repite, así que el costo se amortiza entre todos los fármacos futuros del afiliado.

**Ingreso** (IPS y laboratorio):
- Servicio de genotipado e interpretación a terceros: otras EPS, aseguradoras y medicina prepagada.
- Diferenciador en planes de medicina prepagada.
- Investigación con socios académicos e industria, con consentimiento explícito y datos agregados.

**Lo que hay que calcular antes de comprometer cifras** (requiere datos internos de SURA): volumen anual de prescripción de los fármacos del panel, costo actual de los eventos adversos asociados, y costo por afiliado del genotipado con proveedor local. Con esos tres números se estima el costo por evento adverso evitado, que es la métrica que decide si escalar.

## G2. MVP del programa de genómica preventiva

**Principio de selección:** entra al MVP lo que tiene **evidencia de que cambia desenlaces** y **una acción clínica definida**. Lo demás espera, no porque sea falso, sino porque todavía no se sabe qué hacer con el resultado.

| Componente | Evidencia | ¿En el MVP? | Razón |
|---|---|---|---|
| **Farmacogenómica preventiva** (7 genes) | Ensayo PREPARE: reacciones adversas de 27,7 % a 21,0 %, OR 0,70 (IC 95 % 0,54–0,91), 6.944 pacientes. Guías CPIC nivel A | **Sí** | Es el único componente con un ensayo aleatorizado que demuestra reducción de daño. Acción clínica inmediata y protocolizada |
| **Condiciones monogénicas accionables** (HBOC, Lynch, hipercolesterolemia familiar) | Las tres son **CDC Tier 1**: evidencia suficiente para aplicación poblacional | **Sí, pero con alcance acotado** | Ver abajo |
| **Puntajes poligénicos (PRS)** | Sin ensayos que demuestren mejora de desenlaces. Precisión cae ~2 veces en latinoamericanos y ~4,5 en africanos (Martin et al., 2019) | **No** | Un resultado que no se sabe traducir en conducta clínica genera ansiedad y consultas sin beneficio. Línea de investigación |
| **Nutrigenómica** | Food4Me (n = 1.607, 4 brazos): *"no hubo evidencia de que incluir información fenotípica ni fenotípica más genotípica aumentara la efectividad del consejo nutricional personalizado"* | **No** | La evidencia disponible indica que el componente genético no aporta sobre la personalización sin genotipo |

**Alcance acotado de las condiciones monogénicas.** En el MVP se implementan las dos estrategias de mayor rendimiento y menor costo, no el tamizaje poblacional indiscriminado:

1. **Prueba en cascada:** cuando se identifica un portador, se ofrece la prueba a sus familiares de primer grado. Es la estrategia de mayor rendimiento por prueba realizada y ya está aceptada clínicamente.
2. **Prueba por criterios:** pacientes que cumplen criterios de historia familiar o personal.

El **tamizaje poblacional no seleccionado** queda para la fase 2, y **este trabajo aporta un argumento propio para esa decisión**: en la sección C se detectaron 2 portadores de variantes patogénicas de BRCA1/BRCA2 en 2.504 personas, frente a los 5 a 8 esperados. La razón es que solo se detectan variantes **ya catalogadas**, y la mayoría de las variantes patogénicas de BRCA son privadas de una familia. Un tamizaje poblacional con catálogos actuales tendría un rendimiento menor al que promete la literatura, y **desigual entre ancestrías**. Antes de escalarlo hay que construir capacidad de clasificación de novo con criterios ACMG/AMP.

**Lo que el MVP debe incluir aunque no sea "genómica":**

- **Integración con la prescripción.** Si la alerta no aparece cuando el médico formula, el programa no genera valor. Es el punto donde más programas preventivos fracasan.
- **Política de indeterminados.** La sección D3 mostró que **87 personas con variantes DPYD de función disminuida quedaron sin alerta** en el reporte automático, con una subdetección de 34 veces en ancestría africana. Un reporte automático no se entrega sin auditoría.
- **Política de reinterpretación.** Las clasificaciones cambian: en este trabajo, el valor de actividad de CYP2D6\*41 pasó de 0,5 en la guía de 2018 a 0,25 en la tabla vigente.
- **Métricas de equidad.** Tasa de indeterminados y de resultados accionables **por ancestría**, como indicador de calidad del servicio.

## G3. Respuesta a la propuesta comercial de un test nutrigenético directo al consumidor

> *El área comercial propone vender directamente al consumidor un test nutrigenético que promete "la dieta ideal según tu ADN".*

### Recomendación: no, en esos términos. Pero la intuición comercial es correcta y se puede capturar.

**Por qué no, en tres frentes:**

**1. La promesa no está respaldada.** El ensayo Food4Me asignó al azar a 1.607 personas en cuatro brazos: consejo estándar, consejo personalizado por dieta, por dieta más fenotipo, y por dieta más fenotipo más genotipo. Conclusión textual: *"no hubo evidencia de que incluir información fenotípica, ni fenotípica más genotípica, aumentara la efectividad del consejo nutricional personalizado"*. Es decir, **la personalización funcionó, pero el genotipo no agregó nada**.

**2. Riesgo regulatorio y de responsabilidad.** Vender una promesa de resultado en salud sin evidencia expone a SURA ante la autoridad sanitaria y ante reclamaciones de consumidores. Y como el test se vendería directo al consumidor, sin médico de por medio, SURA queda como responsable único de la interpretación.

**3. Riesgo reputacional, que es el mayor.** SURA no es una empresa de bienestar: es una organización de salud con credibilidad clínica. Vender un producto que la comunidad científica considera sin sustento **daña la credibilidad del resto del programa genómico**, incluido el componente farmacogenético que sí funciona. Se pondría en riesgo lo que tiene evidencia para vender lo que no la tiene.

### La alternativa que sí captura la oportunidad

El área comercial detectó algo real: **existe demanda de personalización en salud y la gente está dispuesta a pagar por ella**. La respuesta no es negar la demanda, sino atenderla con lo que sí funciona.

**Producto propuesto: un programa de prevención personalizada**, con tres componentes, dos de ellos genéticos y uno no:

| Componente | Qué promete | Respaldo |
|---|---|---|
| **Farmacogenómica preventiva** | "Cómo responde tu cuerpo a los medicamentos que podrías necesitar" | Ensayo PREPARE, guías CPIC |
| **Riesgo hereditario accionable** | "Condiciones familiares que se pueden prevenir si se detectan a tiempo" | CDC Tier 1 |
| **Nutrición personalizada sin genotipo** | "Un plan alimentario construido sobre tus hábitos, tus mediciones y tu perfil metabólico" | Food4Me: **la personalización sí funcionó**; lo que no aportó fue el genotipo |

El tercer componente es la clave de la respuesta comercial: **se conserva la personalización, que es lo que el cliente quiere y lo que tiene evidencia, y se retira el genotipo, que es lo que no la tiene**. El producto sigue siendo vendible y además es honesto.

**Cómo comunicarlo a la Gerencia, en una frase:**

> *"La demanda que identificó el área comercial es real, y la podemos atender. Lo que no podemos sostener es la promesa específica de que el ADN define la dieta óptima: el ensayo de referencia mostró que el genotipo no aporta sobre la personalización sin genotipo. Propongo vender la personalización, que sí funciona, junto con los dos componentes genéticos que sí tienen evidencia. Ganamos el mismo mercado sin arriesgar la credibilidad del programa."*

## Referencias

- Swen JJ et al. A 12-gene pharmacogenetic panel to prevent adverse drug reactions: an open-label, multicentre, controlled, cluster-randomised crossover implementation study. *The Lancet* 2023;401:347–356.
- Celis-Morales C et al. Effect of personalized nutrition on health-related behaviour change: evidence from the Food4Me European randomized controlled trial. *Int J Epidemiol* 2017;46(2):578–588.
- Martin AR et al. Clinical use of current polygenic risk scores may exacerbate health disparities. *Nature Genetics* 2019;51:584–591.
- CPIC, guías vigentes: cpicpgx.org. CDC, Tier 1 Genomics Applications.
