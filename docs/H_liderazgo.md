# H. Liderazgo técnico

## H1. Acompañar el crecimiento de un bioinformático junior

**El principio:** un junior no necesita que le resuelvan los problemas, sino que le enseñen **a decidir**. En bioinformática la parte difícil casi nunca es correr la herramienta: es saber qué filtro aplicar, con qué umbral y por qué.

- **Enseñar el porqué, no el comando.** Cada decisión técnica se documenta con su justificación. En este mismo trabajo, cada script lleva al inicio las razones de sus filtros: por qué no se aplica Hardy-Weinberg sobre poblaciones mezcladas, por qué se excluye la inversión 17q21.31, por qué `--absent-to-ref`. Esas líneas son material de formación, no burocracia.
- **Revisión de código como conversación, no como control.** La pregunta que más enseña es "¿qué pasaría si este umbral fuera otro?". Revisar preguntando, no corrigiendo.
- **Darle una pieza completa y propia**, de extremo a extremo, aunque sea pequeña: desde los datos hasta la figura y el párrafo de interpretación. Aprender a cerrar algo enseña más que ejecutar diez pasos sueltos.
- **Entrenar la desconfianza productiva.** La habilidad que más separa a un bioinformático senior es sospechar de su propio resultado. En este trabajo, el 26 % de indeterminados en DPYD parecía un detalle técnico y resultó ser una subdetección de casos accionables. Enseñar a preguntarse "¿este número tiene sentido?" vale más que enseñar diez herramientas.
- **Autonomía gradual y explícita.** Decir qué decisiones ya puede tomar solo y cuáles todavía se consultan, y revisar ese límite periódicamente. La ambigüedad sobre la autonomía es lo que más frustra a un junior.
- **Normalizar el error.** Contar los propios. En este ejercicio cambié dos veces el valor de actividad de CYP2D6\*41 antes de verificarlo en la fuente vigente. Un equipo donde el líder no reconoce errores produce gente que los esconde.

## H2. Aportar a un bioinformático experto

Con un experto el rol se invierte: no se trata de enseñar, sino de **multiplicar su capacidad**.

- **Quitar fricción.** Infraestructura, acceso a datos, permisos, reuniones innecesarias. Cada hora que un experto pierde en trámites es la hora más cara del equipo.
- **Darle los problemas difíciles, no los urgentes.** Su valor está en lo que nadie más puede resolver, no en lo que más grita.
- **Ser un par que revisa, no un jefe que aprueba.** Aportar la pregunta incómoda: "¿esto se sostiene en población colombiana?", "¿qué pasa si el catálogo está sesgado?". Esa pregunta es el aporte, aunque la respuesta técnica sea suya.
- **Conectar su trabajo con la decisión clínica o de negocio.** Traducir el impacto hacia arriba es algo que el líder técnico puede hacer y que suele faltarle al experto.
- **Aprender de él explícitamente.** Pedirle que enseñe al equipo es a la vez desarrollo para los demás y reconocimiento para él.

## H3. Tres solicitudes urgentes al mismo tiempo: clínica, investigación y comercial

**Criterio de priorización, en este orden:**

1. **Irreversibilidad y riesgo para un paciente.** ¿Hay alguien esperando una decisión que afecte su tratamiento? Si la respuesta es sí, eso va primero, siempre.
2. **Urgencia real frente a urgencia declarada.** Preguntar por la fecha límite verdadera y qué ocurre si se entrega un día después. Muchas urgencias se disuelven con esa pregunta.
3. **Costo de la demora.** Una fecha de sometimiento de un artículo o un comité comercial son fechas duras, pero su incumplimiento se recupera; una decisión clínica tardía, no.

**Aplicado al caso:**

| Solicitud | Decisión | Razón |
|---|---|---|
| **Clínica** | **Primera** | Hay un paciente esperando; la decisión es sensible al tiempo y potencialmente irreversible |
| **Investigación** | Segunda, con fecha acordada | Suele tener una fecha límite real, pero de días o semanas, no horas |
| **Comercial** | Tercera, o se entrega una versión parcial | Casi siempre admite un entregable preliminar que desbloquea la conversación |

**Cómo se comunica**, que importa tanto como la priorización:

- **Ese mismo día, a los tres, por escrito.** El daño no lo causa la espera sino el silencio.
- **Decir el orden, el motivo y la hora estimada de entrega.** "Voy a atender primero el caso clínico porque hay un paciente esperando. Tu resultado lo tienes mañana a las 11." Un compromiso concreto tranquiliza; un "estoy en eso" no.
- **Ofrecer un entregable parcial** cuando exista: unos resultados preliminares suelen bastar para que el área comercial avance.
- **Escalar si los tres son realmente inaplazables.** Si no hay forma de cumplir, la decisión sobre qué se sacrifica no es del bioinformático: es de quien asume la consecuencia. Escalar a tiempo no es no saber priorizar, es evitar que la organización se entere tarde.
- **Registrar el patrón.** Si esta colisión se repite, el problema no es de priorización sino de capacidad, y eso es una conversación distinta con la Gerencia, con datos sobre la frecuencia de los conflictos.
