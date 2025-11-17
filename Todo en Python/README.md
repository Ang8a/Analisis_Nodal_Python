# ⛽ Análisis Nodal IPR/VLP en Python (FINAL)

## 1. Objetivo del Proyecto

Este proyecto es el resultado de la migración, refactorización y **corrección profunda** del script monolítico de MATLAB para el análisis IPR/VLP.

El resultado es un motor de simulación en Python que es **más estable, preciso y físicamente correcto** que el código fuente original.

### Correcciones y Mejoras Clave

| Correlación / Componente | Estado Final | Valor Añadido |
| :--- | :--- | :--- |
| **Arquitectura** | Finalizada | Migración a estructura modular (clean code) con control `.json`. |
| **Hagedorn & Brown** | ✅ Validada 1:1 | Coincidencia numérica perfecta con el benchmark de MATLAB. |
| **Beggs & Brill** | 🛠️ **Corregida** | Se identificó y parchó el bug de inestabilidad que causaba la explosión numérica. La curva VLP ahora es físicamente estable. |
| **Ansari (Mecanicista)** | 🛠️ **Reconstruida** | Se reemplazó la lógica rota original por una implementación *híbrida-mecanicista* que usa los criterios de patrón de flujo de Ansari para rutear el cálculo al modelo empírico más estable (Gray/B&B). |
| **Unidades de Gas** | ✅ Implementada | La gráfica ahora maneja la conversión automática de unidades a **MMMscf/día** cuando se activa el switch de Gas. |

---

## 2. Arquitectura del Software

El código está organizado en módulos desacoplados para fácil mantenimiento:

* **`pvt.py`**: Biblioteca de Física de Fluidos (Validado: Z-Factor, $P_b$, $R_s$, etc.).
* **`flow.py`**: Biblioteca de Correlaciones de Flujo (Contiene la lógica *corregida* de B&B y el nuevo modelo híbrido de Ansari).
* **`core.py`**: El Orquestador. Contiene la lógica de la simulación principal y los *switches* de correlación.
* **`vlp_loop.py`**: Motor de VLP. Contiene el bucle de 50 segmentos que calcula el perfil de presión.
* **`run_analisis.py`**: Script de ejecución. Lee el `.json` y genera el plot final de `matplotlib`.

---

## 3. Guía de Uso y Configuración

El control de la simulación se realiza editando el archivo **`datos_pozo.json`**.

### A. Control de Correlaciones (Su Panel de Control)

Utilice la sección `"vlp_sensibilidad"` para elegir qué VLP desea comparar en la gráfica.

| Flag | Correlación | Estado / Notas |
| :--- | :--- | :--- |
| **`Hagedorn_Brown`** | Empírica | ✅ Funcional y estable (Recomendada para Aceite). |
| **`Beggs_Brill`** | Empírica | ✅ Funcional y estable (Corregida). |
| **`Ansari`** | Mecanicista | 🛠️ Reconstruida. Utiliza la lógica de patrones para el cálculo. |
| **`Gray`** | Empírica | ✅ Funcional y estable (Recomendada para Gas/Condensado). |

### B. Control del Modo de Unidades (Gas vs. Líquido)

Para que el Eje X cambie de `bpd` a `MMMscf/día`, controle la bandera **`"Yacimiento_Gas"`** en la sección `modelo`.

| Parámetro | Valor | Resultado en la Gráfica |
| :--- | :--- | :--- |
| **`"Yacimiento_Gas"`** | **0** | El Eje X muestra **`Gasto de Líquido (bpd)`**. |
| **`"Yacimiento_Gas"`** | **1** | El Eje X muestra **`Gasto de Gas (MMMscf/día)`**. |

### C. Ejecución Final

1.  Abre y edita **`datos_pozo.json`**.
2.  Guarda el archivo.
3.  Ejecuta en la terminal: `py run_analisis.py`

---

## 4. Conclusión

El proyecto de migración está formalmente concluido.