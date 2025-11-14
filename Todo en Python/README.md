# Análisis Nodal IPR/VLP en Python

## 1. Objetivo del Proyecto

Este proyecto es una migración y refactorización completa de un script monolítico de análisis nodal de MATLAB a una arquitectura de Python moderna, modular y robusta.

El objetivo principal no era solo traducir el código, sino **validar, corregir y mejorar** la física subyacente, creando un motor de simulación estable y fácil de usar para análisis de sensibilidad.

### Características Principales

* **Arquitectura Modular:** El código está separado en "cerebros" (el motor) y "física" (las bibliotecas), facilitando futuras actualizaciones.
* **Validación Cruzada:** Las correlaciones `Hagedorn_Brown` y `Gray` han sido validadas 1:1 contra los *benchmarks* del script original de MATLAB.
* **Corrección de Física:** Se identificaron y repararon implementaciones numéricamente inestables y físicamente incorrectas de `Beggs_Brill` y `Ansari` que existían en el código original. El motor de Python es ahora **más preciso y estable**.
* **Control por Configuración:** Toda la simulación se controla desde un archivo `datos_pozo.json` externo. Ya no es necesario editar el código fuente de Python para cambiar los datos de un pozo o los modelos a ejecutar.

---

## 2. Arquitectura del Software

El proyecto está organizado en 9 archivos principales, cada uno con una responsabilidad única:

### Archivos de Ejecución y Configuración

* **`run_analisis.py` (Su Aplicación Principal)**
    * **Propósito:** Este es el **único script que necesitas ejecutar**.
    * **Función:** 1) Lee el `json`, 2) Llama al `core.py` para correr la simulación, y 3) Usa `matplotlib` para graficar los resultados.

* **`datos_pozo.json` (Su Panel de Control)**
    * **Propósito:** Este archivo reemplaza el "Bloque 1" de MATLAB. Aquí es donde se **ingresan todos los datos** del pozo y se controla qué simulaciones correr.

* **`test_engine.py` (Banco de Pruebas)**
    * **Propósito:** Este es el script que usamos para validar la migración. Compara los resultados de Python (PVT, Puntos de Operación) contra los *benchmarks* conocidos de MATLAB para asegurar que la física es idéntica.

### Archivos del "Motor" (El Cerebro)

* **`core.py` (El Orquestador)**
    * **Propósito:** Contiene la lógica principal (`correr_simulacion_completa`). Es el "cerebro" que sabe *cuándo* llamar al IPR, *cuándo* llamar al VLP y *cómo* encontrar la intersección.

* **`vlp_loop.py` (El Bucle VLP)**
    * **Propósito:** Contiene la función `calcular_vlp_para_q`. Esta es la traducción directa del bucle de 50 segmentos de MATLAB. Itera desde el cabezal (`Pcabezal`) hasta el fondo, llamando al `pvt.py` y `flow.py` en cada segmento para calcular el gradiente de presión.

### Archivos de "Física" (Las Bibliotecas)

* **`pvt.py` (La Biblioteca PVT)**
    * **Propósito:** Contiene **toda la física de fluidos**.
    * **Funciones Clave:** `prop_oil_eos` (Bo, Rs, etc.), `dens_xfzg_tot` (densidades), `funcion_zg_exp` (Z-Factor). Está validado 1:1 con MATLAB.

* **`ipr.py` (La Biblioteca IPR)**
    * **Propósito:** Contiene todas las ecuaciones de **Oferta** (IPR).
    * **Funciones Clave:** `calcular_ipr_vogel`, `calcular_ipr_darcy_gas`, `calcular_ipr_lineal`, etc.

* **`flow.py` (La Biblioteca VLP)**
    * **Propósito:** Contiene todas las ecuaciones de **Demanda** (VLP) y correlaciones de flujo.
    * **Funciones Clave:** `calcular_vlp_hagedorn_brown`, `calcular_vlp_beggs_brill` (¡Corregida!), `calcular_vlp_gray`, y `calcular_vlp_ansari` (¡Reconstruida!).

* **`utils.py` (Utilidades)**
    * **Propósito:** Almacena el diccionario de `const` (constantes de conversión) para que sea accesible por todos los demás archivos.

---

## 3. Instalación (Setup)

Para correr este proyecto en una nueva computadora, solo necesitas 3 cosas:

1.  **Instalar Python:** Descarga Python 3.10 o más reciente desde `python.org`.
    * **¡Importante!** Durante la instalación, asegúrate de marcar la casilla que dice **"Add Python to PATH"**.

2.  **Instalar Bibliotecas:** Abre una terminal (`cmd` o `PowerShell`) y ejecuta los siguientes 3 comandos uno por uno:
    ```bash
    py -m pip install numpy
    py -m pip install matplotlib
    py -m pip install scipy
    ```

3.  **Copiar los Archivos:** Coloca todos los 9 archivos de Python (`.py`) y el `.json` en la misma carpeta.

---

## 4. ¿Cómo Usar la Herramienta? (Guía Rápida)

Este es tu nuevo flujo de trabajo:

**Paso 1: Configurar la Simulación**
Abre el archivo `datos_pozo.json` en un editor de texto o VS Code.

**Para cambiar los datos del pozo:**
Edita cualquier valor en las secciones `general`, `fluidos`, `yacimiento` o `calibracion`.

**Para seleccionar el IPR:**
En la sección `modelo`, cambia el `IPR_Tipo`:
* `1`: Lineal
* `2`: Vogel (usa `Ql_cal` y `Pwf_cal`)
* `3`: Backpressure (usa `Ql_cal` y `Pwf_cal`)
* `4`: Darcy Gas p² (ignora calibración, usa `k`, `h`, etc.)

**Para seleccionar las VLP (Tu Panel de Control):**
Ve a la sección `vlp_sensibilidad` y usa `1` (encendido) o `0` (apagado) para cada correlación que quieras comparar.

```json
    "modelo": {
        "IPR_Tipo": 2,
        ...
        "vlp_sensibilidad": {
            "Hagedorn_Brown": 1,
            "Beggs_Brill": 1,
            "Ansari": 1,
            "Gray": 0
        }
    }
    