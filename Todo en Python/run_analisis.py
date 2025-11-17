# =====================================================================
# SCRIPT DE ANÁLISIS DE SENSIBILIDAD (VLP-Sens)
# =====================================================================

import json
import numpy as np
import matplotlib.pyplot as plt
from core import correr_simulacion_completa # Importamos nuestro motor validado
import copy 

def cargar_datos(json_file):
    """ Carga los datos desde un archivo JSON y los 'aplana' en un solo dict """
    with open(json_file, 'r') as f:
        datos_anidados = json.load(f)
    
    datos = {}
    datos.update(datos_anidados['general'])
    datos.update(datos_anidados['fluidos'])
    datos.update(datos_anidados['yacimiento'])
    datos.update(datos_anidados['calibracion'])
    datos.update(datos_anidados['modelo'])
    datos['API_o_SG'] = datos['API']
    
    return datos

def plot_sensibilidad_vlp(resultados):
    """ 
    Genera la gráfica IPR vs MÚLTIPLES VLPs con lógica de unidades de gas 
    basada únicamente en la bandera 'Yacimiento_Gas'.
    """
    
    # 1. Extraer datos del IPR Base
    curva_ipr = resultados['curva_ipr']
    ipr_usado_tipo = resultados['ipr_usado']
    vlp_resultados_lista = resultados['lista_vlps']
    
    # --- LOGICA DE UNIDADES CORREGIDA ---
    es_yacimiento_gas = resultados.get('yacimiento_gas', 0)
    RGL_scfstb = resultados['RGL']
    
    if es_yacimiento_gas == 1: # Si la bandera de Gas está en 1, SIEMPRE usar unidades de gas
        # Conversion: bpd * RGL (scf/stb) = scf/d
        # Escala: scf/d / 1e9 = MMMscf/d (Mil Millones de scf/día)
        escala_a_plotear = RGL_scfstb / 1e9 
        x_label_str = 'Gasto de Gas (MMMscf/día)'
        leyenda_unidad = 'MMMscf/día'
    else: # Si es Líquido/Aceite
        escala_a_plotear = 1.0
        x_label_str = 'Gasto de Líquido (bpd)'
        leyenda_unidad = 'bpd'
        
    # Definición de etiqueta IPR (Solo estética)
    if ipr_usado_tipo == 2: ipr_label = 'IPR (Vogel)'
    elif ipr_usado_tipo == 3: ipr_label = 'IPR (Backpressure)'
    elif ipr_usado_tipo == 4: ipr_label = 'IPR (Darcy Gas p^2)'
    else: ipr_label = f'IPR (Tipo {ipr_usado_tipo})'
    # --- FIN LOGICA DE UNIDADES ---
        
    q_ipr = np.array([p[0] for p in curva_ipr])
    p_ipr = np.array([p[1] for p in curva_ipr])
    
    # 2. Definir colores y estilos para las curvas VLP
    colores = ['r', 'g', 'm', 'c', 'orange']
    estilos = ['-', '--', ':', '-.', '-']
    marcadores = ['o', 's', '^', 'd', 'v']
    
    # 3. Crear la gráfica
    plt.figure(figsize=(10, 7))
    
    # Graficar el IPR Base (una sola vez)
    plt.plot(q_ipr * escala_a_plotear, p_ipr, 'b-', linewidth=2.5, label=ipr_label)
    
    # 4. Iterar y graficar cada resultado de VLP
    for i, resultado in enumerate(vlp_resultados_lista):
        
        color = colores[i % len(colores)]
        estilo = estilos[i % len(estilos)]
        marcador = marcadores[i % len(marcadores)]
        
        curva_vlp = resultado['curva_vlp']
        (q_op, p_op) = resultado['punto_op']
        etiqueta_vlp = resultado['nombre']
        
        if not curva_vlp: 
            print(f"Advertencia: No hay datos de VLP para '{etiqueta_vlp}'.")
            continue
            
        q_vlp = np.array([p[0] for p in curva_vlp])
        p_vlp = np.array([p[1] for p in curva_vlp])

        # Construir la etiqueta de leyenda
        if np.isfinite(q_op):
            label_completa = f'VLP ({etiqueta_vlp}) (Op: {q_op * escala_a_plotear:.2f} {leyenda_unidad})'
        else:
            label_completa = f'VLP ({etiqueta_vlp}) (Sin cruce)'

        # Graficar Curva VLP con escala aplicada
        plt.plot(q_vlp * escala_a_plotear, p_vlp, color=color, linestyle=estilo, linewidth=2.0, label=label_completa)
        
        # Graficar Punto de Operación con escala aplicada
        if np.isfinite(q_op):
            plt.plot(q_op * escala_a_plotear, p_op, marker=marcador, markersize=10, 
                     markerfacecolor=color, markeredgecolor='k', 
                     linestyle='none') 
    
    # 5. Formato Final
    plt.title('Análisis de Sensibilidad IPR vs VLP - Python')
    plt.xlabel(x_label_str) # <-- Usar la nueva etiqueta de Gas/Líquido
    plt.ylabel('Pwf (psi)')
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    
    # 6. Mostrar la gráfica
    print("Mostrando gráfica de sensibilidad...")
    plt.show()

# --- Función Principal ---
def main():
    archivo_base = 'datos_pozo.json'
    
    print(f"Cargando Caso Base desde: {archivo_base}")
    datos_base = cargar_datos(archivo_base)
    
    print("Corriendo simulación de sensibilidad VLP...")
    
    # 1. Correr la simulación completa
    resultados_finales = correr_simulacion_completa(datos_base)
    
    print("Simulaciones completadas. Generando gráfica...")
    # 2. Llamar a la función de ploteo
    plot_sensibilidad_vlp(resultados_finales)
    
    print("--- Proceso Finalizado ---")

if __name__ == "__main__":
    main()