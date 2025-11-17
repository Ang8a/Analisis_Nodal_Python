# =====================================================================
# ESTE ES EL CÓDIGO COMPLETO PARA 'core.py' (CON ANSARI ACTIVADO)
# Bórralo todo y pega esto.
# =====================================================================
import numpy as np
from pvt import prop_oil_eos
from ipr import (calcular_ipr_lineal, calcular_ipr_vogel, calcular_ipr_backpressure, calcular_ipr_darcy_gas)
from vlp_loop import calcular_vlp_para_q
from utils import const
import copy

def _set_vlp_flags(datos, nombre_corr):
    """ Función helper para activar la correlación VLP correcta """
    datos['h_rb_MFF_HB'] = 0
    datos['h_radiobutBeggsBrill'] = 0
    datos['h_radiobutton_Ansari'] = 0
    datos['h_CFM_GRAY'] = 0
    
    # Encender la seleccionada
    if nombre_corr == 'Hagedorn_Brown':
        datos['h_rb_MFF_HB'] = 1
    elif nombre_corr == 'Beggs_Brill':
        datos['h_radiobutBeggsBrill'] = 1
    elif nombre_corr == 'Ansari':
        datos['h_radiobutton_Ansari'] = 1
    elif nombre_corr == 'Gray':
        datos['h_CFM_GRAY'] = 1
    
    return datos

def calcular_pb_yacimiento(datos):
    """ Calcula la Pb usando Standing Inverso """
    temp_F_yac = datos['Tfondo']
    api_yac = datos['API']
    gg_yac = datos['Sg']
    glr_yac = datos['RGL']
    
    X_standing = (0.0125 * api_yac - 0.00091 * ((temp_F_yac - 32) * 5 / 9))
    F_standing = 10**(X_standing * 1.2048)
    
    Pb_yacimiento_calculada = ((glr_yac / (gg_yac * F_standing)) - 1.4) * 18.2
    return max(Pb_yacimiento_calculada, 14.7)

def encontrar_interseccion(curva_IPR, curva_VLP):
    """ Traducción de la función 'encontrar_interseccion' de MATLAB """
    if not curva_IPR or not curva_VLP or len(curva_IPR) < 2 or len(curva_VLP) < 2:
        return (np.nan, np.nan)
        
    Q_ipr = np.array([p[0] for p in curva_IPR])
    P_ipr = np.array([p[1] for p in curva_IPR])
    Q_vlp = np.array([p[0] for p in curva_VLP])
    P_vlp = np.array([p[1] for p in curva_VLP])

    try:
        # Interpolar VLP a los puntos Q del IPR
        P_vlp_interp = np.interp(Q_ipr, Q_vlp, P_vlp)
    except Exception:
        return (np.nan, np.nan)

    diferencia = P_ipr - P_vlp_interp
    signo = np.sign(diferencia)
    idx_change = np.where(np.diff(signo) != 0)[0]

    if idx_change.size == 0:
        idxmin = np.argmin(np.abs(diferencia))
        min_diff = np.abs(diferencia[idxmin])
        if min_diff > 50: 
            return (np.nan, np.nan)
        else:
            return (Q_ipr[idxmin], P_ipr[idxmin])
    else:
        i1 = idx_change[-1] 
        i2 = i1 + 1
        
        Q1, Q2 = Q_ipr[i1], Q_ipr[i2]
        F1, F2 = diferencia[i1], diferencia[i2]
        
        if abs(F2 - F1) < 1e-9:
            return (Q_ipr[np.argmin(np.abs(diferencia))], P_ipr[np.argmin(np.abs(diferencia))])
        else:
            qop = Q1 + (Q2 - Q1) * (0 - F1) / (F2 - F1)
            
            # --- CORRECCIÓN FINAL DE COORDENADAS ---
            # En lugar de interpolar en el IPR, donde la P está invertida, 
            # interpolamos en la curva VLP (que es localmente estable y creciente).
            # Luego usamos el IPR para la validación interna (si es necesario).
            pwf_op_vlp = np.interp(qop, Q_vlp, P_vlp)
            
            return (qop, pwf_op_vlp) # Usamos la coordenada Pwf de la curva VLP

# ... (Código existente arriba) ...

def correr_simulacion_completa(datos):
    """ 
    Versión Python del 'motor' de simulación.
    Refactorizado para leer banderas 1/0 desde el JSON.
    """
    
    # 0. Unpack de datos
    pws = datos['Pws']
    pwf_cal = datos['Pwf_cal']
    ql_cal = datos['Ql_cal']
    twf = datos['Tfondo']
    api = datos['API']
    gg = datos['Sg']
    glr = datos['RGL']
    wc = datos['Porcentaje_Agua']
    YCO2, YH2s, YN2 = datos['YCO2'], datos['YH2s'], datos['YN2']
    sg = 141.5 / (131.5 + api)
    
    # 1. Calcular Pb (Solo una vez)
    Pb_yacimiento = calcular_pb_yacimiento(datos)
    
    # 2. Calcular IPR (Solo una vez)
    IPR_tipo = datos['IPR_Tipo']
    if datos['h_IPR_Automatico'] == 1:
        IPR_tipo = 1 if pws > Pb_yacimiento else 2
        
    if pwf_cal <= 0 or ql_cal <= 0:
        pwf_cal_ipr = pws * 0.9
        ql_cal_ipr = ql_cal 
    else:
        pwf_cal_ipr = pwf_cal
        ql_cal_ipr = ql_cal

    pwf_vec = np.linspace(0, pws, 101)
    q_vec = np.zeros(101)

    if IPR_tipo == 1:
        q_vec = calcular_ipr_lineal(datos, pws, pwf_cal_ipr, ql_cal_ipr, pwf_vec)
    elif IPR_tipo == 2:
        q_vec = calcular_ipr_vogel(datos, pws, pwf_cal_ipr, ql_cal_ipr, pwf_vec)
    elif IPR_tipo == 3:
        q_vec = calcular_ipr_backpressure(datos, pws, pwf_cal_ipr, ql_cal_ipr, pwf_vec)
    elif IPR_tipo == 4:
        q_vec = calcular_ipr_darcy_gas(datos, pws, twf, api, gg, sg, glr, wc, YCO2, YH2s, YN2, Pb_yacimiento, pwf_vec)
    
    q_vec = np.maximum(0, q_vec)
    curva_IPR = list(zip(q_vec, pwf_vec))
    
    # 3. Lógica de VLP (Bucle de Sensibilidad)
    qmax_ipr = np.max(q_vec)
    if qmax_ipr <= 0: qmax_ipr = 1000.0
    q_eje_vlp = np.linspace(qmax_ipr * 0.01, qmax_ipr * 1.1, 40)
    
    lista_resultados_vlp = []
    
    # Obtener la bandera de Gas (Nueva Lógica)
    es_yacimiento_gas = datos.get('Yacimiento_Gas', 0)
    
    correlaciones_a_correr = []
    vlp_flags = datos.get('vlp_sensibilidad', {})
    
    if vlp_flags.get('Hagedorn_Brown', 0) == 1:
        correlaciones_a_correr.append('Hagedorn_Brown')
    if vlp_flags.get('Beggs_Brill', 0) == 1:
        correlaciones_a_correr.append('Beggs_Brill')
    if vlp_flags.get('Ansari', 0) == 1:
        correlaciones_a_correr.append('Ansari')
    if vlp_flags.get('Gray', 0) == 1:
        correlaciones_a_correr.append('Gray')
    

    for nombre_corr in correlaciones_a_correr:
        print(f"  Calculando VLP para: {nombre_corr}...")
        
        datos_loop = copy.deepcopy(datos)
        datos_loop = _set_vlp_flags(datos_loop, nombre_corr)
        
        pwf_vlp_calculada = []
        for ql_bpd in q_eje_vlp:
            pwf = calcular_vlp_para_q(ql_bpd, datos_loop, Pb_yacimiento)
            pwf_vlp_calculada.append(pwf)
        
        q_vlp_validos = q_eje_vlp[np.isfinite(pwf_vlp_calculada)]
        p_vlp_validos = np.array(pwf_vlp_calculada)[np.isfinite(pwf_vlp_calculada)]
        
        curva_VLP_loop = list(zip(q_vlp_validos, p_vlp_validos))
        
        (q_op, p_op) = encontrar_interseccion(curva_IPR, curva_VLP_loop)
        
        lista_resultados_vlp.append({
            "nombre": nombre_corr,
            "curva_vlp": curva_VLP_loop,
            "punto_op": (q_op, p_op)
        })

    return {
        "ipr_usado": IPR_tipo,
        "curva_ipr": curva_IPR,
        "lista_vlps": lista_resultados_vlp,
        "RGL": glr,                            # <-- Necesario para la escala de gas
        "yacimiento_gas": es_yacimiento_gas    # <-- Bandera de Gas/Líquido
    }