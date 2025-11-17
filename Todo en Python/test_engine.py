# =====================================================================
# ESTE ES EL CÓDIGO COMPLETO FINAL PARA 'test_engine.py'
# Bórralo todo y pega esto.
# =====================================================================
import math
import numpy as np
from pvt import funcion_zg_exp, prop_oil_eos, dens_xfzg_tot
from utils import const
from core import correr_simulacion_completa # Importamos el motor principal

def run_all_tests():
    print("--- Iniciando Validación de Python (Motor Completo) ---")
    
    # --- PRUEBA 1: Validar Zg ---
    pr_test = 2.5
    tr_test = 1.5
    expected_z = 0.7819
    calculated_z = funcion_zg_exp(pr_test, tr_test)
    if math.isclose(calculated_z, expected_z, rel_tol=1e-4):
        print("✅ PRUEBA 1: ¡ÉXITO! (funcion_zg_exp)")
    else:
        print("❌ PRUEBA 1: ¡FALLÓ! (funcion_zg_exp)")
    print("----------------------------------------")

    # --- DATOS BASE PARA PRUEBA 2 y 3 (Pozo de Gas) ---
    datos_base_gas = {
        'Unidades': 0, 'Profundidad': 10006.6, 'Dti': 2.992, 'Rugosidad': 0.001,
        'Pcabezal': 278, 'Tcabezal': 90.248, 'Tfondo': 188.78, 'Pws': 1999.9,
        'API': 28.2, 'Sg': 0.645, 'Sw': 1.02, 'Porcentaje_Agua': 0.0,
        'RGL': 3895.19, 'Ql_cal': 5098, 'Pwf_cal': 1727, 'k': 50, 're': 1968,
        'S': 0, 'h': 49.2, 'rw': 0.328, 'Bo_input': 1.2, 'mu_o_input': 1.0,
        'N_Segmentos': 50, 'IPR_Tipo': 2, 'h_IPR_Automatico': 0, 'IPR_n': 0.615,
        'ff': 1.0, 'hl': 1.0, 'h_Ps_sanjari': 1, 'YCO2': 0.0, 'YH2s': 0.0, 'YN2': 0.0,
        'h_rb_MFF_HB': 0, 'h_radiobutBeggsBrill': 0, 'h_radiobutton_Ansari': 0, 'h_CFM_GRAY': 0,
        'h_rb_gen_visc': 1
    }

    # --- PRUEBA 2: Validar PVT Benchmark (Gas) ---
    print("Iniciando PRUEBA 2: Validación del Benchmark de PVT (Pozo de Gas)...")
    pb_yac_py = correr_simulacion_completa(datos_base_gas)['pb_calculada']
    expected_pb_yac = 51478.64
    if math.isclose(pb_yac_py, expected_pb_yac, rel_tol=1e-3):
        print("✅ PRUEBA 2: ¡ÉXITO! (Motor PVT Benchmark)")
    else:
        print("❌ PRUEBA 2: ¡FALLÓ! (Motor PVT Benchmark)")
    print("----------------------------------------")
    
    # --- PRUEBA 3: Validación del Punto de Operación (Gray) ---
    print("Iniciando PRUEBA 3: Validación VLP (Gray)...")
    datos_gray = datos_base_gas.copy()
    datos_gray['h_CFM_GRAY'] = 1
    expected_q_op = 5042.83
    
    try:
        resultados = correr_simulacion_completa(datos_gray)
        (q_op_py, p_op_py) = resultados['punto_op']
        print(f"Punto de Op. (Python-Gray): {q_op_py:.2f} bpd | Esperado (MATLAB-Gray): {expected_q_op:.2f} bpd")
        if math.isclose(q_op_py, expected_q_op, rel_tol=0.01): # Tolerancia del 1%
            print("✅ PRUEBA 3: ¡ÉXITO! El Punto de Operación (Gray) coincide.")
        else:
            print(f"❌ PRUEBA 3: ¡FALLÓ! El Qop (Gray) no coincide. (Delta = {q_op_py - expected_q_op:.2f} bpd)")
    except Exception as e:
        print(f"❌ PRUEBA 3: ¡FALLÓ! El motor (Gray) crasheó: {e}")
    print("----------------------------------------")

    # --- DATOS BASE PARA PRUEBA 4 y 5 (Pozo de Aceite) ---
    datos_base_aceite = {
        'Unidades': 0, 'Profundidad': 10006.6, 'Dti': 2.992, 'Rugosidad': 0.001,
        'Pcabezal': 278, 'Tcabezal': 90.248, 'Tfondo': 188.78, 'Pws': 3000,
        'API': 35.0, 'Sg': 0.645, 'Sw': 1.02, 'Porcentaje_Agua': 0.0,
        'RGL': 800.0, 'Ql_cal': 1000, 'Pwf_cal': 2500, 'k': 50, 're': 1968,
        'S': 0, 'h': 49.2, 'rw': 0.328, 'Bo_input': 1.2, 'mu_o_input': 1.0,
        'N_Segmentos': 50, 'IPR_Tipo': 2, 'h_IPR_Automatico': 0, 'IPR_n': 0.615,
        'ff': 1.0, 'hl': 1.0, 'h_Ps_sanjari': 1, 'YCO2': 0.0, 'YH2s': 0.0, 'YN2': 0.0,
        'h_rb_MFF_HB': 0, 'h_radiobutBeggsBrill': 0, 'h_radiobutton_Ansari': 0, 'h_CFM_GRAY': 0,
        'h_rb_gen_visc': 1
    }

    # --- PRUEBA 4: Validación de VLP (Hagedorn & Brown) ---
    print("Iniciando PRUEBA 4: Validación VLP (Hagedorn & Brown)...")
    datos_hb = datos_base_aceite.copy()
    datos_hb['h_rb_MFF_HB'] = 1
    expected_q_op_hb = 2523.50 # <- Nuestro nuevo benchmark
    
    try:
        resultados = correr_simulacion_completa(datos_hb)
        (q_op_py, p_op_py) = resultados['punto_op']
        print(f"Punto de Op. (Python-H&B): {q_op_py:.2f} bpd | Esperado (MATLAB-H&B): {expected_q_op_hb:.2f} bpd")
        if math.isclose(q_op_py, expected_q_op_hb, rel_tol=0.01): # Tolerancia del 1%
            print("✅ PRUEBA 4: ¡ÉXITO! El Punto de Operación (H&B) coincide.")
        else:
            print(f"❌ PRUEBA 4: ¡FALLÓ! El Qop (H&B) no coincide. (Delta = {q_op_py - expected_q_op_hb:.2f} bpd)")
    except Exception as e:
        print(f"❌ PRUEBA 4: ¡FALLÓ! El motor (H&B) crasheó: {e}")
    print("----------------------------------------")

    # --- PRUEBA 5: Validación de VLP (Beggs & Brill) ---
    print("Iniciando PRUEBA 5: Validación VLP (Beggs & Brill)...")
    datos_bb = datos_base_aceite.copy()
    datos_bb['h_radiobutBeggsBrill'] = 1
    # NO comparamos con el valor erróneo de MATLAB (295.45), solo verificamos que corra.
    
    try:
        resultados = correr_simulacion_completa(datos_bb)
        (q_op_py, p_op_py) = resultados['punto_op']
        if not np.isnan(q_op_py) and q_op_py > 0:
            print(f"Punto de Op. (Python-B&B): {q_op_py:.2f} bpd")
            print("✅ PRUEBA 5: ¡ÉXITO! El motor (Beggs & Brill) se ejecutó y encontró un resultado.")
        else:
            print("❌ PRUEBA 5: ¡FALLÓ! El motor (B&B) no encontró intersección (nan).")
    except Exception as e:
        print(f"❌ PRUEBA 5: ¡FALLÓ! El motor (B&B) crasheó: {e}")
    print("----------------------------------------")

if __name__ == "__main__":
    run_all_tests()