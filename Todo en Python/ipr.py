# Este es un nuevo archivo, guárdalo como 'ipr.py'
import numpy as np
from pvt import prop_oil_eos, dens_xfzg_tot
from utils import const

def calcular_ipr_lineal(datos, pws, pwf_cal, ql_cal, pwf_vec):
    """ Calcula IPR Lineal (Tipo 1) """
    if (pws - pwf_cal) < 1e-6:
        Prod = 1e6
    else:
        Prod = ql_cal / (pws - pwf_cal)
    q_vec = Prod * (pws - pwf_vec)
    return q_vec

def calcular_ipr_vogel(datos, pws, pwf_cal, ql_cal, pwf_vec):
    """ Calcula IPR Vogel (Tipo 2) """
    denom = (1 - 0.2 * (pwf_cal / pws) - 0.8 * (pwf_cal / pws)**2)
    if denom <= 1e-6:
        q_max = ql_cal * 1.5
    else:
        q_max = ql_cal / denom
    q_vec = q_max * (1 - 0.2 * (pwf_vec / pws) - 0.8 * (pwf_vec / pws)**2)
    return q_vec

def calcular_ipr_backpressure(datos, pws, pwf_cal, ql_cal, pwf_vec):
    """ Calcula IPR Backpressure (Tipo 3) """
    n = datos['IPR_n']
    if (pws**2 - pwf_cal**2) > 1e-6:
        cp_cal = ql_cal / (pws**2 - pwf_cal**2)**n
    else:
        cp_cal = 0
    q_vec = cp_cal * np.maximum(0, (pws**2 - pwf_vec**2))**n
    return q_vec

def calcular_ipr_darcy_gas(datos, pws, twf, api, gg, sg, glr, wc, YCO2, YH2s, YN2, Pb_yacimiento, pwf_vec):
    """ 
    Calcula IPR Darcy Gas p^2 (Tipo 4)
    Usa el bloque de cálculo de Z y mu_g del script de MATLAB
    """
    T_yac_R = twf + 459.67
    Zgas_yac = 0.9
    mu_g_yac = 0.02
    
    try:
        (Dens_AguaStd, Dens_AguaPT, DensHC_L, DensHC_G, Zv_eos, PM_eos, 
         Ymol_VTP, Rs_VTP, Tb_eos, Tpr_eos, Ppr_eos, Bo_VPT, Pb_eos) = prop_oil_eos(
             api, twf, pws, datos['h_Ps_sanjari'], gg, sg, Pb_yacimiento)
        
        RGA_m3 = glr / const['m3m3_a_scfSTB_orig']
        
        (Ymol_GL, Masa_Gas, Masa_total, Dens_aceite, Dens_gas, Dens_liquido, 
         Xo_masa, Xa_masa, Xg_masa, Z_GL, Zgas, PM_GLib, PMgas_tot, Pr_GL, Dens_GLstd) = dens_xfzg_tot(
             Zv_eos, YCO2, YH2s, YN2, Ymol_VTP, gg, api, PM_eos, twf, pws, 
             Rs_VTP, RGA_m3, wc, DensHC_L, DensHC_G, 0.0)

        Zgas_yac = Zgas
        
        T_yac_R = twf + 459.67
        Rho_g_pcm_yac = Dens_gas / 16.01846
        Rho_g_gcm3_yac = Rho_g_pcm_yac * 0.0160185
        X_visc_yac = 3.5 + (986 / T_yac_R) + (0.01 * PM_GLib)
        Y_visc_yac = 2.4 - (0.2 * X_visc_yac)
        K_visc_yac = (9.4 + 0.02 * PM_GLib) * (T_yac_R**1.5) / (209 + 19 * PM_GLib + T_yac_R)
        mu_g_yac = (K_visc_yac * 1e-4) * np.exp(X_visc_yac * (Rho_g_gcm3_yac**Y_visc_yac))

        if not np.isfinite(mu_g_yac) or mu_g_yac <= 0: mu_g_yac = 0.02
        if not np.isfinite(Zgas_yac) or Zgas_yac <= 0: Zgas_yac = 0.9

    except Exception as e:
        print(f"ADVERTENCIA (Python): Falla al calcular PVT de gas para IPR Tipo 4. {e}")

    # Calcular IPR Darcy Gas p^2
    k = datos['k']
    h = datos['h']
    re = datos['re']
    rw = datos['rw']
    S = datos['S']
    
    denominador_gas_darcy = 1422 * mu_g_yac * Zgas_yac * T_yac_R * (np.log(re / rw) + S)
    
    if denominador_gas_darcy < 1e-6:
        Cg_Mscfd = 1e6
    else:
        Cg_Mscfd = (k * h) / denominador_gas_darcy
        
    q_gas_Mscfd = Cg_Mscfd * np.maximum(0, (pws**2 - pwf_vec**2))
    
    if glr < 1e-3:
        q_vec = 0 * q_gas_Mscfd
    else:
        q_vec = (q_gas_Mscfd * 1000) / glr
        
    return q_vec