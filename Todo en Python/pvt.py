# =====================================================================
# ESTE ES EL CÓDIGO COMPLETO PARA 'pvt.py'
# Bórralo todo y pega esto.
# =====================================================================

import numpy as np
from utils import const # Importamos las constantes

def funcion_zg_exp(pr_g, tr_g):
    """
    Traducción de la función de MATLAB Funcion_ZgExp
    Estimación Explícita de Z (Correlación original del HBC.m)
    """
    # Prevenir división por cero o logs de cero si tr_g es inválido
    if not isinstance(tr_g, (int, float, np.ndarray)): tr_g = float(tr_g)
    if not isinstance(pr_g, (int, float, np.ndarray)): pr_g = float(pr_g)
    
    tr_g = np.clip(tr_g, 1e-6, None) # Asegura que tr_g no sea cero
        
    # Constantes
    Za1 = 0.317842
    Za2 = 0.382216
    Za3 = -7.76835
    Za4 = 14.2905
    Za5 = 2.18363e-06
    Za6 = -0.00469257
    Za7 = 0.0962541
    Za8 = 0.16672
    Za9 = 0.96691
    Za10 = 0.063069
    Za11 = -1.966847
    Za12 = 21.0581
    Za13 = -27.0246
    Za14 = 16.23
    Za15 = 207.783
    Za16 = -488.161
    Za17 = 176.29
    Za18 = 1.88453
    Za19 = 3.05921
    
    # Cálculos
    inv_tr_g = 1 / tr_g 
    
    ZA = Za1 * inv_tr_g * np.exp(Za2 * (1 - inv_tr_g)**2) * pr_g
    ZB = Za3 * inv_tr_g + Za4 * inv_tr_g**2 + Za5 * inv_tr_g**6 * pr_g**6
    ZC = Za9 + Za8 * pr_g * inv_tr_g + Za7 * inv_tr_g**2 * pr_g**2 + Za6 * inv_tr_g**3 * pr_g**3
    ZD = Za10 * inv_tr_g * np.exp(Za11 * (1 - inv_tr_g)**2)
    ZE = Za12 * inv_tr_g + Za13 * inv_tr_g**2 + Za14 * inv_tr_g**3
    ZF = Za15 * inv_tr_g + Za16 * inv_tr_g**2 + Za17 * inv_tr_g**3
    ZG = Za18 + Za19 * inv_tr_g
    
    # Denominador de yZ (protegido contra división por cero)
    den_yZ = (1 + ZA**2) / ZC - (ZA**2 * ZB) / ZC**3
    den_yZ = np.where(den_yZ == 0, 1e-9, den_yZ)
    
    yZ = ZD * pr_g / den_yZ
    
    # Denominador final (protegido contra división por cero)
    den_final = (ZD * pr_g + ZE * yZ**2 - ZF * (yZ + 0j)**ZG) * (1 - yZ)**3
    den_final = np.where(den_final == 0, 1e-9, den_final)
    
    zg_explicit = ZD * pr_g * (1 + yZ + yZ**2 - yZ**3) / den_final
    
    # Manejar NaNs o Infs y devolver solo la parte real (como en MATLAB)
    zg_explicit = np.real(zg_explicit)
    
    # np.where es el equivalente de MATLAB (condicion, valor_si_true, valor_si_false)
    zg_explicit = np.where(np.isfinite(zg_explicit), zg_explicit, 0.9)
        
    return float(zg_explicit)


# =====================================================================
# REEMPLAZA ESTA FUNCIÓN EN 'pvt.py'
# =====================================================================

def prop_oil_eos(api, temp_F, press_psi, Pb_eos_OPtion, gg, sg, Pb_yacimiento):
    """
    Traducción de la función de MATLAB PropOil_EOS
    CORREGIDO: Paréntesis en FNC_eos
    """
    # Inicializar
    Zv_eos = 0.0
    Zl_eos = 0.0
    Pb_eos = 0.0
    
    # =====================================================================
    # [INICIO DE CORRECCIÓN] Paréntesis de FNC_eos
    # La fórmula es (A / B) ** 10
    termino_A = np.log(1.07 - (141.5 / (131.5 + api))) - 3.56073
    termino_B = -2.93886
    FNC_eos = (termino_A / termino_B)**(1 / 0.1)
    # [FIN DE CORRECCIÓN]
    # =====================================================================
    
    Tb_eos = (1080 - np.exp(6.97996 - 0.01964 * FNC_eos**(2/3))) * 9 / 5
    SG_eos = 1.07 - np.exp(3.56073 - 2.93886 * FNC_eos**0.1)
    PM_eos = 0.000045673 * Tb_eos**2.1962 * SG_eos**-1.0164
    Nc_eos = (FNC_eos + 4) / 14.5
    press_psi = np.real(press_psi)
    
    if Nc_eos <= 7:
        W_eos = 0.000000018606083 * FNC_eos**3 - 0.0000110073908 * FNC_eos**2 + 0.00432004994 * FNC_eos - 0.0366640594
    else:
        W_eos = -(0.3 - np.exp(-6.252 + 3.64457 * FNC_eos**0.1))
        
    if Nc_eos < 17 or Nc_eos > 50:
        Pc_eos = np.real((0.032688 + 0.000385 * FNC_eos)**-1.25) * 14.50377
    else:
        Pc_eos = np.real(-(-0.64 - np.exp(6.34492 - 0.7239 * FNC_eos**0.3))) * 14.50377
        
    if Nc_eos > 50:
        Tc_eos = (3 / 7 * Tb_eos * np.log10(Pc_eos / 1.01325) / (W_eos + 1)) * 9 / 5 + Tb_eos
    else:
        Tb_TcR = 1.2 - np.exp(-0.34742 - 0.02327 * FNC_eos**0.55)
        Tc_eos = Tb_eos / Tb_TcR
        
    Zc_eos = 0.29056 - 0.08775 * W_eos
    
    if W_eos >= 0.491:
        Kw_eos = 0.378893464 + 1.4897153 * W_eos - 0.1713184 * W_eos**2 + 0.0196554 * W_eos**3
    else:
        Kw_eos = 0.37464 + 1.54226 * W_eos - 0.26992 * W_eos**2
        
    ACvL_eos = 0.077796074 * 83.14472 * (Pc_eos / 14.50377) / (Tc_eos * 5 / 9) * (1 - 2.258 * PM_eos**-0.1823)
    Tpr_eos = (temp_F + 459.67) / Tc_eos
    Ppr_eos = (press_psi + 14.697) / Pc_eos
    
    Zv_expl = funcion_zg_exp(Ppr_eos, Tpr_eos)
    
    fo_Tr = 14.7114 - 6.7632 / Tpr_eos + 26.5948 * Tpr_eos - 34.5428 * Tpr_eos**0.8
    f1_Tr = 49.1821 - 14.6979 / Tpr_eos + 87.9972 * Tpr_eos - 122.495 * Tpr_eos**0.8
    f2_Tr = 6.6828 - 1.8259 / Tpr_eos + 7.8256 * Tpr_eos - 12.7191 * Tpr_eos**0.8
    Pvr_Tar_San = np.exp(fo_Tr + W_eos * f1_Tr + W_eos**2 * f2_Tr)
    Pvap_Tar_San = Pc_eos * Pvr_Tar_San
    
    fo_LeeKes = 5.92714 - 6.09648 / Tpr_eos - 1.28862 * np.log(Tpr_eos) + 0.169347 * Tpr_eos**6
    f1_LeeKes = 15.2518 - 15.6875 / Tpr_eos - 13.4721 * np.log(Tpr_eos) + 0.43577 * Tpr_eos**6
    Pvr_LeeKes = np.exp(fo_LeeKes + W_eos * f1_LeeKes)
    Pvap_LeeKes = Pc_eos * Pvr_LeeKes
    
    Pb_eos = Pvap_LeeKes if Pb_eos_OPtion == 1 else Pvap_Tar_San
    
    Alfa_eos = (1 + Kw_eos * (1 - Tpr_eos**0.5))**2
    Omega_a = 0.66121 - 0.76106 * Zc_eos
    Omega_b = 0.02207 + 0.20868 * Zc_eos
    K1_eos = (0.57765 - 1.8708 * Zc_eos) / Omega_b + 1
    K2_eos = (1.8708 * Zc_eos - 0.57765) / Omega_b
    A_eos = Omega_a * Alfa_eos * Ppr_eos / Tpr_eos**2
    B_eos = Omega_b * Ppr_eos / Tpr_eos
    a1_vpt = B_eos * (K1_eos - 1) - 1
    a2_vpt = B_eos**2 * (K2_eos - K1_eos) - K1_eos * B_eos + A_eos
    a3_vpt = -B_eos * (K2_eos * B_eos**2 + K2_eos * B_eos + A_eos)
    Poli_Z = [1, a1_vpt, a2_vpt, a3_vpt]
    
    if any(np.isnan(Poli_Z)) or any(np.isinf(Poli_Z)):
        Zl_eos = 0.0
        Zv_eos = 0.0
    else:
        Z_root = np.roots(Poli_Z)
        RaicesRealZ = Z_root[np.isreal(Z_root)].real
        RaicesZ = RaicesRealZ[RaicesRealZ > 0]
        if RaicesZ.size > 0:
            if RaicesZ.size == 1:
                if press_psi > Pb_eos:
                    Zl_eos = RaicesZ[0]
                    Zv_eos = 0.0
                else:
                    Zl_eos = 0.0
                    Zv_eos = Zv_expl
            else:
                Zl_eos = np.min(RaicesZ)
                Zv_eos = np.max(RaicesZ)
                
    if Zv_eos <= 0:
        KeLV = 0.0
    elif Zv_eos >= 1:
        KeLV = Pb_eos / press_psi
    else:
        KeLV = 1 / Ppr_eos * np.exp(5.371 * (1 + W_eos) * (1 - 1 / Tpr_eos))
        
    Xmol_VTP = 1 / (1 + KeLV)
    Ymol_VTP = 1 - Xmol_VTP # Corrección de sintaxis
    
    DensHC_G = 0.0
    if Zv_eos > 0:
        Vg_vtp = 83.14472 * Zv_eos * ((temp_F + 459.67) * 5 / 9) / ((press_psi + 14.697) / 14.50377)
        if Vg_vtp != 0: DensHC_G = PM_eos / Vg_vtp
        
    VL_Sat_Ratch = 83.14474 * (Tc_eos * 5 / 9) / (Pc_eos / 14.50377) * (0.29056 - 0.08775 * W_eos)**(1 + (1 - Tpr_eos)**(2/7))
    Dens_AguaStd = -0.0036 * 14.697 + 999.97495 * (1 - (15.556 - 3.983035)**2 * (15.556 + 301.797) / (522528.9 * (15.556 + 69.34881)))
    TCent = (temp_F - 32) * 5 / 9
    Dens_AguaPT = -0.0036 * (press_psi + 14.697) + 999.97495 * (1 - (TCent - 3.983035)**2 * (TCent + 301.797) / (522528.9 * (TCent + 69.34881)))
    
    VL_vtp = 0.0
    if Zl_eos > 0:
        VL_vtp = 83.14472 * Zl_eos * ((temp_F + 459.67) * 5 / 9) / ((press_psi + 14.697) / 14.50377) - ACvL_eos
        
    VL_SatHC = (VL_Sat_Ratch + A_eos * Omega_a) * Xmol_VTP + VL_vtp * (1 - Xmol_VTP)
    
    DensHC_L = 0.0
    if VL_SatHC != 0: DensHC_L = PM_eos / VL_SatHC
    
    # --- Lógica de Régimen Saturado/Subsaturado ---
    X_standing = (0.0125 * api - 0.00091 * ((temp_F - 32) * 5 / 9))
    F_standing = 10**(X_standing * 1.2048) # Equivale a (10**X)**1.2048

    Rsb_StandingFnc = gg * (Pb_yacimiento / 18.2 + 1.4) * F_standing
    Bo_StandingFnc = 0.9759 + 0.00012 * (Rsb_StandingFnc * (gg / sg)**0.5 + 1.25 * temp_F)**1.2
    
    if press_psi <= Pb_yacimiento:
        # --- Régimen Saturado (P <= Pb) ---
        Rs_VTP = gg * (press_psi / 18.2 + 1.4) * F_standing
        Bo_VPT = 0.9759 + 0.00012 * (Rs_VTP * (gg / sg)**0.5 + 1.25 * temp_F)**1.2
    else:
        # --- Régimen Subsaturado (P > Pb) ---
        Rs_VTP = Rsb_StandingFnc
        Bo_VPT = Bo_StandingFnc
        
    # Fallbacks
    if not np.isfinite(DensHC_L): DensHC_L = 800.0
    if not np.isfinite(DensHC_G): DensHC_G = 1.0
    if not np.isfinite(Zv_eos): Zv_eos = 0.9
    if not np.isfinite(PM_eos): PM_eos = 20.0
    if not np.isfinite(Rs_VTP): Rs_VTP = 100.0
    if not np.isfinite(Bo_VPT): Bo_VPT = 1.1
    if not np.isfinite(Pb_eos): Pb_eos = 0.0
        
    return (Dens_AguaStd, Dens_AguaPT, DensHC_L, DensHC_G, Zv_eos, PM_eos, 
            Ymol_VTP, Rs_VTP, Tb_eos, Tpr_eos, Ppr_eos, Bo_VPT, Pb_eos)


def dens_xfzg_tot(Zv_eos, YCO2, YH2s, YN2, Ymol_VTP, gg, api, PM_eos, 
                  temp_F, press_psi, Rs_VTP, RGA_m3, wc, DensHC_L, 
                  DensHC_G, Dens_AguaPT):
    """
    Traducción de la función de MATLAB DensXfZg_tot
    """
    YGLpuro = 1 - (YCO2 + YH2s + YN2)
    Xmol_VTP = 1 - Ymol_VTP
    PM_VTP = Xmol_VTP * PM_eos + 28.964 * Ymol_VTP * gg
    Ymas_VTP = Ymol_VTP * (28.964 * gg / PM_VTP)
    Rs_m3 = Rs_VTP / const['m3m3_a_scfSTB_orig']
    
    Epsilon_GL = 107.6 * ((YCO2 + YH2s) - (YCO2 + YH2s)**2.2) + 5.9 * (YH2s**0.06 - YH2s**0.68)
    
    if gg >= 1.862 or api > 45:
        Pc_GL = 744 - 125.4 * gg + 5.9 * gg**2
        Tc_GL = 164.3 + 357.7 * gg - 67.7 * gg**2
    else:
        Pc_GL = 671.1 + 14 * gg - 34.3 * gg**2
        Tc_GL = 120.1 + 429 * gg - 62.9 * gg**2
        
    Pc_GL = (Pc_GL * YGLpuro + 1073 * YCO2 + 1307 * YH2s + 492 * YN2) * Tc_GL / ((Tc_GL + Epsilon_GL) + YH2s * (1 - YH2s) * Epsilon_GL)
    Tc_GL = Tc_GL * YGLpuro + 547.7 * YCO2 + 672.4 * YH2s + 227.1 * YN2 - Epsilon_GL
    
    Pr_GL = (press_psi + 14.697) / Pc_GL
    Tr_GL = (temp_F + 459.67) / Tc_GL
    
    # =====================================================================
    # Aquí es donde Pylance daba el error. Ahora 'funcion_zg_exp' sí existe.
    Z_GL = funcion_zg_exp(Pr_GL, Tr_GL)
    
    Pr_GLstd = 14.697 / Pc_GL
    Tr_GLstd = (60 + 459.67) / Tc_GL
    
    Z_GLstd = funcion_zg_exp(Pr_GLstd, Tr_GLstd)
    # =====================================================================
    
    PM_GLib = gg * 28.964 * YGLpuro + YCO2 * 44.01 + YH2s * 34.08 + YN2 * 28.02
    Dens_GLstd = PM_GLib * 14.697 / (Z_GLstd * 10.732 * (459.67 + 60)) * 16.018463
    
    Vol_GL = RGA_m3 - Rs_m3 if RGA_m3 > Rs_m3 else 0.0
    Masa_GL = Dens_GLstd * Vol_GL
    Dens_GL = PM_GLib * (press_psi + 14.697) / (Z_GL * 10.732 * (459.67 + temp_F)) * 16.018463
    
    den_Dens_gas = Ymas_VTP / PM_eos + Masa_GL / PM_GLib
    if den_Dens_gas == 0: den_Dens_gas = 1e-9
    
    Dens_gas = DensHC_G * 1000 * (Ymas_VTP / PM_eos / den_Dens_gas) + Dens_GL * (Masa_GL / PM_GLib / den_Dens_gas)
    Dens_aceite = DensHC_L * 1000
    
    den_Masa_Agua = 1 - wc / 100
    if den_Masa_Agua == 0: den_Masa_Agua = 1e-9
    
    Masa_Agua = (1 / den_Masa_Agua - 1) * (-0.0036 * 14.697 + 999.97495 * (1 - (15.556 - 3.983035)**2 * (15.556 + 301.797) / (522528.9 * (15.556 + 69.34881))))
    Masa_oil = DensHC_L * 1000
    Masa_Gas = DensHC_G * 1000 * Rs_m3 + Masa_GL
    
    den_Moles_gasTOt_1 = PM_VTP if PM_VTP != 0 else 1e-9
    den_Moles_gasTOt_2 = PM_GLib if PM_GLib != 0 else 1e-9
    Moles_gasTOt = DensHC_G * 1000 * Rs_m3 / den_Moles_gasTOt_1 + Masa_GL / den_Moles_gasTOt_2
    
    Ymol_GL = 0.0
    if Masa_GL > 0 and Moles_gasTOt != 0:
        Ymol_GL = (Masa_GL / PM_GLib) / Moles_gasTOt
        
    Masa_total = Masa_Agua + Masa_oil + Masa_Gas
    Zgas = Z_GL * Ymol_GL + (1 - Ymol_GL) * Zv_eos
    
    Xg_masa = 0.0
    Xo_masa = 0.0
    Xa_masa = 0.0
    if Masa_total != 0:
        Xg_masa = Masa_Gas / Masa_total
        Xo_masa = Masa_oil / Masa_total
        Xa_masa = Masa_Agua / Masa_total
        
    Dens_liquido = Dens_aceite
    if (Xo_masa + Xa_masa) != 0:
        Dens_liquido = (Dens_aceite * Xo_masa + Dens_AguaPT * Xa_masa) / (Xo_masa + Xa_masa)
        
    PMgas_tot = 0.0
    if Moles_gasTOt != 0:
        PMgas_tot = Masa_Gas / Moles_gasTOt
        
    # Fallbacks
    if not np.isfinite(Dens_gas): Dens_gas = 1.0
    if not np.isfinite(Dens_liquido): Dens_liquido = 800.0
    if not np.isfinite(Z_GL): Z_GL = 0.9
    if not np.isfinite(Zgas): Zgas = 0.9
    if not np.isfinite(PM_GLib): PM_GLib = 20.0
    if not np.isfinite(PMgas_tot): PMgas_tot = 20.0
        
    return (Ymol_GL, Masa_Gas, Masa_total, Dens_aceite, Dens_gas, Dens_liquido, 
            Xo_masa, Xa_masa, Xg_masa, Z_GL, Zgas, PM_GLib, PMgas_tot, Pr_GL, Dens_GLstd)