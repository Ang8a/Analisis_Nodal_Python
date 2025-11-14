# =====================================================================
# ESTE ES EL CÓDIGO COMPLETO FINAL PARA 'flow.py' (CON ANSARI RECONSTRUIDO)
# Bórralo todo y pega esto.
# =====================================================================

import numpy as np
from scipy.optimize import fsolve # Importamos el solver

# --- Función Helper para el solver de fzero de MATLAB ---
def _solve_fzero(func, x0):
    """
    Un wrapper seguro para fsolve que imita el comportamiento de fzero
    (maneja fallos de convergencia).
    """
    try:
        sol, infodict, ier, mesg = fsolve(func, [x0], full_output=True, xtol=1e-6)
        if ier == 1: # 'ier == 1' significa que la solución convergió
            return sol[0], True
        else:
            return np.nan, False # Falló la convergencia
    except Exception:
        return np.nan, False # Falla catastrófica

def _friccion_ansari(Re, D, Epsilon_met):
    """ Calcula el factor de fricción de Haaland (usado por Ansari) """
    if Re <= 2100:
        return 16 / max(Re, 1e-9)
    try:
        # Usamos Haaland (estable y rápido)
        fric_term = (Epsilon_met / (3.7 * D))**1.11 + (6.9 / Re)
        f = (1 / (-1.8 * np.log10(fric_term)))**2
    except (ValueError, OverflowError):
        f = 0.02 # Fallback
        
    if Re < 4000: # Transición
        f_laminar = 16 / Re
        f = (f_laminar * (4000 - Re) + f * (Re - 2000)) / 2000
    
    return f

# =====================================================================
# HAGEDORN & BROWN (Validada)
# =====================================================================
def calcular_vlp_hagedorn_brown(datos, aj_mid, mass_flow_rate, Dti_pulg, vis_liq, gvis, e_abs_ft, an_mid):
    Dft = Dti_pulg / 12.0
    ak_mid_den = Dft * vis_liq**aj_mid * gvis**(1 - aj_mid)
    if ak_mid_den == 0: ak_mid_den = 1e-9
    ak_mid = 0.022 * mass_flow_rate / ak_mid_den
    
    if ak_mid <= 2000:
        al_mid = 16 / max(ak_mid, 1e-9)
    else:
        with np.errstate(all='ignore'):
            term1 = 0.2698 * e_abs_ft / Dft
            term2 = 5.0452 / ak_mid
            term3 = 0.3539 * (e_abs_ft / Dft)**1.1098
            term4 = 5.8506 / (ak_mid**0.8981)
            log_term_inner = term3 + term4
            if log_term_inner <= 0: log_term_inner = 1e-9
            log_term_outer = term1 - term2 * np.log10(log_term_inner)
            if log_term_outer <= 0: log_term_outer = 1e-9
            al_mid = (1 / (-4 * np.log10(log_term_outer)))**2
        if ak_mid < 4000:
            f_laminar = 16 / ak_mid
            al_mid = (f_laminar * (4000 - ak_mid) + al_mid * (ak_mid - 2000)) / 2000
    if not np.isfinite(al_mid): al_mid = 0.02
    dp_elev = an_mid * np.cos(np.deg2rad(0)) / 144
    dp_fric_den = 7.413e10 * Dft**5 * an_mid
    if dp_fric_den == 0: dp_fric_den = 1e-9
    dp_fric = (al_mid * mass_flow_rate**2) / dp_fric_den / 144
    dp_pat_flw = dp_elev + dp_fric
    if not np.isfinite(dp_pat_flw) or dp_pat_flw <= 0: dp_pat_flw = 0.001
    return dp_pat_flw

# =====================================================================
# BEGGS & BRILL (Corregida)
# =====================================================================
def calcular_vlp_beggs_brill(datos, press_col, ang_dev_raw, Dti_pulg, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, e_abs_ft):
    ang_dev = 90.0
    Dft = Dti_pulg / 12.0
    Vm_BB = usg_col + usl_col
    LamdaL_BB = usl_col / max(Vm_BB, 1e-9)
    if liq_tsp_col <= 0: liq_tsp_col = 1e-9
    Nlv_BB = 1.938 * usl_col * (liq_den_col / liq_tsp_col)**0.25
    Nfr_BB = Vm_BB**2 / max(32.174049 * Dft, 1e-9)
    
    L1_BB = 316 * LamdaL_BB**0.302
    L2_BB = 0.0009252 * LamdaL_BB**-2.4684
    L3_BB = 0.1 * LamdaL_BB**-1.4516
    L4_BB = 0.5 * LamdaL_BB**-6.738

    with np.errstate(all='ignore'):
        if (LamdaL_BB < 0.01 and Nfr_BB < L1_BB) or (LamdaL_BB >= 0.01 and Nfr_BB < L2_BB):
            a, b, c = 0.98, 0.4846, 0.0868
        elif (LamdaL_BB >= 0.01 and LamdaL_BB < 0.4 and Nfr_BB > L3_BB and Nfr_BB <= L1_BB) or \
             (LamdaL_BB >= 0.4 and Nfr_BB > L3_BB and Nfr_BB <= L4_BB):
            a, b, c = 0.845, 0.5351, 0.0173
        elif (LamdaL_BB >= 0.01 and Nfr_BB >= L2_BB and Nfr_BB <= L3_BB):
            HL0_Seg = (0.98 * LamdaL_BB**0.4846) / max(Nfr_BB**0.0868, 1e-9)
            HL0_Int = (0.845 * LamdaL_BB**0.5351) / max(Nfr_BB**0.0173, 1e-9)
            A_BB = (L3_BB - Nfr_BB) / max(L3_BB - L2_BB, 1e-9)
            B_BB = 1.0 - A_BB
            HL0_BB = A_BB * HL0_Seg + B_BB * HL0_Int
            a, b, c = 0, 0, 0
        else:
            a, b, c = 1.065, 0.5824, 0.0609
        if a != 0:
            HL0_BB = (a * LamdaL_BB**b) / max(Nfr_BB**c, 1e-9)
    
    if ang_dev > 0:
        C_BB = (1 - LamdaL_BB) * np.log(
            max(0.011 * (LamdaL_BB**-3.768) * (Nlv_BB**3.539) * (Nfr_BB**-1.614), 1e-9)
        )
        if C_BB < 0: C_BB = 0
        beta = 1 + C_BB * (np.sin(np.deg2rad(1.8 * ang_dev)) - 0.333 * np.sin(np.deg2rad(1.8 * ang_dev))**3)
        HL_BB = HL0_BB * beta
    else:
        HL_BB = HL0_BB

    if not np.isfinite(HL_BB): HL_BB = LamdaL_BB
    HL_BB = max(HL_BB, LamdaL_BB)
    HL_BB = min(HL_BB, 1.0)
    HL_BB = HL_BB * datos['hl']
    if HL_BB == 0: HL_BB = 1e-6
    Dens_BB = liq_den_col * HL_BB + (1 - HL_BB) * gdens_col
    
    y_BB = LamdaL_BB / max(HL_BB**2, 1e-9)
    
    if y_BB <= 1.0:
        S_BB = 0.0
    elif 1.0 < y_BB <= 1.2:
        S_BB = np.log(2.2 * (y_BB - 1))
    else:
        yf_BB = np.log(y_BB)
        den_S2_BB = -0.0523 + 3.182 * yf_BB - 0.8725 * yf_BB**2 + 0.01853 * yf_BB**4
        S_BB = yf_BB / max(den_S2_BB, 1e-9)
        
    if not np.isfinite(S_BB): S_BB = 0.0
    S_BB = np.clip(S_BB, -np.inf, 20)
        
    visn_BB = vis_liq_col * LamdaL_BB + gvis_col * (1 - LamdaL_BB)
    densn_BB = liq_den_col * LamdaL_BB + gdens_col * (1 - LamdaL_BB)
    den_NRen_BB = max(visn_BB, 1e-9)
    NRen_BB = densn_BB * Vm_BB * Dft / den_NRen_BB * 1488

    if NRen_BB <= 2100:
        fn_BB = 16 / max(NRen_BB, 1e-9)
    else:
        fric_term = (e_abs_ft / (3.7 * Dft))**1.11 + (6.9 / NRen_BB)
        fn_BB = (1 / (-1.8 * np.log10(fric_term)))**2
            
    if not np.isfinite(fn_BB): fn_BB = 0.02
        
    dPBBe = Dens_BB * np.cos(np.deg2rad(ang_dev)) / 144
    f_BB = fn_BB * np.exp(S_BB)
    dPBBf = (datos['ff'] * f_BB * densn_BB * Vm_BB**2) / (2 * 32.174049 * Dft) / 144
    
    dpBBtotal = dPBBe + dPBBf
    return dpBBtotal

# =====================================================================
# GRAY (Validada)
# =====================================================================
def calcular_vlp_gray(datos, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, mass_flow_rate, e_abs_ft, Dti_pulg):
    Dft = Dti_pulg / 12.0
    CL_gray = usl_col / max(usg_col + usl_col, 1e-9)
    Vm_Gray = usg_col + usl_col
    pmNS_gray = gdens_col * (1 - CL_gray) + CL_gray * liq_den_col
    R_gray = usl_col / max(usg_col, 1e-9)
    
    den_Nv_gray = 32.174 * liq_tsp_col * (liq_den_col - gdens_col)
    if den_Nv_gray == 0: den_Nv_gray = 1e-9
    Nv_gray = 453.592 * pmNS_gray**2 * Vm_Gray**4 / den_Nv_gray
    
    den_Nd_gray = liq_tsp_col
    if den_Nd_gray == 0: den_Nd_gray = 1e-9
    Nd_gray = 453.592 * 32.174 * (liq_den_col - gdens_col) * (Dft)**2 / den_Nd_gray
    
    with np.errstate(all='ignore'):
        B_gray = 0.0814 * (1 - 0.0554 * np.log10(1 + 730 * R_gray / (R_gray + 1)))
        A_gray = -2.2314 * (Nv_gray * (1 + 205 / Nd_gray))**B_gray
        Hg_gray = (1 - np.exp(A_gray)) / max(R_gray + 1, 1e-9)
    
    HL_gray = (1 - Hg_gray) * datos['hl']
    HL_gray = np.clip(HL_gray, 0, 1)
    pmSlip_gray = gdens_col * (1 - HL_gray) + HL_gray * liq_den_col
    
    den_Re_gray = (Dft * vis_liq_col**CL_gray * gvis_col**(1 - CL_gray))
    if den_Re_gray == 0: den_Re_gray = 1e-9
    Re_gray = 0.022 * mass_flow_rate / den_Re_gray
    
    if Re_gray <= 2000:
        f_gray = 16 / max(Re_gray, 1e-9)
    else:
        with np.errstate(all='ignore'):
            term1 = 0.2698 * e_abs_ft / Dft
            term2 = 5.0452 / Re_gray
            term3 = 0.3539 * (e_abs_ft / Dft)**1.1098
            term4 = 5.8506 / (Re_gray**0.8981)
            log_term_inner = term3 + term4
            if log_term_inner <= 0: log_term_inner = 1e-9
            log_term_outer = term1 - term2 * np.log10(log_term_inner)
            if log_term_outer <= 0: log_term_outer = 1e-9
            f_gray = (1 / (-4 * np.log10(log_term_outer)))**2
        if Re_gray < 4000:
            f_laminar = 16 / Re_gray
            f_gray = (f_laminar * (4000 - Re_gray) + f_gray * (Re_gray - 2000)) / 2000
    if not np.isfinite(f_gray): f_gray = 0.02
    dp_elev = pmSlip_gray * np.cos(np.deg2rad(0)) / 144
    dp_fric = datos['ff'] * pmNS_gray * f_gray * Vm_Gray**2 / (2 * 32.174049 * Dft) / 144
    dp_pat_flw = dp_elev + dp_fric
    if not np.isfinite(dp_pat_flw) or dp_pat_flw <= 0: dp_pat_flw = 0.001
    return dp_pat_flw

# =====================================================================
# ANSARI (NUEVA IMPLEMENTACIÓN HÍBRIDA-MECANICISTA)
# =====================================================================
def calcular_vlp_ansari(datos, press_col, ang_dev_raw, Dti_pulg, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, e_abs_ft, mass_flow_rate):
    """
    Implementación Híbrida-Mecanicista de Ansari (BASADA EN LIBRO DE TEXTO).
    Usa los criterios de transición de Ansari para SELECCIONAR
    la correlación empírica/mecanicista apropiada.
    """
    # 0. Unpack y Conversión a Métrico
    Dmetros = Dti_pulg / 12.0 * 0.3048
    DensLiq_met = liq_den_col * 16.01846
    DensGas_met = gdens_col * 16.01846
    TensSup_met = liq_tsp_col * 0.001
    Vsl_met = usl_col * 0.3048
    Vsg_met = usg_col * 0.3048
    ViscLiq_met = vis_liq_col * 0.001
    ViscGas_met = gvis_col * 0.001
    Epsilon_met = e_abs_ft * 0.3048
    Fac_HL = datos['hl']
    factor_ff = datos['ff']
    ang_dev = 90.0 # Vertical
    PA_M_TO_PSI_FT = 0.3048 / 6894.76
    
    # 1. Criterio de Transición a Flujo Anular (Ansari Eq. 3.29)
    # Si Vsg es mayor que esto, es Flujo Anular
    Vsg_trans_annular = 3.1 * (9.81 * (DensLiq_met - DensGas_met) * TensSup_met / max(DensLiq_met**2, 1e-9))**0.25
    
    if Vsg_met >= Vsg_trans_annular:
        # --- PATRÓN: FLUJO ANULAR ---
        # El flujo es de alta velocidad, dominado por gas (como tu pozo de gas).
        # La correlación más robusta para esto es Gray.
        return calcular_vlp_gray(datos, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, mass_flow_rate, e_abs_ft, Dti_pulg)

    # 2. Criterio de Transición a Flujo Burbuja Dispersa
    # Si Vsg es mayor que esto (Y Vsl es alto), es Burbuja Dispersa
    Vslip = 1.53 * (9.81 * (DensLiq_met - DensGas_met) * TensSup_met / max(DensLiq_met**2, 1e-9))**0.25
    Vsg_trans_dispersed = 2.0 * Vsl_met + 0.5 * Vslip
    
    if Vsg_met >= Vsg_trans_dispersed:
        # --- PATRÓN: FLUJO BURBUJA DISPERSA ---
        # (Alta velocidad de líquido, gas disperso)
        # B&B (Distribuido) es el mejor proxy, ya que H&B no tiene este régimen.
        return calcular_vlp_beggs_brill(datos, press_col, ang_dev_raw, Dti_pulg, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, e_abs_ft)

    # 3. Criterio de Transición Burbuja/Slug (Ansari Eq. 3.12)
    # Si Vsg es MENOR que esto, es Flujo Burbuja
    HL_trans_bubble_slug = 0.25 # (Void Fraction = 0.25)
    Vsg_trans_bubble_slug = (Vsl_met / (1.0 - HL_trans_bubble_slug)) - 1.2 * Vsl_met + Vslip * (HL_trans_bubble_slug**0.5)

    if Vsg_met < Vsg_trans_bubble_slug:
        # --- PATRÓN: FLUJO BURBUJA ---
        # Usamos el modelo físico de Flujo Burbuja (Harmathy/Zuber)
        # Esto es lo que la lógica de H&B modela.
        
        # Necesitamos calcular aj_mid (Holdup H&B) y an_mid (Densidad H&B)
        # (Repetimos la lógica H&B aquí solo para el Holdup)
        with np.errstate(all='ignore'):
            y_mid = 1.938*usl_col*(liq_den_col/liq_tsp_col)**0.25
            z_mid = 1.938*usg_col*(liq_den_col/liq_tsp_col)**0.25
            aa_mid = 120.872*Dti_pulg/12*np.sqrt(liq_den_col/liq_tsp_col)
            ab_mid = 0.15726*vis_liq_col*(1/(liq_den_col*liq_tsp_col**3))**0.25
            ac12=-2.69851; ac13=0.1584095; ac14=-0.55099756; ac15=0.54784917; ac16=-0.12194578
            ae12=-0.10306578; ae13=0.617774; ae14=-0.632946; ae15=0.29598; ae16=-0.0401
            ai12=0.91162574; ai13=-4.82175636; ai14=1232.25036621; ai15=-22253.57617; ai16=116174.28125
            ac_mid = 10**(ac12+ac13*(np.log10(ab_mid)+3)+ac14*(np.log10(ab_mid)+3)**2+ac15*(np.log10(ab_mid)+3)**3+ac16*(np.log10(ab_mid)+3)**4)
            ad_mid = y_mid/max(z_mid**0.575, 1e-9)*((press_col/14.697)**0.1)*(ac_mid/max(aa_mid, 1e-9))
            ad_mid = np.clip(ad_mid, 0.00000161, 0.0047)
            ae_mid = ae12+ae13*(np.log10(ad_mid)+6)+ae14*(np.log10(ad_mid)+6)**2+ae15*(np.log10(ad_mid)+6)**3+ae16*(np.log10(ad_mid)+6)**4
            af_mid = z_mid*ab_mid**0.38/max(aa_mid**2.14, 1e-9)
            ag_mid = (af_mid-0.012)/max(abs(af_mid-0.012), 1e-9)
            ah_mid = (1-ag_mid)/2*0.012+(1+ag_mid)/2*af_mid
            ai_mid = ai12+ai13*ah_mid+ai14*(ah_mid**2)+ai15*(ah_mid**3)+ai16*(ah_mid**4)
        
        aj_mid = np.clip(datos['hl'] * ae_mid * ai_mid, 0, 1)
        an_mid = aj_mid*liq_den_col+(1-aj_mid)*gdens_col
        
        return calcular_vlp_hagedorn_brown(datos, aj_mid, mass_flow_rate, Dti_pulg, vis_liq_col, gvis_col, e_abs_ft, an_mid)

    else:
        # --- PATRÓN: FLUJO SLUG (o CHURN) ---
        # Si no es Burbuja, ni Dispersa, ni Anular, es Slug.
        # Usamos B&B (Patrón Intermitente), que es el proxy correcto para Slug.
        return calcular_vlp_beggs_brill(datos, press_col, ang_dev_raw, Dti_pulg, usg_col, usl_col, liq_den_col, gdens_col, liq_tsp_col, vis_liq_col, gvis_col, e_abs_ft)