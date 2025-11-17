# =====================================================================
# ESTE ES EL CÓDIGO COMPLETO PARA 'vlp_loop.py' (CON ANSARI RECONSTRUIDO)
# Bórralo todo y pega esto.
# =====================================================================

import numpy as np
from pvt import prop_oil_eos, dens_xfzg_tot, funcion_zg_exp
from flow import (calcular_vlp_gray, calcular_vlp_hagedorn_brown, 
                  calcular_vlp_beggs_brill, calcular_vlp_ansari) # <-- Ansari AÑADIDO
from utils import const

# Variable persistente para advertencias (como en MATLAB)
advertencia_pvt_mostrada = False

def calcular_vlp_para_q(ql_bpd, datos, Pb_yacimiento):
    """
    Traducción completa del bucle VLP de 'calcular_vlp_para_q'
    Ahora con el 'switch' de correlaciones VLP completo.
    """
    global advertencia_pvt_mostrada
    if ql_bpd <= 1e-6:
        advertencia_pvt_mostrada = False
        try:
            api_o_sg = datos['API'] if datos['Unidades'] == 0 else datos['API_o_SG']
            sg_liq = ( (141.5 / (131.5 + api_o_sg)) * (1-datos['Porcentaje_Agua']/100) + 
                       datos['Sw'] * (datos['Porcentaje_Agua']/100))
        except:
             sg_liq = 1.0
        dens_liq_estatica = sg_liq * 62.4
        gradiente_estatico = dens_liq_estatica / 144 # psi/ft
        return datos['Pcabezal'] + gradiente_estatico * datos['Profundidad']
    
    # Unpack de datos
    Pcab_psi = datos['Pcabezal']
    Tcab_F = datos['Tcabezal']
    Tfon_F = datos['Tfondo']
    Prof_ft = datos['Profundidad']
    Dti_pulg = datos['Dti']
    e_abs_ft = datos['Rugosidad'] / 12.0
    glr_scfstb = datos['RGL']
    wc_pct = datos['Porcentaje_Agua']
    api = datos['API']
    gg = datos['Sg']
    gw = datos['Sw']
    num_seg = datos['N_Segmentos']
    sg = 141.5 / (131.5 + api)
    YCO2, YH2s, YN2 = datos['YCO2'], datos['YH2s'], datos['YN2']

    # Constantes
    aera = np.pi / 4 * (Dti_pulg / 12)**2
    
    # Constantes H&B
    ac12=-2.69851; ac13=0.1584095; ac14=-0.55099756; ac15=0.54784917; ac16=-0.12194578
    ae12=-0.10306578; ae13=0.617774; ae14=-0.632946; ae15=0.29598; ae16=-0.0401
    ai12=0.91162574; ai13=-4.82175636; ai14=1232.25036621; ai15=-22253.57617; ai16=116174.28125
    
    # Constantes Viscosidad
    Za10=0.063069; Za11=-1.966847; Za12=21.0581; Za13=-27.0246
    
    # Cálculos iniciales de masa
    gas = glr_scfstb * ql_bpd
    qw = ql_bpd * wc_pct / 100
    aceite = ql_bpd * (1 - wc_pct / 100)
    
    RGA_m3 = 0.0
    if aceite > 1e-6:
        RGA_m3 = (gas / aceite) / const['m3m3_a_scfSTB_orig']
        
    liquid = 1.0
    if ql_bpd > 1e-6:
        liquid = ((ql_bpd - qw) * sg + qw * gw) / ql_bpd
        
    mass_gas_lpd = 0.0765 * gg * gas
    mass_liq_lpd = liquid * 62.37 * ql_bpd * (5.61458) # bbl/d a lb/d (aprox)
    
    mass_flow_rate = mass_gas_lpd + mass_liq_lpd
    Pb_eos_OPtion = datos['h_Ps_sanjari']

    # tpc/ppc
    if gg>=1.862 or api>=45: 
        tpc=164.3+357.7*gg-67.7*gg**2; ppc=744-125.4*gg+5.9*gg**2
    else: 
        tpc=120.1+429*gg-62.9*gg**2; ppc=671.1+14*gg-34.3*gg**2
        
    if YCO2>0 or YH2s>0 or YN2>0:
        YHCg=1-YCO2-YH2s-YN2
        Ppca=ppc*YHCg+1073*YCO2+1307*YH2s+492*YN2; Tpca=tpc*YHCg+547.7*YCO2+672.4*YH2s+227.1*YN2
        epsilon=107.6*((YCO2+YH2s)-(YCO2+YH2s)**2.2)+5.9*(YH2s**0.06-YH2s**0.68)
        tpc=Tpca-epsilon; ppc=Ppca*(Tpca-epsilon)/(Tpca+YH2s*(1-YH2s)*epsilon)
        
    # Inicializar bucle
    depth_col = np.linspace(0, Prof_ft, num_seg)
    press_col = np.zeros(num_seg)
    dp_pat_flw = 0.0 # dP del segmento ANTERIOR
    
    for i in range(num_seg):
        if i == 0:
            press_col_i = Pcab_psi
            temper_col_i = Tcab_F
        else:
            delta_L = depth_col[i] - depth_col[i-1]
            press_col_i = press_col[i-1] + dp_pat_flw * delta_L
            temper_col_i = Tcab_F + (Tfon_F - Tcab_F) / Prof_ft * depth_col[i]
        
        press_col[i] = press_col_i
        
        if press_col_i < 0: press_col_i = 14.7
        tem_Cel_col_i = (temper_col_i - 32) * 5 / 9

        try:
            (Dens_AguaStd, Dens_AguaPT, DensHC_L, DensHC_G, Zv_eos, PM_eos, 
             Ymol_VTP, Rs_VTP, Tb_eos, Tpr_eos, Ppr_eos, Bo_VPT, Pb_eos) = prop_oil_eos(
                 api, temper_col_i, press_col_i, Pb_eos_OPtion, gg, sg, Pb_yacimiento)
            
            (Ymol_GL, Masa_Gas, Masa_total, Dens_aceite, Dens_gas, Dens_liquido, 
             Xo_masa, Xa_masa, Xg_masa, Z_GL, Zgas, PM_GLib, PMgas_tot, Pr_GL, Dens_GLstd) = dens_xfzg_tot(
                 Zv_eos, YCO2, YH2s, YN2, Ymol_VTP, gg, api, PM_eos, temper_col_i, press_col_i, 
                 Rs_VTP, RGA_m3, wc_pct, DensHC_L, DensHC_G, Dens_AguaPT)
        except Exception as e:
             if not advertencia_pvt_mostrada: print(f"ADVERTENCIA (Py): Falla de PVT en P={press_col_i:.0f} psi. {e}")
             return np.nan
        
        liq_den_col_i = Dens_liquido / 16.01846
        gdens_col_i = Dens_gas / 16.01846

        if not np.isfinite(liq_den_col_i) or liq_den_col_i <= 0:
            if ql_bpd > 0 and not advertencia_pvt_mostrada: print("ADVERTENCIA (Py): PVT (Dens. Liq.) fuera de rango.")
            liq_den_col_i = 62.4
        if not np.isfinite(gdens_col_i) or gdens_col_i <= 0:
            if ql_bpd > 0 and not advertencia_pvt_mostrada: print("ADVERTENCIA (Py): PVT (Dens. Gas.) fuera de rango.")
            gdens_col_i = 0.1
        
        sup_ten_col = (1.11591-0.00305*tem_Cel_col_i)*(38.085-0.259*api)*np.exp(-0.00086306*press_col_i)*(0.056379+0.94362*np.exp(-0.0038491*Rs_VTP))
        wat_sup_col = (79.485-0.1922*tem_Cel_col_i)-0.1048*press_col_i**0.637
        liq_tsp_col_i = sup_ten_col*(1-wc_pct/100) + wat_sup_col*(wc_pct/100)
        if not np.isfinite(liq_tsp_col_i) or liq_tsp_col_i <= 0: liq_tsp_col_i = 30.0
            
        vis_oil_col = 1.0
        if datos['h_rb_gen_visc'] == 1:
            try:
                A_vis = 10.715 * (Rs_VTP + 100)**-0.515
                B_vis = 5.441 * (Rs_VTP + 150)**-0.338
                vis_oil_dead = A_vis * (temper_col_i**-B_vis)
                if vis_oil_dead < 0.1: vis_oil_dead = 0.1
                if press_col_i <= Pb_yacimiento:
                    vis_sat_col = vis_oil_dead
                else:
                    m_vazquez = 2.6 * (press_col_i**1.187) * np.exp(-11.513 - 8.98e-5 * press_col_i)
                    vis_sat_col = vis_oil_dead * (press_col_i / max(Pb_yacimiento, 1e-6))**m_vazquez
                vis_oil_col = vis_sat_col
            except:
                vis_oil_col = 1.0
        
        if not np.isfinite(vis_oil_col) or vis_oil_col <= 0: vis_oil_col = 1.0
        
        T_R = temper_col_i + 459.67
        Rho_g_pcm = gdens_col_i
        Rho_g_gcm3 = Rho_g_pcm * 0.0160185
        X_visc = 3.5 + (986 / T_R) + (0.01 * PM_GLib)
        Y_visc = 2.4 - (0.2 * X_visc)
        K_visc = (9.4 + 0.02 * PM_GLib) * (T_R**1.5) / (209 + 19 * PM_GLib + T_R)
        
        try:
            gvis_col_i = (K_visc * 1e-4) * np.exp(X_visc * (Rho_g_gcm3**Y_visc))
        except OverflowError:
            gvis_col_i = 0.02
            
        if not np.isfinite(gvis_col_i) or gvis_col_i <= 0: gvis_col_i = 0.01
            
        vis_agua = 1.0
        fw = wc_pct / 100.0
        vis_liq_col_i = (vis_oil_col * (1 - fw)) + (vis_agua * fw)
        if not np.isfinite(vis_liq_col_i) or vis_liq_col_i <= 0: vis_liq_col_i = 1.0
        
        k_mid_col = (temper_col_i + 459.67) / tpc
        j_mid_col = (press_col_i + 14.697) / ppc
        v_mid_col_i = funcion_zg_exp(j_mid_col, k_mid_col)
        if not np.isfinite(v_mid_col_i) or v_mid_col_i <= 0: v_mid_col_i = 0.9
            
        usg_col_i = 1/aera * gas * v_mid_col_i * (459.67 + temper_col_i) / (459.67 + 60) * (14.697 / max(press_col_i, 1e-6)) / 86400
        usl_col_i = mass_liq_lpd / max(aera * liq_den_col_i * 24 * 3600, 1e-9)

        try:
            y_mid = 1.938*usl_col_i*(liq_den_col_i/liq_tsp_col_i)**0.25
            z_mid = 1.938*usg_col_i*(liq_den_col_i/liq_tsp_col_i)**0.25
            aa_mid = 120.872*Dti_pulg/12*np.sqrt(liq_den_col_i/liq_tsp_col_i)
            ab_mid = 0.15726*vis_liq_col_i*(1/(liq_den_col_i*liq_tsp_col_i**3))**0.25
        except (ValueError, ZeroDivisionError):
             if not advertencia_pvt_mostrada: print(f"ADVERTENCIA (Py): Falla en Holdup (y_mid/z_mid). P={press_col_i:.0f} psi")
             return np.nan
        
        with np.errstate(all='ignore'):
            ac_mid = 10**(ac12+ac13*(np.log10(ab_mid)+3)+ac14*(np.log10(ab_mid)+3)**2+ac15*(np.log10(ab_mid)+3)**3+ac16*(np.log10(ab_mid)+3)**4)
            ad_mid = y_mid/max(z_mid**0.575, 1e-9)*((press_col_i/14.697)**0.1)*(ac_mid/max(aa_mid, 1e-9))
            ad_mid = np.clip(ad_mid, 0.00000161, 0.0047)
            ae_mid = ae12+ae13*(np.log10(ad_mid)+6)+ae14*(np.log10(ad_mid)+6)**2+ae15*(np.log10(ad_mid)+6)**3+ae16*(np.log10(ad_mid)+6)**4
            af_mid = z_mid*ab_mid**0.38/max(aa_mid**2.14, 1e-9)
            ag_mid = (af_mid-0.012)/max(abs(af_mid-0.012), 1e-9)
            ah_mid = (1-ag_mid)/2*0.012+(1+ag_mid)/2*af_mid
            ai_mid = ai12+ai13*ah_mid+ai14*(ah_mid**2)+ai15*(ah_mid**3)+ai16*(ah_mid**4)
        
        aj_mid_col_i = datos['hl'] * ae_mid * ai_mid
        aj_mid_col_i = np.clip(aj_mid_col_i, 0, 1)
        
        an_mid_i = aj_mid_col_i*liq_den_col_i+(1-aj_mid_col_i)*gdens_col_i
        if not np.isfinite(an_mid_i) or an_mid_i <= 0: an_mid_i = liq_den_col_i

        # =====================================================================
        # [INICIO DE CORRECCIÓN] Switch de VLP completo (CON ANSARI)
        # =====================================================================
        dp_pat_flw = 0.0
        if datos['h_rb_MFF_HB'] == 1:
            dp_pat_flw = calcular_vlp_hagedorn_brown(datos, aj_mid_col_i, mass_flow_rate, Dti_pulg, vis_liq_col_i, gvis_col_i, e_abs_ft, an_mid_i)
        elif datos['h_radiobutBeggsBrill'] == 1:
            dp_pat_flw = calcular_vlp_beggs_brill(datos, press_col_i, 0, Dti_pulg, usg_col_i, usl_col_i, liq_den_col_i, gdens_col_i, liq_tsp_col_i, vis_liq_col_i, gvis_col_i, e_abs_ft)
        elif datos['h_radiobutton_Ansari'] == 1:
            dp_pat_flw = calcular_vlp_ansari(datos, press_col_i, 0, Dti_pulg, usg_col_i, usl_col_i, liq_den_col_i, gdens_col_i, liq_tsp_col_i, vis_liq_col_i, gvis_col_i, e_abs_ft, mass_flow_rate)
        elif datos['h_CFM_GRAY'] == 1:
            dp_pat_flw = calcular_vlp_gray(datos, usg_col_i, usl_col_i, liq_den_col_i, gdens_col_i, liq_tsp_col_i, vis_liq_col_i, gvis_col_i, mass_flow_rate, e_abs_ft, Dti_pulg)
        else: # Default H&B
            dp_pat_flw = calcular_vlp_hagedorn_brown(datos, aj_mid_col_i, mass_flow_rate, Dti_pulg, vis_liq_col_i, gvis_col_i, e_abs_ft, an_mid_i)
        # =====================================================================
        # [FIN DE CORRECCIÓN]
        # =====================================================================

        if not np.isfinite(dp_pat_flw) or dp_pat_flw <= 0:
            dp_pat_flw = 0.001
            
        if ql_bpd > 0 and not advertencia_pvt_mostrada:
            advertencia_pvt_mostrada = True # Solo muestra advertencias una vez por Q

    return press_col[num_seg-1] # Devolver la Pwf final (último elemento)