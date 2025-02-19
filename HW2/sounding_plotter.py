""" 
    Script Created By Roz Roberts on 02/17/25
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from metpy.units import units as u
import metpy.calc as mc

current_directory = "."

sounding_files = [os.path.join(root, file) for root, dirs, files in os.walk(current_directory) for file in files if file.endswith(".csv")]
    
# print(sounding_files)

### Constants
R = 287 * u.joule / (u.kilogram * u.kelvin)
L_v = 2.501* (10**6) * u.joule / u.kilogram
C_p = 1004 * u.joule / (u.kilogram * u.kelvin)
eps = 0.622
dry_rate = (9.8/1000) * u.kelvin / u.meter
g = 9.8 * u.meters/(u.second ** 2) 



for sound in sounding_files:
    
    df = pd.read_csv(sound)
    
    # print(df)
    
    df.replace(-9999, np.nan, inplace=True)  # replaces empty values (-9999) with numpy nan values which can then be removed eaiser
    
    df.dropna(inplace=True)
    
    df.drop_duplicates(keep='first', inplace=True)
    
    press = df["pres[hPa]"].values * u.hPa  # Converts each datafram column into a numpy list for plotting, and sets the correct units for the value
    height = df["height[m]"].values * u.meters
    temp = df["temp[degC]"].values * u.degC
    dewT = df["dewpoint[degC]"].values * u.degC
    relh = df["rh[%]"].values
    mixr = df["mr[g/kg]"].values * u.grams / u.kilograms
    wddir = df["direction[deg]"].values * u.degrees
    wdspkt = df["speed[kt]"].values * u.knots
    theta = df["theta[K]"].values * u.kelvin
    thetaE = df["thetaE[K]"].values * u.kelvin
    thetaV = df["thetaV[K]"].values * u.kelvin
    
    
    ### HAND MADE PARCEL PROFILE
    lcl_height_aprox = (125 * u.meters) * (temp[0].m-dewT[0].m)  # This is our height (z) in meters
    
    # z = temp[0]/dry_rate (1 - (p/p0)^(R*dry_rate/g)) Hypsometric Equation for non-isothermal stmosphere 
    
    
    lcl_press_apox = press[0] * (1 - (lcl_height_aprox * dry_rate / temp[0]))**(g/(R.to('meter**2/(second**2 * kelvin)')*dry_rate))  # Hypsometric Equation in a non-isothermal atmosphere, rearanged to solve for pressure
    
    parcel_profile = [temp[0].m]
    
    lcl_lim = sum(press > lcl_press_apox)

    mixr_s_lst = []
    
    i = 0
    n=0
    for p in press:
        if n < lcl_lim:
            T_next = ((parcel_profile[i] * u.degC) * (press[i+1]/press[i])**(R/C_p)).to('degC').m
            parcel_profile.append(T_next)
            
            vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa
            mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))   # Saturation mixing ratio, unitless
            
            mixr_s_lst.append(mixr_s)
            
            n+=1
        else:
            if i+1 >= press.size:
                vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa
                mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))   # Saturation mixing ratio, unitless
                mixr_s_lst.append(mixr_s)
                continue
            delta_p =  press[i] - press[i+1]
            
            vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa
            mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))   # Saturation mixing ratio, unitless
            mixr_s_lst.append(mixr_s)
            
            T_top = R * ((parcel_profile[i] * u.degC).to("degK")) + L_v * mixr_s
            T_bottom = C_p + (((L_v**2) * mixr_s * eps)/(R*((parcel_profile[i] * u.degC).to("degK"))**2))
            
            T_delta = (delta_p / press[i]) * ((T_top/T_bottom))
            
            
            T_next = parcel_profile[i] - T_delta.m
            parcel_profile.append(T_next)
            
        i+=1
    
    parcel_profile = np.array(parcel_profile) * u.degC
    
    ### Handmade Parcel Path END
    
    u_wind, v_wind = mc.wind_components(wdspkt, wddir)  # Takes the wind speed in knots and wind direction in degrees and gets the UV components
    
    wind_barb_cuttoff = sum(press > 100 * u.hPa)
    
    u_wind = u_wind[:wind_barb_cuttoff]
    v_wind = v_wind[:wind_barb_cuttoff]
    
    
    parcel_prof_metpy = mc.parcel_profile(press, temp[0], dewT[0]).to('degC')
    lfc_press, lfc_temp = mc.lfc(press, temp, dewT)
    el_press, el_temp = mc.el(press, temp, dewT)
    lcl_press, lcl_temp = mc.lcl(press[0], temp[0], dewT[0])
    cape, cin = mc.cape_cin(press, temp, dewT, parcel_prof_metpy)
    
    
    ### Hand Generatred EL, LFC, CAPE, AND CIN VALUES
    
    path_diff_lcl = parcel_profile.to('degK') - temp.to('degK')
    
    l = 0 
    for j in path_diff_lcl[lcl_lim:]:
        if j < 0 and l==0:
            LFC_index = np.where(path_diff_lcl == j)[0][0]
            l+=1
    
    for j in np.flip(path_diff_lcl[lcl_lim:]):
        if j > 0 and l==1:
            EL_index = np.where(path_diff_lcl == j)[0][0] + 1
            l+=1
    
    EL_press_HG = press[EL_index]
    LFC_press_HG = press[LFC_index]
    
    EL_height = height[EL_index]
    LFC_height = height[LFC_index]
    
    
    temp_V_env = temp.to('degK') * (1 + 0.61*mixr.to('kg/kg'))
    temp_V_parce = parcel_profile.to('degK') * (1 + 0.61 * np.array(mixr_s_lst))
    
    
    
    CAPE = 0
    CIN = 0
    
    heights_LFC = height[height <= LFC_height]
    
    heights_EL = height[~(height < LFC_height)]
    heights_EL = heights_EL[heights_EL <= EL_height]
    
    
    for h in heights_LFC:
        h_ind = np.where(heights_LFC == h)[0][0] + 1
        if h_ind >= heights_LFC.size:
            continue
        else:
            delta_Z = heights_LFC[h_ind] - h
            
            CIN_partial = g * delta_Z * ((temp_V_parce[h_ind] - temp_V_env[h_ind])/temp_V_env[h_ind])
            
            CIN += CIN_partial
    
    CIN = CIN.to('J/kg')
    
    
    for h in heights_EL:
        h_ind = np.where(heights_EL == h)[0][0] + 1
        if h_ind >= heights_EL.size:
            continue
        else:
            delta_Z = heights_EL[h_ind] - h
            
            CAPE_partial = g * delta_Z * ((temp_V_parce[h_ind] - temp_V_env[h_ind])/temp_V_env[h_ind])
            
            CAPE += CAPE_partial
    
    CAPE = CAPE.to('J/kg')
    
    
    fig = plt.figure(figsize=(19,9))
    
    skew = SkewT(fig, rotation=45)
    
    skew.plot(press, temp, c='red', label='Temperature')  # Profile Lines
    skew.plot(press, dewT, c='green', label='Dewpoint')
    skew.plot(press, parcel_prof_metpy, c='black', linestyle='-.', label='Parcel Profile (MetPy Generated)')
    skew.plot(press, parcel_profile, c='fuchsia', linestyle='-.', label='Parcel Profile (Hand Generated)')
    
    
    
    skew.plot_barbs(press[:wind_barb_cuttoff], u_wind, v_wind)  # Wind Barbs

    
    t_sep = np.arange(-40, 60, 10)*u.degC  # Seperation Values
    m_sep = np.arange(1, 42, 2) * u.grams / u.kilograms
    p_sep = np.arange(400, 1000, 1) * u.hPa
    
    skew.plot_dry_adiabats(t0=t_sep, alpha=0.5, linestyle="--",color='orange')  # Normal SkewT lines
    skew.plot_mixing_lines(mixing_ratio=m_sep, pressure = p_sep, alpha=0.5, linestyle=':', color='teal')
    skew.plot_moist_adiabats(t0=t_sep, alpha=0.3, linestyle="-",color='g')
    
    
    # skew.ax.axhline(lfc_press, color='yellow', linestyle=':', label="Level of Free Convection (MetPy Generated)")
    skew.ax.axhline(LFC_press_HG,color='yellow', linestyle='--', label='Level of Free Convection (Hand Generated)')
    # skew.ax.axhline(el_press,color='black', linestyle='--', label='Equilibrium Level (MetPy Generated)')
    skew.ax.axhline(EL_press_HG,color='black', linestyle='--', label='Equilibrium Level (Hand Generated)')
    # skew.ax.axhline(lcl_press,color='purple', linestyle='--', label='Lifting Condensation Level (MetPy Generated)')
    skew.ax.axhline(lcl_press_apox,color='purple', linestyle='--', label='Lifting Condensation Level (Hand Generated)')
    
    
    text = (f'CAPE (MP): {cape:.2f} J/kg, CAPE (HG): {CAPE:.2f} J/kg\n'
                 f'CIN (MP): {cin:.2f} J/kg, CIN (HG): {CIN:.2f} J/kg\n'
                 f'LCL (MP): {lcl_press:.2f}; LCL (HG): {lcl_press_apox:.2f}\n'
                 f'LFC (MP): {lfc_press:.2f} hPa, LFC (HG): {LFC_press_HG:.2f} hPa\n'
                 f'EL (MP): {el_press:.2f} hPa, EL (HG): {EL_press_HG:.2f} hPa')
    skew.ax.text(0.02, 0.02, text, transform=skew.ax.transAxes, fontsize=8, 
                 bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
    
    
    plt.legend(loc='upper right', fontsize='small')
    plt.title(f'Skew-T Log-P Diagram: {sound[30:-4]}')
    skew.ax.set_xlabel('Temperature (°C)')
    skew.ax.set_ylabel('Pressure (hPa)')
    skew.ax.grid(True)





    # ax_hodo = fig.add_axes([0.68, 0.6, 0.28, 0.28])
    # h = Hodograph(ax_hodo, component_range=90)
    # h.add_grid()
    # h.plot_colormapped(u_wind, v_wind, press[:wind_barb_cuttoff], label='Hodograph')
    
    plt.savefig(f"{sound[:-4]}_results.pdf")
    plt.show()
    
    
    