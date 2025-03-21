""" 
    Script Created By Roz Roberts on 02/17/25
    
    Last Updated: Roz Roberts - 02/19/25
    
    WARNING: THIS SCRIPT SCANS CURRENT WORKING DIRECTORY FOR .CSV FILES, INCORRECTLY FORMATTED CSV's WILL CRASH THIS SCRIPT
    
    Packages Used: OS, numpy, pandas, matplotlib, metpy
    
    Formulas Used:
    Lifting Condensation Level Height Approximation; Z = 125 * (T-T_d)
    Non-isothermal Hypsometric Equation* (Rearranged to solve for Pressure); P = P0 * (1 - (Z * Gamma / T))^(g / (R * Gamma))
    Dry Parcel "Next Temperature" (DPNT); T2 = T1 * (P2/P1)^(R/C_p)
    Saturation Vapor Pressure Equation; V_s = 6.112 * exp(17.67 * T0 / (T0 + 243.5))
    Saturation Mixing Ratio Equation; W_s = epsilon * (V_s / (P0 - V_s))
    Saturated Approximate "Next Temperature" (SPANT); D_T = (D_P/P0) * (R * T0 + L_v * W_s) / (C_p + (((L_v^2) * W_s * epsilon) / (R * (T0)^2)))
    Virtual Temperature Equation; T_v = T * (1 + 0.61 * W_s)
    CIN/CAPE "Integration" Equation; C = g * D_Z * ((T_v_Parcel - T_v_Env)/T_v_Env)
    
    * Non-isothermal Hypsometric Equation Collected from Wallace and Hobbs 1977 pages 60-61
    
    Outputs: PDF Figures of the results, as well as a CSV of each of the values calculated. Both CSV and PDF end with "_results" for easy identification.
"""


import os  # All imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from metpy.units import units as u
import metpy.calc as mc

current_directory = "."  # The Current Working Directory of this script

# Scans all files in working directory then collects the CSV files that are used in the following For loop, excludes results from sounding_plotter.py
sounding_files = [os.path.join(root, file) for root, dirs, files in os.walk(current_directory) for file in files if file.endswith(".csv") and not file.endswith('_results.csv')]



### List of all Constants used in this script
R = 287 * u.joule / (u.kilogram * u.kelvin)  # Specific gas constant for Dry Air
L_v = 2.501* (10**6) * u.joule / u.kilogram  # Latent heat of vaporization of water
C_p = 1004 * u.joule / (u.kilogram * u.kelvin)  # Heat Capacity (ISOBARIC) of water
eps = 0.622  # Ratio of gas constant of water to air, unitless version (kg/kg)
dry_rate = (9.8/1000) * u.kelvin / u.meter  # Dry Adiabatic Lapse Rate
g = 9.8 * u.meters/(u.second ** 2)  # Normal gravity on Earth



def interpolate_intersection(i, pressure_lim, temp_diff, temp_env_lim, height_lim):  # The function for finding the interpolated pressure, temperature, and height
    p1, p2 = pressure_lim[i], pressure_lim[i+1]  # Initializations for the current values to interpolate between
    t1, t2 = temp_diff[i], temp_diff[i+1]
    h1, h2 = height_lim[i], height_lim[i+1]
    
    # Linear interpolation for pressure level, temperature, and altitude
    p_intersect = p1 - t1 * (p2 - p1) / (t2 - t1)
    t_intersect = temp_env_lim[i] + (temp_env_lim[i+1] - temp_env_lim[i]) * (p_intersect - p1) / (p2 - p1)
    h_intersect = h1 + ((p_intersect - p1)/(p2 - p1)) * (h2 - h1)
    
    return p_intersect, t_intersect, h_intersect  # Returns desired values



for sound in sounding_files:  # Loop for to process each .CSV file individually
    
    df = pd.read_csv(sound)  # Gathers the data from the .CSV's and imports it into a Pandas dataframe
    
    
    df.replace(-9999, np.nan, inplace=True)  # Replaces empty values (-9999 in the CSV) with numpy nan values which can then be removed easier
    
    df.dropna(inplace=True)  # Removes nan values
    
    df.drop_duplicates(subset=df.columns[0], keep='first', inplace=True)  # Removes duplicate values, while retaining the first instance of a value - this fixes graphing issues later on
    
    # Below converts each dataframe column into a smaller dataframe for plotting and manipulation, and sets the correct units for the value (units from CSV headers)
    press = df["pres[hPa]"].values * u.hPa  # Pressure in Hectopascals
    height = df["height[m]"].values * u.meters  # Altitude in Meters
    temp = df["temp[degC]"].values * u.degC  # Environmental Temperature in Celsius
    dewT = df["dewpoint[degC]"].values * u.degC  # Environmental Dewpoint in Celsius
    relh = df["rh[%]"].values  # Relative Humidity in Percentage
    mixr = df["mr[g/kg]"].values * u.grams / u.kilograms  # Mixing Ratio in Grams per Kilogram (NOT UNITLESS VERSION)
    wddir = df["direction[deg]"].values * u.degrees  # Wind Direction in Degrees
    wdspkt = df["speed[kt]"].values * u.knots  # Wind Speed in Knots
    theta = df["theta[K]"].values * u.kelvin  # Potential Temperature in Kelvin
    thetaE = df["thetaE[K]"].values * u.kelvin  # Equivalent Potential Temperature in Kelvin
    thetaV = df["thetaV[K]"].values * u.kelvin  # Virtual Potential Temperature in Kelvin
    
    
    ### HAND MADE PARCEL PROFILE
    lcl_height_aprox = (125 * u.meters) * (temp[0].m-dewT[0].m)  # The Lifting Condensation Level (LCL) in meters
    lcl_press_aprox = press[0] * (1 - (lcl_height_aprox * dry_rate / temp[0]))**(g/(R.to('meter**2/(second**2 * kelvin)')*dry_rate))  # Hypsometric Equation in a non-isothermal atmosphere, rearranged to solve for pressure of the LCL
    
    parcel_profile = [temp[0].m]  # Initializes our hand made Parcel Profile, and sets the first Temperature to be equal to the initial Environmental Temperature
    
    lcl_lim = sum(press > lcl_press_aprox)  # Finds all pressure values less than the Lifting Condensation Level, as we need to treat each area differently 

    mixr_s_lst = []  # Initializes the Saturated Mixing Ratio profile 
    
    i = 0  # Both i (pressure) and n (LCL) are indexes
    n=0
    for p in press:  # This loops goes through each pressure level and depending on if we are below or above the LCL we calculate the next Parcel Profile Temperature, and Saturation Mixing Ratio
        if n < lcl_lim:  # Using the Dry Adiabatic lapse rate we calculate the next temperature in the profile, BELLOW THE LCL
            T_next = ((parcel_profile[i] * u.degC) * (press[i+1]/press[i])**(R/C_p)).to('degC').m  # The Next Temperature in the profile
            parcel_profile.append(T_next)  # Adds the next temperature to the Parcel Profile
            
            vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa  # Calculates the Saturation Vapor Pressure (SVP) of the parcel
            mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))   # Calculates the Saturation Mixing Ratio (unitless form) (SMR) from the SVP
            
            mixr_s_lst.append(mixr_s)  # Adds SMR to SMR Profile
            
            n+=1  # Increments n index until we are greater than the LCL limit value specified earlier 
        
        else:  # Now we treat every value above the LCL
            if i+1 >= press.size:  # Handles the last value in the profile - special case
                vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa
                mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))
                mixr_s_lst.append(mixr_s)
                continue  # Ends Loop as we just calculated the last mixing ratio values
            
            delta_p =  press[i] - press[i+1]  # Calculates the difference in pressure form the next pressure level, and the current pressure level
            
            vap_press_s = 6.112 * np.exp(17.67*parcel_profile[i]/(parcel_profile[i]+243.5)) * u.hPa
            mixr_s = eps * (vap_press_s / (press[i]-vap_press_s))
            mixr_s_lst.append(mixr_s)
            
            T_top = R * ((parcel_profile[i] * u.degC).to("degK")) + L_v * mixr_s  # Calculates the numerator of the Saturated Parcel Approximate "Next Temperature" (SPANT)
            T_bottom = C_p + (((L_v**2) * mixr_s * eps)/(R*((parcel_profile[i] * u.degC).to("degK"))**2))  # Calculates the denominator of the SPANT
            
            T_delta = (delta_p / press[i]) * ((T_top/T_bottom))  # Combines the difference in temperatures and the numerator and denominator to get the full SPANT expression
            
            
            T_next = parcel_profile[i] - T_delta.m  # Calculates the next temperature from the temperature changed calculated by the SPANT
            parcel_profile.append(T_next)  # Adds Next Temperature to Parcel Profile
            
        i+=1  # Increases i index for pressure index tracking
    
    parcel_profile = np.array(parcel_profile) * u.degC  # Changes parcel profile into Celsius values, and finishes the parcel profile
    
    
    
    u_wind, v_wind = mc.wind_components(wdspkt, wddir)  # Takes the wind speed in knots and wind direction in degrees and gets the UV components
    
    wind_barb_cutoff = sum(press > 100 * u.hPa)  # Limits wind-barbs to 100 hPa minimum
    
    u_wind = u_wind[:wind_barb_cutoff]
    v_wind = v_wind[:wind_barb_cutoff]
    
    
    parcel_prof_metpy = mc.parcel_profile(press, temp[0], dewT[0]).to('degC')  # Metpy generated values: Parcel Path, LFC, EL, LCL, CAPE, CIN
    lfc_press, lfc_temp = mc.lfc(press, temp, dewT)
    el_press, el_temp = mc.el(press, temp, dewT)
    lcl_press, lcl_temp = mc.lcl(press[0], temp[0], dewT[0])
    cape, cin = mc.cape_cin(press, temp, dewT, parcel_prof_metpy)
    
    
    ### Hand Generated EL, LFC, CAPE, AND CIN VALUES
    
    temp_diff = temp.to('degK') - parcel_profile.to('degK')  # Calculates total temperature difference between the parcel and the environment
    
    temp_diff = temp_diff.m  # Gets only the magnitude, without the "delta_degrees_celsius"
    
    temp_diff_lim = temp_diff[lcl_lim:]  # Truncates data to what is needed for interpolation for EL and LFC
    
    
    sign_changes = np.where(np.diff(np.sign(temp_diff_lim)))[0]  # Finds every index where the sign of the difference changes, indicating an intersection
    
    if sign_changes.size == 0:  # If there are no sign changes then there is no EL or LFC, we set these values to Zero or Nan due to this, we will use the None indexes for future if-statements
        EL_index = LFC_index = None
        EL_press_HG = 0 * u.hPa
        EL_height_HG = np.nan * u.meters
        LFC_press_HG = 0 * u.hPa
        LFC_height_HG = np.nan * u.meters
        
    else:  # Otherwise we get the indexes of the EL and LCL; then interpolate to find the exact pressure level for the EL and LCL
        LFC_index = sign_changes[0] + lcl_lim
        EL_index = sign_changes[-1] + lcl_lim
        
        LFC_press_HG, LFC_temp_HG,LFC_height_HG = interpolate_intersection(LFC_index, press, temp_diff, temp, height)  # Calculating the exact intersection point
        EL_press_HG, EL_temp_HG, EL_height_HG = interpolate_intersection(EL_index, press, temp_diff, temp, height)
    
    
    
    temp_V_env = temp.to('degK') * (1 + 0.61*mixr.to('kg/kg'))  # Calculates the virtual temperature of the environment
    temp_V_parcel = parcel_profile.to('degK') * (1 + 0.61 * np.array(mixr_s_lst))  # Calculates the virtual temperature of the parcel
    
    
    CAPE = [0 * (u.meters**2 / u.sec**2)] # Initializes CAPE and CIN values (will be added to while we "integrate")
    CIN = [0 * (u.meters**2 / u.sec**2)]
    
    heights_LFC = height[height <= LFC_height_HG]  # Finds all height values less than the LFC
    
    heights_EL = height[~(height < LFC_height_HG)]  # Finds all height values more than the LFC
    heights_EL = heights_EL[heights_EL <= EL_height_HG]  # Cuts off the height values to less than the EL
    
    
    for h in heights_LFC:  # "Integration" to find CIN
        h_ind = np.where(heights_LFC == h)[0][0] + 1  # A check to determine if we are on the last value or not
        if h_ind >= heights_LFC.size:
            continue
        
        else:  # After Check - Integration!
            delta_Z = heights_LFC[h_ind] - h  # Calculates the change in height between the current height and the next
            
            CIN_partial = g * delta_Z * ((temp_V_parcel[h_ind] - temp_V_env[h_ind])/temp_V_env[h_ind])  # Calculates CIN
            
            CIN.append(CIN_partial)  # Adds the CIN to the total CIN
    
    CIN_val = sum(CIN).to('J/kg')  # Converts CIN to the correct Values, and finishes the CIN Calculation
    
    
    for h in heights_EL:  # Exact Same Setup as for CIN just different Bounds
        h_ind = np.where(heights_EL == h)[0][0] + 1
        if h_ind >= heights_EL.size:
            continue
        
        
        else:
            delta_Z = heights_EL[h_ind] - h
            
            CAPE_partial = g * delta_Z * ((temp_V_parcel[h_ind] - temp_V_env[h_ind])/temp_V_env[h_ind])
            
            CAPE.append(CAPE_partial)
    
    CAPE_val = sum(CAPE).to('J/kg')  # Finishes off the CAPE Calculation
    
    if EL_index == None:  # If there is no LFC/EL then there is no CAPE and CIN, we set these to zero to indicate this fact (like earlier!)
        CIN_val = 0.0 * u.joule / u.kilogram
        CAPE_val = 0.0 * u.joule / u.kilogram
    
    
    
    temp_gradient = np.gradient(temp)  # Calculates the temperature gradient of the environment
    
    inversion_points = np.where(temp_gradient >= 0)[0]  # Finds all locations where the gradient is zero or positive, indicating an inversion is occurring
    
    
    if len(inversion_points) > 0:  # If there is more than one inversion location, we need to select the first instance
        inversion_point = inversion_points[0]  # Selects first inversion
        inversion_pressure = press[inversion_point]  # Finds the inversion pressure
        lft_inversion_pressure = None
        
        if inversion_pressure == press[0]:  # If, however, the first "inversion" is at the surface we also want to indicate the first inversion aloft as well, so we check
            inversion_point = inversion_points[1]  # Finds next inversion location
            lft_inversion_pressure = press[inversion_point]  # Finds the next inversion location pressure
        
    else:
        inversion_pressure == np.nan * u.hPa  # If there is NO inversion we indicate this by setting the inversion pressure to NAN
    
    # 1-Kilometer Wind Shear Estimation
    
    km_wnd_indx = np.where(height == height[height > 1000 * u.meter][0])[0][0]  # Finds the first pressure level above 1 kilometer to calculate shear from
    u_shear, v_shear = u_wind[km_wnd_indx] - u_wind[0], v_wind[km_wnd_indx] - v_wind[0]  # Finds the U and V components of the shear
    
    shear_spdkt = np.sqrt(u_shear**2 + v_shear**2)  # Calculates Shear speed in knots
    shear_dir = np.degrees(np.arctan(v_shear/u_shear))  # Calculates Shear direction in degrees
    if shear_dir < 0:  # Makes sure that negative degree values are correctly positive
        shear_dir = 360 + shear_dir
    
    
    fig = plt.figure(figsize=(19,9))  # Sets plotting figure
    
    skew = SkewT(fig, rotation=45)  # Creates the Skew-T
    
    skew.plot(press, temp, c='red', label='Temperature')  # Plots all the Profiles
    skew.plot(press, dewT, c='green', label='Dewpoint')
    skew.plot(press, parcel_prof_metpy, c='black', linestyle='-', label='Parcel Profile (MetPy Generated)')
    skew.plot(press, parcel_profile, c='fuchsia', linestyle='-.', label='Parcel Profile (Hand Generated)')
    
    
    
    skew.plot_barbs(press[:wind_barb_cutoff], u_wind, v_wind)  # Plots the Wind Barbs

    
    t_sep = np.arange(-40, 60, 10)*u.degC  # Formatting Separation Values
    m_sep = np.arange(1, 42, 2) * u.grams / u.kilograms
    p_sep = np.arange(400, 1000, 1) * u.hPa
    
    skew.plot_dry_adiabats(t0=t_sep, alpha=0.5, linestyle="--",color='orange')  # Plots the Base SkewT lines based on separation values
    skew.plot_mixing_lines(mixing_ratio=m_sep, pressure = p_sep, alpha=0.5, linestyle=':', color='teal')
    skew.plot_moist_adiabats(t0=t_sep, alpha=0.3, linestyle="-",color='g')
    
    
    # skew.ax.axhline(lfc_press, color='yellow', linestyle=':', label="Level of Free Convection (MetPy Generated)")
    skew.ax.axhline(LFC_press_HG,color='yellow', linestyle='--', label='Level of Free Convection (Hand Generated)')  # LFC Line Plotted
    # skew.ax.axhline(el_press,color='black', linestyle='--', label='Equilibrium Level (MetPy Generated)')
    skew.ax.axhline(EL_press_HG,color='black', linestyle='--', label='Equilibrium Level (Hand Generated)')  # EL Line Plotted
    # skew.ax.axhline(lcl_press,color='purple', linestyle='--', label='Lifting Condensation Level (MetPy Generated)')
    skew.ax.axhline(lcl_press_aprox,color='purple', linestyle='--', label='Lifting Condensation Level (Hand Generated)')  # LCL Line Plotted
    
    if lft_inversion_pressure != None:
        skew.ax.axhline(inversion_pressure,color='Orange', linestyle='--', label='Lowest Inversion (Hand Generated)')  # Inversion Level Plotted
        skew.ax.axhline(lft_inversion_pressure,color='goldenrod', linestyle='-.', label='First Inversion Aloft (Hand Generated)')  # Inversion Level Plotted
    else:
        skew.ax.axhline(inversion_pressure,color='Orange', linestyle='--', label='Lowest Inversion (Hand Generated)')
    
    
    
    text = (f'CAPE (MP): {cape:.2f}, CAPE (HG): {CAPE_val:.2f}\n'  # Formats all result values we want to show on the final figure
                 f'CIN (MP): {cin:.2f}, CIN (HG): {CIN_val:.2f}\n'
                 f'LCL (MP): {lcl_press:.2f}; LCL (HG): {lcl_press_aprox:.2f}\n'
                 f'LFC (MP): {lfc_press:.2f}, LFC (HG): {LFC_press_HG:.2f}\n'
                 f'EL (MP): {el_press:.2f}, EL (HG): {EL_press_HG:.2f}\n'
                 f'LFC Alt. (HG): {LFC_height_HG:.2f}, EL Alt. (HG): {EL_height_HG:.2f}\n'
                 f'Surface Pressure {press[0]:.2f}\n'
                 f'Lowest Inversion Pressure {inversion_pressure:.2f}\n'
                 f'Wind Shear Speed [kts]: {shear_spdkt:.2f}, Direction [deg]: {shear_dir:.2f}\n')
    
    if lft_inversion_pressure != None:  # Text considerations for when we do and don't have two plotted inversion layers
        text = text + f"First Lifted Inversion Pressure {lft_inversion_pressure:.2f}\n" + f'MP means MetPy Generated and HG means Hand Generated\nThis script is in part Powered by MetPy'
    else:
        text = text + f'MP means MetPy Generated and HG means Hand Generated\nThis script is in part Powered by MetPy'
    
    skew.ax.text(0.02, 0.02, text, transform=skew.ax.transAxes, fontsize=8,  # Plots the formatted results onto our final figure
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
        
    
    plt.legend(loc='upper right', fontsize='small')  # Plots a Legend for all lines on the plot
    plt.title(f'Skew-T Log-P Diagram: {sound[-18:-4]}')  # Plots the title
    skew.ax.set_xlabel('Temperature (°C)')  # Plots the x-axis title
    skew.ax.set_ylabel('Pressure (hPa)')  # Plots the y-axis title
    skew.ax.grid(True)  # Plots the background grid


    # THE BELLOW COMMENTED SECTION IS EXTRA CODE TO PLOT A HODOGRAPH TO THE FINAL FIGURE AND COLOR MAP IT
    # ax_hodo = fig.add_axes([0.68, 0.6, 0.28, 0.28])
    # h = Hodograph(ax_hodo, component_range=90)
    # h.add_grid()
    # h.plot_colormapped(u_wind, v_wind, press[:wind_barb_cutoff], label='Hodograph')
    
    
    CIN_CAPE = CIN.copy()
    CIN_CAPE.extend(CAPE)
    
    
    output = {'pressure[hPa]' : press.m,  # The dictionary of all profiles produced by this code
              'height[m]' : height.m,
              'Temperature[degC]' : temp.m,
              'Dewpoint[degC]' : dewT.m,
              'Relative Humidity[%]' : relh,
              'Environmental Mixing Ratio[g/kg]' : mixr.m,
              'Wind Direction[deg]' : wddir.m,
              'Wind Speed[kt]' : wdspkt.m,
              'Theta[K]' : theta.m,
              'ThetaE[K]' : thetaE.m,
              'ThetaV[K]' : thetaV.m,
              'Parcel Temperature [degC]' : parcel_profile.m,
              'Saturation Mixing Ratio for the Parcel [unit-less]' : [mixs.m for mixs in mixr_s_lst],
              'U Wind Component' : u_wind.m,
              'V Wind Component' : v_wind.m,
              'Parcel-Environment Temperature Difference [degC]' : temp_diff,
              'Environmental Virtual Temperature [degC]' : temp_V_env.m,
              'Parcel Virtual Temperature [degC]' : temp_V_parcel.m,
              'CIN/CAPE (CIN/CAPE Changes past the LFC)[J/kg]' : CIN_CAPE,
              'CIN (To Show where the LFC is)[J/kg]' : CIN,
              'Temp Gradient [degC]' : temp_gradient.m
              }
    
    
    
    
    profiles = {key: pd.Series(value) for key,value in output.items()}  # Turns dictionary into a list of pandas series objects
    
    fdf = pd.concat(profiles, axis=1)  # Turns the series objects into one dataframe
    csv_filename = f"{sound[:-4]}_results.csv"  # Name of CSV
    
    fdf.to_csv(csv_filename, index=False, mode='w')  # Writes the different profiles out to a CSV
    
    plt.savefig(f"{sound[:-4]}_results.pdf")  # Automatically Saves the Final Figures as a PDF
    
    plt.show()  # Shows each figure once it is finished being created
    
    