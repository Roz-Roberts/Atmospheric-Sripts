import numpy as np
import matplotlib.pyplot as plt
from metpy.plots import SkewT
from metpy.calc import parcel_profile, lcl, cape_cin, lfc, el
from metpy.units import units as u
import pandas as pd
from siphon.simplewebservice.wyoming import WyomingUpperAir
from datetime import datetime
import os


def fetch_wyoming_data(date:datetime, station:str, file_path: str):
    # Ensure the directory exists
    os.makedirs(file_path, exist_ok=True)

    # Construct the full path for the CSV file
    csv_filename = f"sounding_{station}_{date.strftime('%Y%m%d_%H')}.csv"
    csv_filepath = os.path.join(file_path, csv_filename)



    # Fetch data
    try:
        print(f"Attempting to download sounding data for {station} on {date}...")
        df = WyomingUpperAir.request_data(date, station)
        df.to_csv(csv_filepath, index=False)
        print(f"Successfully downloaded and saved as {csv_filename}")

    except Exception as e:
        print(f"Failed to download data. Reason: {e}")
        # If download fails, check if a local CSV file exists
        if os.path.exists(csv_filepath):
            try:
                print(f"Reading local CSV file: {csv_filename}")
                df = pd.read_csv(csv_filepath)
            except Exception as e:
                print(f"Error reading local file. Reason: {e}")
                df = None
            
        else:
            print("No local data available. Exiting script.")
            df = None
        
    if df is not None:
        print("Returning Pandas Array")
        return df
    else:
        exit()

def plotter(df, date, station):
    # Convert data to appropriate units
    pressure = df['pressure'].values * u.hPa
    temperature = df['temperature'].values * u.degC
    dewpoint = df['dewpoint'].values * u.degC
    wind_speed = df['speed'].values * u.knots
    wind_direction = df['direction'].values * u.degrees

    wind_u = df['u_wind'].values * u.knots
    wind_v = df['v_wind'].values * u.knots

    # Masking and Skip values for Wind Barb Plotting
    skip_value = 5
    mask_100hPa = pressure >= 100 * u.hPa

    # Surface Parcel Path Plotting
    parcel_temps = parcel_profile(pressure, temperature[0], dewpoint[0])



    # Plot using Skew-T
    fig = plt.figure(figsize=(8, 10))
    skew = SkewT(fig, rotation=45)
    skew.plot(pressure, temperature, 'r', linewidth=2, label="Temperature")
    skew.plot(pressure, dewpoint, 'g', linewidth=2, label="Dewpoint")
    skew.plot_barbs(pressure[mask_100hPa][::skip_value], wind_u[mask_100hPa][::skip_value], wind_v[mask_100hPa][::skip_value])
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()


    skew.plot(pressure, parcel_temps, 'b', linestyle='dashed', linewidth=2, label='Parcel Path')

    lcl_press, lcl_temp = lcl(pressure[0], temperature[0], dewpoint[0])
    cape, cin = cape_cin(pressure, temperature, dewpoint, parcel_temps)
    lfc_pressure, lfc_temp = lfc(pressure, temperature, dewpoint)
    el_pressure, el_temp = el(pressure, temperature, dewpoint)


    skew.ax.plot(lcl_temp, lcl_press, 'ko', markersize=8, label="LCL")   # LCL in Black
    skew.ax.plot(lfc_temp, lfc_pressure, 'mo', markersize=8, label="LFC")   # LFC in Magenta
    skew.ax.plot(el_temp, el_pressure, 'co', markersize=8, label="EL")      # EL in Cyan


    skew.ax.axhline(lcl_press.m, color='black', linestyle=':', linewidth=2)
    x_cord = -40
    skew.ax.text(x_cord, lcl_press.m, f"{lcl_press:.1f}", ha = "left", va="center", fontsize=10, color= 'black', bbox=dict(facecolor='white', alpha=1))
    
    skew.ax.axhline(lfc_pressure.m, color='magenta', linestyle=':', linewidth=2)
    x_cord = -50
    skew.ax.text(x_cord, lfc_pressure.m, f"{lfc_pressure:.1f}", ha = "left", va="center", fontsize=10, color= 'magenta', bbox=dict(facecolor='white', alpha=1))
    
    skew.ax.axhline(el_pressure.m, color='cyan', linestyle=':', linewidth=2)
    x_cord = -90
    skew.ax.text(x_cord, el_pressure.m, f"{el_pressure:.1f}", ha = "left", va="center", fontsize=10, color= 'cyan', bbox=dict(facecolor='white', alpha=1))
    
    x_cord = -40
    skew.ax.text(x_cord, 200, f"CAPE and CIN {cape.m:.1f} & {cin.m:.1f} J/Kg", ha = "left", va="center", fontsize=10, color= 'black', bbox=dict(facecolor='white', alpha=1))
    

    plt.legend()
    plt.title(f"Skew-T Log-P Diagram for {station} on {date}")
    plt.show()



def main():
    # Define date and station
    date = datetime(2021, 12, 11, 0)  # YYYY, M, D, H (UTC time) - Converted Date 12/10/21, 5:00:00 PM (JUST BEFORE THE OUTBREAK)
    station = "LZK"  # Station code
    save_directory = "C:/Users/roswe/Desktop/WAF2/Atmospheric-Sripts/HW3"


    df = fetch_wyoming_data(date, station, save_directory)
    plotter(df, date, station)

main()
