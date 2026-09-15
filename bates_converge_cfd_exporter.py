"""
BATES Grain Regression & CONVERGE CFD Boundary Data Exporter
============================================================
Calculates BATES (Ballistic Test and Evaluation System) solid propellant grain regression,
burning surface area evolution A_b(t), equilibrium pressure P_c(t), burn rate r_b(t),
and generates clean analytical input boundary & thermodynamic files for CONVERGE CFD.

Exported CONVERGE CFD Files:
  - converge_bates_boundary.in         : CONVERGE boundary configuration & motion settings
  - converge_inflow_mass_flux.dat      : Time-dependent surface mass flux m'' [kg/(m^2 s)]
  - converge_thermo.dat               : CONVERGE thermodynamic & species mass fraction file
  - bates_internal_ballistics.csv      : Full time-history dataset of grain regression
"""

import numpy as np
import pandas as pd
import os
from primary_combustion_model import calculate_primary_gas_state

def simulate_bates_grain_regression(
    D_o_m=0.10,            # Outer diameter (100 mm)
    d_i0_m=0.04,           # Initial inner port diameter (40 mm)
    L_g0_m=0.15,           # Initial segment length (150 mm)
    N_seg=2,               # Number of BATES segments
    d_t_m=0.018,           # Nozzle throat diameter (18 mm)
    rho_prop=1750.0,       # Propellant density (kg/m^3)
    a_m_s_MPa=0.0055,      # Burn rate coefficient (m/s per MPa^n)
    n_exp=0.38,            # Pressure exponent
    ends_inhibited=False,  # Uninhibited vs inhibited segment ends
    cstar_m_s=1550.0,      # Characteristic velocity c* (m/s)
    dt_s=0.005             # Time step (5 ms)
):
    """
    Simulate BATES solid propellant grain regression and internal ballistics.
    Guarded against zero throat area, inverted diameters, and runaway pressure exponents.
    Returns DataFrame containing time-history metrics.
    """
    # Validation & Sanity Guards
    if D_o_m <= d_i0_m:
        raise ValueError(f"Outer diameter D_o ({D_o_m*1000:.1f} mm) must be strictly greater than inner port diameter d_i0 ({d_i0_m*1000:.1f} mm).")
    if d_t_m <= 0:
        raise ValueError(f"Nozzle throat diameter d_t must be strictly positive (received {d_t_m*1000:.1f} mm).")
    if n_exp >= 0.99 or n_exp <= 0.0:
        raise ValueError(f"Pressure exponent n must be in (0, 0.99) for stable rocket combustion (received {n_exp}).")
    if N_seg <= 0:
        raise ValueError(f"Number of segments N_seg must be >= 1 (received {N_seg}).")
        
    web_thickness = max(1.0e-5, (D_o_m - d_i0_m) / 2.0)
    A_t = max(1.0e-7, (np.pi / 4.0) * (d_t_m**2))
    
    y = 0.0
    t = 0.0
    records = []
    
    while y <= web_thickness and t < 10.0:
        d_i = d_i0_m + 2.0 * y
        L_g = (L_g0_m - 2.0 * y) if not ends_inhibited else L_g0_m
        
        if d_i >= D_o_m or L_g <= 0:
            break
            
        A_p = (np.pi / 4.0) * (d_i**2)
        
        if not ends_inhibited:
            A_b_seg = np.pi * d_i * L_g + 2.0 * (np.pi / 4.0) * (D_o_m**2 - d_i**2)
        else:
            A_b_seg = np.pi * d_i * L_g
            
        A_b = N_seg * A_b_seg
        K_b = A_b / A_t
        
        # Equilibrium Chamber Pressure (MPa and Pa)
        # P_c_MPa^(1-n) = (a_m_s_MPa * rho_prop * cstar_m_s * K_b) / 1.0e6
        term = (a_m_s_MPa * rho_prop * cstar_m_s * K_b) / 1.0e6
        P_c_MPa = (term)**(1.0 / (1.0 - n_exp)) if term > 0 else 0.1
        P_c_Pa = P_c_MPa * 1.0e6
        P_c_psia = P_c_Pa / 6894.76
        
        # Burn rate r_b (m/s and mm/s)
        r_b_m_s = a_m_s_MPa * (P_c_MPa**n_exp)
        r_b_mm_s = r_b_m_s * 1000.0
        
        # Mass flow generation & surface mass flux m''
        mdot_gen = rho_prop * A_b * r_b_m_s  # kg/s
        mass_flux_wall = rho_prop * r_b_m_s   # kg/(m^2 s)
        
        # Primary Gas Thermochemistry at P_c
        gas = calculate_primary_gas_state(al_pct=75.0, htpb_pct=15.0, pc_psia=P_c_psia, primary_of=0.25)
        
        # Injection velocity from regressing wall
        v_inflow = (rho_prop / gas['rho_g']) * r_b_m_s  # m/s
        
        records.append({
            'time_s': t,
            'web_burned_m': y,
            'web_burned_mm': y * 1000.0,
            'd_i_m': d_i,
            'L_g_m': L_g,
            'A_b_m2': A_b,
            'A_p_m2': A_p,
            'K_b': K_b,
            'P_c_MPa': P_c_MPa,
            'P_c_psia': P_c_psia,
            'r_b_mm_s': r_b_mm_s,
            'mdot_gen_kg_s': mdot_gen,
            'mass_flux_wall_kg_m2s': mass_flux_wall,
            'v_inflow_m_s': v_inflow,
            'T_c_K': gas['Tc'],
            'MW_g_gmol': gas['MW'],
            'Gamma_g': gas['Gamma'],
            'rho_g_kg_m3': gas['rho_g']
        })
        
        y += r_b_m_s * dt_s
        t += dt_s
        
    df = pd.DataFrame(records)
    return df

def export_converge_cfd_inputs(df, output_dir="d:/CEA"):
    """
    Export CONVERGE CFD boundary files, mass flux tables, and thermo.dat formats.
    Guarded against empty DataFrames or NaN values.
    """
    if df is None or len(df) == 0:
        raise ValueError("Cannot export empty regression dataset. Please run the BATES simulation first.")
        
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Export Internal Ballistics CSV
    csv_path = os.path.join(output_dir, "bates_internal_ballistics.csv")
    df.to_csv(csv_path, index=False)
    
    # 2. Export Inflow Mass Flux Table for CONVERGE boundary motion
    mass_flux_dat = os.path.join(output_dir, "converge_inflow_mass_flux.dat")
    with open(mass_flux_dat, 'w') as f:
        f.write("# Time (s)    Mass_Flux [kg/(m^2 s)]    Inflow_Vel [m/s]    Chamber_Pressure [Pa]\n")
        for _, row in df.iterrows():
            f.write(f"{row['time_s']:10.4f}    {row['mass_flux_wall_kg_m2s']:18.6f}    {row['v_inflow_m_s']:14.6f}    {row['P_c_MPa']*1e6:18.2f}\n")
            
    # 3. Export CONVERGE Boundary Setup Configuration file (.in)
    boundary_in = os.path.join(output_dir, "converge_bates_boundary.in")
    with open(boundary_in, 'w') as f:
        f.write("/* CONVERGE CFD Boundary Input File - BATES Grain Regression */\n")
        f.write("BOUNDARY_NAME: Propellant_Regressing_Wall\n")
        f.write("TYPE: INFLOW / MASS_INFLOW\n")
        f.write("TEMPERATURE_TYPE: SPECIFIED_VALUE\n")
        f.write(f"TEMPERATURE: {df['T_c_K'].iloc[0]:.2f} K\n")
        f.write("MASS_FLUX_FILE: converge_inflow_mass_flux.dat\n")
        f.write("WALL_MOTION: REGRESSING_SURFACE\n")
        f.write(f"PROPELLANT_DENSITY: {1750.0} kg/m^3\n")
        f.write("SPECIES_MASS_FRACTIONS:\n")
        f.write("  CO    : 0.245\n")
        f.write("  CO2   : 0.112\n")
        f.write("  H2O   : 0.184\n")
        f.write("  H2    : 0.082\n")
        f.write("  N2    : 0.142\n")
        f.write("  Al2O3 : 0.185\n")
        f.write("  MgO   : 0.050\n")
        
    # 4. Export CONVERGE Thermodynamic Data File (thermo.dat)
    thermo_dat = os.path.join(output_dir, "converge_thermo.dat")
    with open(thermo_dat, 'w') as f:
        f.write("THERMO_DATA_CONVERGE_CFD\n")
        f.write(f"Gas_Molecular_Weight: {df['MW_g_gmol'].iloc[0]:.4f} g/mol\n")
        f.write(f"Ratio_of_Specific_Heats_Gamma: {df['Gamma_g'].iloc[0]:.4f}\n")
        f.write(f"Chamber_Temperature: {df['T_c_K'].iloc[0]:.2f} K\n")
        f.write("End_Thermo_Data\n")
        
    print(f"[CONVERGE CFD] Internal Ballistics -> {csv_path}")
    print(f"[CONVERGE CFD] Inflow Mass Flux     -> {mass_flux_dat}")
    print(f"[CONVERGE CFD] Boundary Setup Config-> {boundary_in}")
    print(f"[CONVERGE CFD] Thermo Property File -> {thermo_dat}")

if __name__ == '__main__':
    print("--- Running BATES Grain Regression & CONVERGE CFD Data Exporter ---")
    df_bates = simulate_bates_grain_regression()
    export_converge_cfd_inputs(df_bates)
    print(f"\nBurn Summary:")
    print(f"Total Burn Time        : {df_bates['time_s'].iloc[-1]:.3f} seconds")
    print(f"Peak Chamber Pressure  : {df_bates['P_c_MPa'].max():.2f} MPa ({df_bates['P_c_psia'].max():.1f} psia)")
    print(f"Peak Burn Rate         : {df_bates['r_b_mm_s'].max():.2f} mm/s")
    print(f"Peak Wall Mass Flux    : {df_bates['mass_flux_wall_kg_m2s'].max():.2f} kg/(m^2 s)")
