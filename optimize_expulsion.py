"""
Optimization & Safe Operating Envelope Engine for Primary Metal Expulsion
==========================================================================
Executes parametric multi-variable sweeps over:
  - Al:Mg Mass Percentage (0% to 100% Al)
  - Initial Metal Particle Diameter d_p0 (10 to 100 um)
  - Primary Chamber Length L_c (0.05 to 0.50 m)
  - Primary Chamber Pressure P_c (150 to 500 psia)

Finds design boundaries satisfying the research objective:
  - Metal Expulsion Efficiency eta_expulsion >= 90%
  - Minimization of Slag Deposition
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from primary_combustion_model import simulate_1d_primary_combustor

def run_parametric_expulsion_sweep(
    al_pct_list=[25.0, 50.0, 75.0, 90.0],
    dp_um_list=[15.0, 30.0, 50.0, 75.0, 100.0],
    length_m_list=np.linspace(0.05, 0.50, 10),
    pc_psia=200.0,
    primary_of=0.25
):
    """
    Execute a grid search across propellant formulations and combustor geometries.
    """
    results = []
    print("Executing Parametric Expulsion & Slag Optimization Sweep...")
    
    total_runs = len(al_pct_list) * len(dp_um_list) * len(length_m_list)
    run_count = 0
    
    for al_pct in al_pct_list:
        for dp_um in dp_um_list:
            for L_c in length_m_list:
                run_count += 1
                try:
                    _, metrics = simulate_1d_primary_combustor(
                        al_pct=al_pct,
                        d_p0_um=dp_um,
                        pc_psia=pc_psia,
                        primary_of=primary_of,
                        chamber_length_m=L_c
                    )
                    
                    results.append({
                        'Al_pct': al_pct,
                        'Mg_pct': 100.0 - al_pct,
                        'd_p0_um': dp_um,
                        'L_c_m': L_c,
                        'Pc_psia': pc_psia,
                        'Primary_OF': primary_of,
                        'Expulsion_Eff_pct': metrics['exit_expulsion_eff_pct'],
                        'Slag_Pct': metrics['total_slag_mass_pct'],
                        'Residence_Time_ms': metrics['residence_time_ms'],
                        'Target_>=90_Met': metrics['target_met']
                    })
                except Exception as e:
                    print(f"Error at Al={al_pct}%, dp={dp_um}um, Lc={L_c}m: {e}")
                    
    df = pd.DataFrame(results)
    return df

def plot_safe_operating_envelope(df, output_img="expulsion_optimization_map.png"):
    """
    Generate diagnostic plots showing safe operating regions (Expulsion Efficiency >= 90%).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Expulsion Efficiency vs Chamber Length for various Particle Sizes (at Al 75%)
    df_sub = df[df['Al_pct'] == 75.0]
    for dp_um in sorted(df_sub['d_p0_um'].unique()):
        sub_dp = df_sub[df_sub['d_p0_um'] == dp_um]
        ax1.plot(sub_dp['L_c_m'], sub_dp['Expulsion_Eff_pct'], '-o', label=f'd_p0 = {dp_um} µm')
        
    ax1.axhline(90.0, color='r', linestyle='--', linewidth=2, label='Target Threshold (90%)')
    ax1.set_xlabel('Primary Chamber Length L_c (m)')
    ax1.set_ylabel('Metal Expulsion Efficiency (%)')
    ax1.set_title('Expulsion Efficiency vs Combustor Length (75% Al / 25% Mg)')
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.legend()
    
    # 2. Slag Accumulation vs Chamber Length
    for dp_um in sorted(df_sub['d_p0_um'].unique()):
        sub_dp = df_sub[df_sub['d_p0_um'] == dp_um]
        ax2.plot(sub_dp['L_c_m'], sub_dp['Slag_Pct'], '-s', label=f'd_p0 = {dp_um} µm')
        
    ax2.set_xlabel('Primary Chamber Length L_c (m)')
    ax2.set_ylabel('Predicted Slag Mass Fraction (%)')
    ax2.set_title('Slag Accumulation vs Combustor Length')
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    print(f"[SAVED] Optimization Contour Plot -> {output_img}")
    plt.close()

if __name__ == '__main__':
    df_results = run_parametric_expulsion_sweep()
    csv_out = "expulsion_optimization_results.csv"
    df_results.to_csv(csv_out, index=False)
    print(f"[SAVED] Optimization Sweep Results -> {csv_out}")
    
    plot_safe_operating_envelope(df_results)
    
    # Summary of design recommendations
    valid_designs = df_results[df_results['Target_>=90_Met'] == True]
    print("\n--- DESIGN RECOMMENDATIONS (Expulsion Efficiency >= 90%) ---")
    print(f"Total Valid Operating Points Found: {len(valid_designs)} out of {len(df_results)}")
    if len(valid_designs) > 0:
        print("\nSample Design Configurations:")
        print(valid_designs[['Al_pct', 'd_p0_um', 'L_c_m', 'Expulsion_Eff_pct', 'Slag_Pct']].head(10).to_string(index=False))
