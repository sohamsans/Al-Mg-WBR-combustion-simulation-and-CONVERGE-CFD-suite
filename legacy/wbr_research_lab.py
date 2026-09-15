import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
import numpy as np
import pandas as pd # For data export
from primary_combustion_model import simulate_1d_primary_combustor

def run_research_analysis():
    try:
        # 1. UI Parameters
        al_pct = float(al_entry.get()) / 100.0
        mg_pct = 1.0 - al_pct
        pc = float(pc_entry.get())
        eps = float(eps_entry.get())
        sub_frozen = bool(frozen_var.get())
        
        # 2. Manual Propellant Definitions (Bypasses all library errors)
        add_new_oxidizer('MyH2O', "ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-3788.5 t(k)=298.15")
        add_new_fuel('AlMg_Fuel', f"fuel AlMg_Fuel Al {al_pct} Mg {mg_pct} wt%=100.0 h,cal=0.0 t(k)=298.15")
        
        # 3. Initialize CEA
        isp_obj = CEA_Obj(oxName='MyH2O', fuelName='AlMg_Fuel')
        
        # 4. Sweep MR and Collect "Deep Data"
        mr_range = np.linspace(1.5, 10.0, 40)
        data = []

        for mr in mr_range:
            # 1. Performance and Temp 
            isp_vac = isp_obj.get_Isp(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
            temps = isp_obj.get_Temperatures(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
            
            # 2. Extract MW and Gamma
            mw_gam = isp_obj.get_Chamber_MolWt_gamma(Pc=pc, MR=mr, eps=eps)
            mw = mw_gam[0]
            gamma = mw_gam[1]

           # 3. Slag Tracking Logic (Refined)
            species_data = isp_obj.get_SpeciesMassFractions(pc, mr, eps)
            
            # Ensure we only take the Chamber dictionary (index 0)
            if isinstance(species_data, tuple):
                chamber_species = species_data[0]
            else:
                chamber_species = species_data
            
            slag_sum = 0.0
            # Identify specific Al-Mg slag products
            slag_targets = ['AL2O3(L)', 'AL2O3(s)', 'MGO(s)', 'MG(L)', 'AL(L)']
            
            for name, fraction in chamber_species.items():
                # Convert name to uppercase to match CEA standard
                if any(target in name.upper() for target in slag_targets):
                    # Extract float if it's trapped in a list
                    val = fraction[0] if isinstance(fraction, (list, np.ndarray)) else fraction
                    slag_sum += val

            # Final sanity check: 
            # If the sum is still > 1, it is definitely a percentage (0-100)
            if slag_sum > 1.0:
                slag_fraction = slag_sum / 100.0
            else:
                slag_fraction = slag_sum
                
            # Store results
            data.append({
                'MR': mr,
                'Isp': isp_vac,
                'T_Chamber': temps[0],
                'T_Exit': temps[2],
                'MW': mw,
                'Gamma': gamma,
                'Slag_Fraction': slag_fraction
            })

        df = pd.DataFrame(data)

        # 5. Visualizing the Fluid/Gas Properties
        for widget in plot_frame.winfo_children():
            widget.destroy()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        
        # Performance Graph
        ax1.plot(df['MR'], df['Isp'], 'b-', label='Isp (s)')
        ax1.set_ylabel('Vacuum Isp (s)', color='b')
        ax_temp = ax1.twinx()
        ax_temp.plot(df['MR'], df['T_Chamber'], 'r--', label='Chamber T')
        ax_temp.set_ylabel('Temp (K)', color='r')
        ax1.set_title("Performance & Thermal Profile")

        # Fluid Properties Graph (For CFD setup)
        ax2.plot(df['MR'], df['MW'], 'g-', label='Mol Weight')
        ax2.set_ylabel('Molecular Weight', color='g')
        ax_gamma = ax2.twinx()
        ax_gamma.plot(df['MR'], df['Gamma'], 'm:', label='Gamma')
        ax_gamma.set_ylabel('Ratio of Specific Heats (Gamma)', color='m')
        ax2.set_xlabel('Water-to-Fuel Ratio')
        ax2.set_title("Gas/Fluid Properties for CFD")

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Save for later
        global global_df
        global_df = df
        res_label.config(text=f"CEA Sweep Complete. {len(df)} points.", foreground="green")

    except Exception as e:
        messagebox.showerror("Error", str(e))

def run_expulsion_simulation():
    try:
        al_pct = float(al_entry.get())
        pc = float(pc_entry.get())
        dp_um = float(dp_entry.get())
        lc_m = float(lc_entry.get())
        
        df_traj, metrics = simulate_1d_primary_combustor(
            al_pct=al_pct,
            d_p0_um=dp_um,
            pc_psia=pc,
            primary_of=0.25,
            chamber_length_m=lc_m
        )
        
        for widget in plot_frame.winfo_children():
            widget.destroy()
            
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
        
        # Plot 1: Particle Diameter and Expulsion Efficiency vs Distance
        ax1.plot(df_traj['x_m'], df_traj['d_p_um'], 'b-', label='Droplet Diameter (µm)')
        ax1.set_ylabel('d_p (µm)', color='b')
        ax_eff = ax1.twinx()
        ax_eff.plot(df_traj['x_m'], df_traj['expulsion_eff_pct'], 'g--', label='Expulsion Eff (%)')
        ax_eff.axhline(90.0, color='r', linestyle=':', label='Target Threshold (90%)')
        ax_eff.set_ylabel('Expulsion Efficiency (%)', color='g')
        ax1.set_title('Primary Combustion: Droplet Size & Expulsion Efficiency')
        
        # Plot 2: Velocities vs Distance
        ax2.plot(df_traj['x_m'], df_traj['u_p_ms'], 'm-', label='Particle Vel u_p (m/s)')
        ax2.plot(df_traj['x_m'], df_traj['u_g_ms'], 'c--', label='Gas Vel u_g (m/s)')
        ax2.set_xlabel('Combustor Position x (m)')
        ax2.set_ylabel('Velocity (m/s)')
        ax2.set_title('Particle Acceleration & Velocity Tracking')
        ax2.legend()
        
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        
        eff_val = metrics['exit_expulsion_eff_pct']
        status_text = f"Expulsion Eff @ Exit: {eff_val:.1f}% | Slag: {metrics['total_slag_mass_pct']:.2f}%"
        color = "green" if eff_val >= 90.0 else "orange"
        res_label.config(text=status_text, foreground=color)
        
        global global_df
        global_df = df_traj
    except Exception as e:
        messagebox.showerror("Error", str(e))

def export_csv():
    if 'global_df' in globals():
        global_df.to_csv("WBR_Research_Data.csv", index=False)
        messagebox.showinfo("Success", "Data exported to WBR_Research_Data.csv")

# --- UI Setup ---
root = tk.Tk()
root.title("WBR Primary Combustion & Expulsion Physics Lab")
root.geometry("1100x750")

# Input Side
input_frame = ttk.LabelFrame(root, text="Configuration", padding=15)
input_frame.pack(side="left", fill="y", padx=10, pady=10)

ttk.Label(input_frame, text="Al Mass % (0-100):").pack(anchor="w")
al_entry = ttk.Entry(input_frame); al_entry.insert(0, "75"); al_entry.pack(fill="x", pady=4)

ttk.Label(input_frame, text="Chamber Pressure (psia):").pack(anchor="w")
pc_entry = ttk.Entry(input_frame); pc_entry.insert(0, "200"); pc_entry.pack(fill="x", pady=4)

ttk.Label(input_frame, text="Exit Area Ratio (eps):").pack(anchor="w")
eps_entry = ttk.Entry(input_frame); eps_entry.insert(0, "8"); eps_entry.pack(fill="x", pady=4)

ttk.Label(input_frame, text="Particle Size d_p0 (µm):").pack(anchor="w")
dp_entry = ttk.Entry(input_frame); dp_entry.insert(0, "50"); dp_entry.pack(fill="x", pady=4)

ttk.Label(input_frame, text="Primary Chamber Length L_c (m):").pack(anchor="w")
lc_entry = ttk.Entry(input_frame); lc_entry.insert(0, "0.20"); lc_entry.pack(fill="x", pady=4)

frozen_var = tk.IntVar()
ttk.Checkbutton(input_frame, text="Use Frozen Flow", variable=frozen_var).pack(pady=8)

ttk.Button(input_frame, text="Run CEA Water Sweep", command=run_research_analysis).pack(fill="x", pady=4)
ttk.Button(input_frame, text="Run 1D Primary Expulsion Model", command=run_expulsion_simulation).pack(fill="x", pady=4)
ttk.Button(input_frame, text="Export CSV for CFD", command=export_csv).pack(fill="x", pady=4)

res_label = ttk.Label(input_frame, text="Ready", foreground="gray")
res_label.pack(pady=15)

# Plot Side
plot_frame = ttk.Frame(root)
plot_frame.pack(side="right", fill="both", expand=True)

if __name__ == '__main__':
    root.mainloop()