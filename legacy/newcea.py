import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
import numpy as np
import pandas as pd

def run_research_analysis():
    try:
        # 1. UI Parameters
        al_pct = float(al_entry.get()) / 100.0
        mg_pct = 1.0 - al_pct
        pc = float(pc_entry.get())
        eps = float(eps_entry.get())
        mdot = float(mdot_entry.get())
        sub_frozen = bool(frozen_var.get())
        
       # 2. Propellant Blending (DRDO Baseline AP + HTPB)
        # Explicitly define AP Oxidizer
        add_new_oxidizer('AP_Oxidizer', "ox NH4CLO4(I) N 1.0 H 4.0 CL 1.0 O 4.0 wt%=100.0 h,cal=-70690.0 t(k)=298.15")
        
        # Calculate mass weights (15% HTPB binder, 85% Metals)
        htpb_wt = 15.0
        al_wt = al_pct * 85.0  
        mg_wt = mg_pct * 85.0  
        
        # Calculate Moles of each element per 100g of the fuel mixture
        # HTPB (C 7.07 H 10.12 O 0.20) has a Molecular Weight of ~98.32 g/mol
        moles_htpb = htpb_wt / 98.32
        moles_C = 7.07 * moles_htpb
        moles_H = 10.12 * moles_htpb
        moles_O = 0.20 * moles_htpb
        
        # Moles of Metals
        moles_Al = al_wt / 26.98
        moles_Mg = mg_wt / 24.31
        
        # Calculate Bulk Enthalpy (HTPB is -12.5 cal/g, Metals are 0)
        bulk_h_cal = (htpb_wt * -12.5) / 100.0
        
        # Build ONE bulletproof chemical string balancing all atoms
        fuel_string = (f"fuel DRDO_Fuel C {moles_C:.4f} H {moles_H:.4f} O {moles_O:.4f} "
                       f"AL {moles_Al:.4f} MG {moles_Mg:.4f} "
                       f"wt%=100.0 h,cal={bulk_h_cal:.2f} t(k)=298.15")
        
        add_new_fuel('DRDO_Fuel', fuel_string)
        
        # 3. Initialize CEA
        isp_obj = CEA_Obj(oxName='AP_Oxidizer', fuelName='DRDO_Fuel')
        
        # 4. Sweep MR and Collect "Deep Data"
        mr_range = np.linspace(0.1, 2.5, 40)
        data = []

        for mr in mr_range:
            # 1. Performance and Temp 
            isp_vac = isp_obj.get_Isp(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
            
            # Calculate Thrust in Newtons
            g0 = 9.80665
            thrust_N = mdot * isp_vac * g0
            thrust_kN = thrust_N / 1000.0 
            
            temps = isp_obj.get_Temperatures(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
            
            # 2. Extract MW and Gamma
            mw_gam = isp_obj.get_Chamber_MolWt_gamma(Pc=pc, MR=mr, eps=eps)
            mw = mw_gam[0]
            gamma = mw_gam[1]

            # 3. Slag Tracking Logic (Bulletproofed)
            molWtD, massFracD = isp_obj.get_SpeciesMassFractions(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
            
            slag_sum = 0.0
            unburned_metal_sum = 0.0
            
            for name, fraction_list in massFracD.items():
                clean_name = name.strip().upper().replace('*', '')
                
                try:
                    val = float(fraction_list[0])
                except (ValueError, TypeError, IndexError):
                    val = 0.0
                    
                # Identify True Slag (Oxides) safely
                if 'AL2O3' in clean_name or 'MGO' in clean_name:
                    slag_sum += val
                    
                # Identify Unburned Metal (Expulsion Carryover) safely
                elif clean_name in ['AL(L)', 'AL(CR)', 'AL(S)', 'MG(L)', 'MG(CR)', 'MG(S)']:
                    unburned_metal_sum += val
                    
            # Store results
            data.append({
                'MR': mr,
                'Isp': isp_vac,
                'Thrust_kN': thrust_kN,
                'T_Chamber': temps[0],
                'T_Exit': temps[2],
                'MW': mw,
                'Gamma': gamma,
                'Slag_Fraction': slag_sum,
                'Unburned_Metal_Fraction': unburned_metal_sum
            })

        df = pd.DataFrame(data)

        # 5. Visualizing the Fluid/Gas Properties
        for widget in plot_frame.winfo_children():
            widget.destroy()

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 9))
        
        # --- Plot 1: Performance (Isp and Temp) ---
        ax1.plot(df['MR'], df['Isp'], 'b-', label='Isp (s)')
        ax1.set_ylabel('Vacuum Isp (s)', color='b')
        ax_temp = ax1.twinx()
        ax_temp.plot(df['MR'], df['T_Chamber'], 'r--', label='Chamber T (K)')
        ax_temp.set_ylabel('Temp (K)', color='r')
        ax1.set_title("1. Specific Impulse & Thermal Profile")

        # --- Plot 2: Primary Motor Thrust ---
        ax2.plot(df['MR'], df['Thrust_kN'], 'm-', linewidth=2, label='Thrust (kN)')
        ax2.set_ylabel('Thrust (kN)', color='m')
        ax2.set_title(f"2. Thrust Output (@ {mdot} kg/s Mass Flow)")
        ax2.grid(True, linestyle=':', alpha=0.7) 

        # --- Plot 3: Condensed Phase (Slag vs Unburned Metal) ---
        ax3.plot(df['MR'], df['Slag_Fraction'], 'r-', label='Slag (Oxides)')
        ax3.plot(df['MR'], df['Unburned_Metal_Fraction'], 'g-', linewidth=2, label='Unburned Metal')
        ax3.set_ylabel('Mass Fraction', color='k')
        ax3.set_xlabel('Oxidizer-to-Fuel Ratio (O/F)')
        ax3.set_title("3. Condensed Phase Tracking (For CFD)")
        ax3.legend()

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        global global_df
        global_df = df
        res_label.config(text=f"Sweep Complete. {len(df)} points analyzed.")

    except Exception as e:
        messagebox.showerror("Error", str(e))

def export_csv():
    if 'global_df' in globals():
        global_df.to_csv("WBR_Research_Data.csv", index=False)
        messagebox.showinfo("Success", "Data exported to WBR_Research_Data.csv\nIncludes MW, Gamma, Slag, and Unburned Metal.")

# --- UI Setup ---
root = tk.Tk()
root.title("WBR Research & Slag Tracking Lab")
root.geometry("1100x700")

# Input Side
input_frame = ttk.LabelFrame(root, text="Configuration", padding=15)
input_frame.pack(side="left", fill="y", padx=10, pady=10)

ttk.Label(input_frame, text="Al Mass % (0-100):").pack(anchor="w")
al_entry = ttk.Entry(input_frame); al_entry.insert(0, "75"); al_entry.pack(fill="x", pady=5)

ttk.Label(input_frame, text="Chamber Pressure (psia):").pack(anchor="w")
pc_entry = ttk.Entry(input_frame); pc_entry.insert(0, "150"); pc_entry.pack(fill="x", pady=5)

ttk.Label(input_frame, text="Exit Area Ratio (eps):").pack(anchor="w")
eps_entry = ttk.Entry(input_frame); eps_entry.insert(0, "10"); eps_entry.pack(fill="x", pady=5)

ttk.Label(input_frame, text="Propellant Mass Flow (kg/s):").pack(anchor="w")
mdot_entry = ttk.Entry(input_frame); mdot_entry.insert(0, "0.5"); mdot_entry.pack(fill="x", pady=5)

frozen_var = tk.IntVar()
ttk.Checkbutton(input_frame, text="Use Frozen Flow", variable=frozen_var).pack(pady=10)

ttk.Button(input_frame, text="Run Analysis", command=run_research_analysis).pack(fill="x", pady=5)
ttk.Button(input_frame, text="Export CSV for CFD", command=export_csv).pack(fill="x", pady=5)

res_label = ttk.Label(input_frame, text="Ready", foreground="gray")
res_label.pack(pady=20)

# Plot Side
plot_frame = ttk.Frame(root)
plot_frame.pack(side="right", fill="both", expand=True)

root.mainloop()