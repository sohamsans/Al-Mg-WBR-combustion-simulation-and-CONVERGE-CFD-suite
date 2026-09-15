import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import traceback
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

def run_research_analysis():
    if species_data: chamber_species = species_data[0] if isinstance(species_data, (list, tuple)) else species_data
    slag_sum = 0.0
    slag_targets = ['AL2O3', 'MGO', 'AL2O3(L)', 'AL2O3(S)', 'MGO(S)', 'MG(S)', 'AL(L)']
    for name, val in chamber_species.items():
        try:
            name_u = str(name).upper()
            v = val[0] if isinstance(val, (list, tuple, np.ndarray)) else float(val)  # ensure v is a float
            if any(target in name_u for target in slag_targets):
                slag_sum += v
        except Exception as e:
            print(f"Warning extracting slag data for {name}: {e}")  # Informative error message
    try:
        al_pct = float(al_entry.get())
        if not (0 <= al_pct <= 100):
            raise ValueError("Al mass % must be between 0 and 100")
        al_frac = al_pct / 100.0
        mg_frac = 1.0 - al_frac
        Pc_psia = float(pc_entry.get())
        eps = float(eps_entry.get())
        use_frozen = bool(frozen_var.get())


        # Register fuel/oxidizer
        add_new_oxidizer('MyH2O', "ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-68.32 t(k)=298.15")
        add_new_fuel('AlMg_Fuel', f"fuel AlMg_Fuel Al {al_frac*100:.6g} Mg {mg_frac*100:.6g} wt%=100.0 h,cal=0.0 t(k)=298.15")

        cea = CEA_Obj(oxName='MyH2O', fuelName='AlMg_Fuel')

        mr_range = np.linspace(1.5, 10.0, 40)
        rows = []
        for mr in mr_range:
            isp_vac = cea.get_Isp(Pc=Pc_psia, MR=mr, eps=eps, frozen=use_frozen)
            temps = cea.get_Temperatures(Pc=Pc_psia, MR=mr, eps=eps, frozen=use_frozen)
            mw_gam = cea.get_Chamber_MolWt_gamma(Pc=Pc_psia, MR=mr, eps=eps)
            species_data = cea.get_SpeciesMassFractions(Pc_psia, mr, eps)

            # Robust extraction of data
            slag_sum = sum([float(val[0]) for name, val in species_data[0].items() if any(target in name.upper() for target in ['AL2O3', 'MGO', 'AL'])])
            slag_fraction = slag_sum / 100.0 if slag_sum > 1.0 else slag_sum

            rows.append({
                'MR': mr,
                'Isp_vac_s': isp_vac,
                'T_chamber_K': temps[0] if len(temps) > 0 else np.nan,
                'T_exit_K': temps[2] if len(temps) > 2 else np.nan,
                'MolecularWeight': mw_gam[0],
                'Gamma': mw_gam[1],
                'SlagFraction': slag_fraction
            })

        df = pd.DataFrame(rows)

        for widget in plot_frame.winfo_children():
            widget.destroy()
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), tight_layout=True)
        ax1.plot(df['MR'], df['Isp_vac_s'], 'b-', label='Isp (s)')
        ax1.set_ylabel('Isp (s)', color='b')
        ax_temp = ax1.twinx()
        ax_temp.plot(df['MR'], df['T_chamber_K'], 'r--', label='Chamber T')
        ax_temp.set_ylabel('T (K)', color='r')
        ax1.set_title('Performance & Chamber Temp')

        ax2.plot(df['MR'], df['MolecularWeight'], 'g-', label='MolWeight')
        ax2.set_ylabel('Molecular Weight', color='g')
        ax_gamma = ax2.twinx()
        ax_gamma.plot(df['MR'], df['Gamma'], 'm:', label='Gamma')
        ax_gamma.set_ylabel('Gamma', color='m')
        ax2.set_xlabel('MR')
        ax2.set_title('Gas Properties (for CFD)')

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        global global_df
        global_df = df
        res_label.config(text=f"Sweep complete: {len(df)} points.", foreground="green")
    except Exception as exc:
        tb = traceback.format_exc()
        messagebox.showerror("Error running analysis", f"{exc}\n\n{tb}")
        res_label.config(text="Error", foreground="red")

def export_csv():
    if 'global_df' not in globals():
        messagebox.showwarning("No data", "Run the sweep first.")
        return
    fname = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files","*.csv")])
    if not fname:
        return
    global_df.to_csv(fname, index=False)
    messagebox.showinfo("Saved", f"Exported CSV to: {fname}")

# --- UI Setup ---
root = tk.Tk()
root.title("WBR Research & Slag Tracking Lab")
root.geometry("1100x700")

input_frame = ttk.LabelFrame(root, text="Configuration", padding=12)
input_frame.pack(side="left", fill="y", padx=8, pady=8)

ttk.Label(input_frame, text="Al mass % (0-100):").pack(anchor="w")
al_entry = ttk.Entry(input_frame); al_entry.insert(0, "75"); al_entry.pack(fill="x", pady=4)
ttk.Label(input_frame, text="Chamber Pressure (psia):").pack(anchor="w")
pc_entry = ttk.Entry(input_frame); pc_entry.insert(0, "200"); pc_entry.pack(fill="x", pady=4)
ttk.Label(input_frame, text="Exit area ratio (eps):").pack(anchor="w")
eps_entry = ttk.Entry(input_frame); eps_entry.insert(0, "8"); eps_entry.pack(fill="x", pady=4)

frozen_var = tk.IntVar()
ttk.Checkbutton(input_frame, text="Use frozen flow", variable=frozen_var).pack(pady=8)
ttk.Button(input_frame, text="Run analysis", command=run_research_analysis).pack(fill="x", pady=4)
ttk.Button(input_frame, text="Export CSV for CFD", command=export_csv).pack(fill="x", pady=4)

res_label = ttk.Label(input_frame, text="Ready", foreground="gray")
res_label.pack(pady=20)

plot_frame = ttk.Frame(root)
plot_frame.pack(side="right", fill="both", expand=True)

root.mainloop()