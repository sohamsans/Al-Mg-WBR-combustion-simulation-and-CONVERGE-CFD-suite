import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import traceback

# rocketcea import - ensure you have rocketcea installed and accessible
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

def run_research_analysis():
    try:
        # --- 1. Read UI parameters (user enters Pc in psia as before) ---
        al_pct = float(al_entry.get())
        if not (0 <= al_pct <= 100):
            raise ValueError("Al mass % must be between 0 and 100")
        al_frac = al_pct / 100.0
        mg_frac = 1.0 - al_frac

        Pc_psia = float(pc_entry.get())        # user-provided chamber pressure in psia
        eps = float(eps_entry.get())          # expansion ratio (Aexit/At)
        use_frozen = bool(frozen_var.get())

        # NOTE: rocketcea historically uses Pc in psia for many methods; keep consistent.
        # If you later convert to SI, document conversions explicitly.

        # --- 2. (Optional) Register simple custom fuel/oxidizer using rocketcea helpers ---
        # The string format used by rocketcea's add_new_fuel/add_new_oxidizer can be strict.
        # Wrap in try/except to avoid crashing if the format isn't accepted.
        try:
            add_new_oxidizer('MyH2O', "ox     H2O(L)   H 2 O      1  wt%=100.0  h,cal=-68.32  t(k)=298.15")
        except Exception as e:
            # If add_new_oxidizer fails, it's not fatal — continue if 'H2O' already available
            print("Warning: add_new_oxidizer raised:", e)

        try:
            # Create a simple "fuel" alias. The exact rocketcea string format may vary by version.
            add_new_fuel('AlMg_Fuel', f"fuel   AlMg_Fuel   Al {al_frac*100:.6g}  Mg {mg_frac*100:.6g}  wt%=100.0  h,cal=0.0  t(k)=298.15")
        except Exception as e:
            print("Warning: add_new_fuel raised:", e)

        # --- 3. Initialize CEA object ---
        # If custom names failed, fallback to a generic CEA_Obj using standard fuel/ox names:
        try:
            cea = CEA_Obj(oxName='MyH2O', fuelName='AlMg_Fuel')
        except Exception:
            # Fallback: try using water and aluminum alone — user may need to provide real fuel/ox
            cea = CEA_Obj(oxName='H2O(L)', fuelName='AL')

        # --- 4. Sweep mixture ratio (MR) and collect data ---
        mr_range = np.linspace(1.5, 10.0, 40)
        rows = []
        for mr in mr_range:
            # call get_Isp with named arguments where possible
            try:
                isp_vac = cea.get_Isp(Pc=Pc_psia, MR=mr, eps=eps, frozen=use_frozen)
            except TypeError:
                # some rocketcea versions have different signature
                isp_vac = cea.get_Isp(Pc=Pc_psia, MR=mr, eps=eps)

            # Temperatures: try common method names and keep defensive
            try:
                temps = cea.get_Temperatures(Pc=Pc_psia, MR=mr, eps=eps, frozen=use_frozen)
            except Exception:
                # fallback - many implementations return (Tc, Tt, Te) or similar
                try:
                    temps = cea.get_Temperatures(Pc=Pc_psia, MR=mr, eps=eps)
                except Exception:
                    temps = (np.nan, np.nan, np.nan)

            # Mol weight & gamma (defensive)
            try:
                mw_gam = cea.get_Chamber_MolWt_gamma(Pc=Pc_psia, MR=mr, eps=eps)
                mw = mw_gam[0]
                gamma = mw_gam[1]
            except Exception:
                # fallback: try separate calls or set NaN
                mw = np.nan
                gamma = np.nan

            # Species mass fractions -> robust extraction
            slag_fraction = np.nan
            try:
                # Try common method names/arg orders
                try:
                    species_data = cea.get_SpeciesMassFractions(Pc_psia, mr, eps)
                except TypeError:
                    species_data = cea.get_SpeciesMassFractions(Pc=Pc_psia, MR=mr, eps=eps)
            except Exception:
                species_data = None

            # parse species_data if available (supports dicts or tuple(dict,...))
            if species_data:
                chamber_species = species_data[0] if isinstance(species_data, (list, tuple)) else species_data
                # robust iterate
                slag_sum = 0.0
                slag_targets = ['AL2O3', 'MGO', 'AL2O3(L)', 'AL2O3(S)', 'MGO(S)', 'MG(S)', 'AL(L)']
                for name, val in chamber_species.items():
                    try:
                        name_u = str(name).upper()
                        # If returned value is list-like, pick first element
                        v = val[0] if isinstance(val, (list, tuple, np.ndarray)) else float(val)
                        if any(target in name_u for target in slag_targets):
                            slag_sum += v
                    except Exception:
                        continue
                # normalize if percent
                slag_fraction = slag_sum/100.0 if slag_sum > 1.0 else slag_sum

            rows.append({
                'MR': float(mr),
                'Isp_vac_s': float(isp_vac) if np.isfinite(isp_vac) else np.nan,
                'T_chamber_K': float(temps[0]) if len(temps) > 0 else np.nan,
                'T_exit_K': float(temps[2]) if len(temps) > 2 else np.nan,
                'MolecularWeight': float(mw) if np.isfinite(mw) else np.nan,
                'Gamma': float(gamma) if np.isfinite(gamma) else np.nan,
                'SlagFraction': float(slag_fraction) if not np.isnan(slag_fraction) else np.nan
            })

        df = pd.DataFrame(rows)
        # Plot and display
        for w in plot_frame.winfo_children():
            w.destroy()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), tight_layout=True)
        ax1.plot(df['MR'], df['Isp_vac_s'], 'b-', label='Isp (s)')
        ax1.set_ylabel('Isp (s)', color='b')
        ax1.set_xlabel('MR')
        ax_temp = ax1.twinx()
        ax_temp.plot(df['MR'], df['T_chamber_K'], 'r--', label='Chamber T')
        ax_temp.set_ylabel('T (K)', color='r')
        ax1.set_title('Performance & Chamber Temp')

        ax2.plot(df['MR'], df['MolecularWeight'], 'g-', label='MolWt')
        ax2.set_ylabel('Mol Weight', color='g')
        ax_gamma = ax2.twinx()
        ax_gamma.plot(df['MR'], df['Gamma'], 'm:', label='Gamma')
        ax_gamma.set_ylabel('Gamma', color='m')
        ax2.set_xlabel('MR')
        ax2.set_title('Gas Properties (for CFD)')

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        # store & update UI
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

# --- UI setup ---
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