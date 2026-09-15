"""
Unified Final Research Laboratory & Parametric Sweep Suite (CEA & WBR Physics)
==============================================================================
Monolithic Tkinter GUI Application providing four specialized tabs:
  1. Thermochemical Equilibrium (NASA CEA) & Water-Ramjet Performance Sweeps
  2. 1D Droplet Trajectory & Expulsion Efficiency (eta_expulsion) Modeling
  3. Multi-Variable Parametric Sweeps (Al:Mg, d_p0, L_c, P_c) & Optimization
  4. CFD Data Exporter & Technical Visualizer
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import pandas as pd
import traceback
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

# Import high-fidelity 1D primary combustion physics model
from primary_combustion_model import simulate_1d_primary_combustor, get_alloy_properties
from bates_converge_cfd_exporter import simulate_bates_grain_regression, export_converge_cfd_inputs

class UnifiedWBRResearchLab(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Unified Al–Mg Hydro-Reactive WBR Research & Simulation Suite")
        self.geometry("1280x850")
        
        self.global_df = None
        self.sweep_df = None
        self.bates_df = None
        
        self._build_layout()
        
    def _build_layout(self):
        # Top Control & Parameter Panel
        control_frame = ttk.LabelFrame(self, text="Propellant & Combustor Parameters", padding=10)
        control_frame.pack(side="top", fill="x", padx=10, pady=5)
        
        # Row 1: Formulation & Geometry
        ttk.Label(control_frame, text="Al Mass % (0-100):").grid(row=0, column=0, sticky="w", padx=5)
        self.al_entry = ttk.Entry(control_frame, width=8)
        self.al_entry.insert(0, "75")
        self.al_entry.grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(control_frame, text="HTPB Binder %:").grid(row=0, column=2, sticky="w", padx=5)
        self.htpb_entry = ttk.Entry(control_frame, width=8)
        self.htpb_entry.insert(0, "15")
        self.htpb_entry.grid(row=0, column=3, padx=5, pady=2)

        ttk.Label(control_frame, text="Chamber P_c (psia):").grid(row=0, column=4, sticky="w", padx=5)
        self.pc_entry = ttk.Entry(control_frame, width=8)
        self.pc_entry.insert(0, "200")
        self.pc_entry.grid(row=0, column=5, padx=5, pady=2)

        ttk.Label(control_frame, text="Area Ratio (eps):").grid(row=0, column=6, sticky="w", padx=5)
        self.eps_entry = ttk.Entry(control_frame, width=8)
        self.eps_entry.insert(0, "8")
        self.eps_entry.grid(row=0, column=7, padx=5, pady=2)

        # Row 2: Particle & Primary Geometry
        ttk.Label(control_frame, text="Particle d_p0 (µm):").grid(row=1, column=0, sticky="w", padx=5)
        self.dp_entry = ttk.Entry(control_frame, width=8)
        self.dp_entry.insert(0, "50")
        self.dp_entry.grid(row=1, column=1, padx=5, pady=2)

        ttk.Label(control_frame, text="Chamber Length L_c (m):").grid(row=1, column=2, sticky="w", padx=5)
        self.lc_entry = ttk.Entry(control_frame, width=8)
        self.lc_entry.insert(0, "0.20")
        self.lc_entry.grid(row=1, column=3, padx=5, pady=2)

        ttk.Label(control_frame, text="Primary O/F Ratio:").grid(row=1, column=4, sticky="w", padx=5)
        self.of_entry = ttk.Entry(control_frame, width=8)
        self.of_entry.insert(0, "0.25")
        self.of_entry.grid(row=1, column=5, padx=5, pady=2)

        self.frozen_var = tk.IntVar()
        ttk.Checkbutton(control_frame, text="Frozen Flow", variable=self.frozen_var).grid(row=1, column=6, columnspan=2, padx=5)

        # Status Bar
        self.status_label = ttk.Label(self, text="Ready", foreground="gray", font=("Arial", 10, "italic"))
        self.status_label.pack(side="bottom", fill="x", padx=10, pady=5)

        # Main Tabbed Interface
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(side="top", fill="both", expand=True, padx=10, pady=5)

        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        self.tab3 = ttk.Frame(self.notebook)
        self.tab4 = ttk.Frame(self.notebook)
        self.tab5 = ttk.Frame(self.notebook)

        self.notebook.add(self.tab1, text="1. NASA CEA Equilibrium & Water Sweep")
        self.notebook.add(self.tab2, text="2. 1D Droplet Trajectory & Expulsion Model")
        self.notebook.add(self.tab3, text="3. Parametric Sweeps & Optimization")
        self.notebook.add(self.tab4, text="4. CFD Exporter & Data Viewer")
        self.notebook.add(self.tab5, text="5. BATES Grain & CONVERGE CFD Setup")

        self._setup_tab1()
        self._setup_tab2()
        self._setup_tab3()
        self._setup_tab4()
        self._setup_tab5()

    def _setup_tab1(self):
        btn_frame = ttk.Frame(self.tab1, padding=5)
        btn_frame.pack(side="top", fill="x")
        
        ttk.Button(btn_frame, text="Run Water-Ramjet Sweep (H2O)", command=self.run_cea_water_sweep).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Export CSV Data", command=self.export_csv).pack(side="left", padx=5)

        self.fig1, (self.ax1_1, self.ax1_2) = plt.subplots(2, 1, figsize=(8, 6))
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.tab1)
        self.canvas1.get_tk_widget().pack(fill="both", expand=True)

    def _setup_tab2(self):
        btn_frame = ttk.Frame(self.tab2, padding=5)
        btn_frame.pack(side="top", fill="x")

        ttk.Button(btn_frame, text="Run 1D Primary Expulsion Model", command=self.run_primary_expulsion_model).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Export Trajectory Data", command=self.export_csv).pack(side="left", padx=5)

        self.fig2, (self.ax2_1, self.ax2_2) = plt.subplots(2, 1, figsize=(8, 6))
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.tab2)
        self.canvas2.get_tk_widget().pack(fill="both", expand=True)

    def _setup_tab3(self):
        btn_frame = ttk.Frame(self.tab3, padding=5)
        btn_frame.pack(side="top", fill="x")

        ttk.Button(btn_frame, text="Run Multi-Variable Parametric Sweep", command=self.run_parametric_optimization_sweep).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Export Optimization CSV", command=self.export_sweep_csv).pack(side="left", padx=5)

        self.fig3, (self.ax3_1, self.ax3_2) = plt.subplots(1, 2, figsize=(9, 5))
        self.canvas3 = FigureCanvasTkAgg(self.fig3, master=self.tab3)
        self.canvas3.get_tk_widget().pack(fill="both", expand=True)

    def _setup_tab4(self):
        btn_frame = ttk.Frame(self.tab4, padding=5)
        btn_frame.pack(side="top", fill="x")

        ttk.Button(btn_frame, text="Generate CFD Propellant_Data.dat Header", command=self.generate_cfd_header).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Export Complete CSV", command=self.export_csv).pack(side="left", padx=5)

        self.text_preview = tk.Text(self.tab4, wrap="none", font=("Courier", 10))
        self.text_preview.pack(fill="both", expand=True, padx=5, pady=5)

    def _setup_tab5(self):
        btn_frame = ttk.Frame(self.tab5, padding=5)
        btn_frame.pack(side="top", fill="x")

        ttk.Button(btn_frame, text="Run BATES Grain Regression Simulation", command=self.run_bates_simulation).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Export CONVERGE CFD Boundary & Thermo Files", command=self.export_converge_files).pack(side="left", padx=5)

        self.fig5, (self.ax5_1, self.ax5_2) = plt.subplots(2, 1, figsize=(8, 6))
        self.canvas5 = FigureCanvasTkAgg(self.fig5, master=self.tab5)
        self.canvas5.get_tk_widget().pack(fill="both", expand=True)

    def _get_float(self, entry, name, min_val, max_val, default_val):
        try:
            val = float(entry.get())
            if val < min_val or val > max_val:
                messagebox.showwarning("Input Out of Range", f"{name} ({val}) is outside recommended range [{min_val}, {max_val}]. Clamping value.")
                val = float(np.clip(val, min_val, max_val))
                entry.delete(0, tk.END)
                entry.insert(0, str(val))
            return val
        except (ValueError, TypeError):
            messagebox.showwarning("Invalid Input", f"Invalid input for {name}. Using default {default_val}.")
            entry.delete(0, tk.END)
            entry.insert(0, str(default_val))
            return default_val

    def run_cea_water_sweep(self):
        try:
            al_pct = self._get_float(self.al_entry, "Al Mass %", 0.0, 100.0, 75.0)
            mg_pct = 100.0 - al_pct
            pc = self._get_float(self.pc_entry, "Chamber Pressure", 10.0, 5000.0, 200.0)
            eps = self._get_float(self.eps_entry, "Area Ratio", 1.0, 100.0, 8.0)
            sub_frozen = bool(self.frozen_var.get())

            add_new_oxidizer('MyH2O', "ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-3788.5 t(k)=298.15")
            add_new_fuel('AlMg_Fuel', f"fuel AlMg_Fuel Al {al_pct/100:.4f} Mg {mg_pct/100:.4f} wt%=100.0 h,cal=0.0 t(k)=298.15")

            cea = CEA_Obj(oxName='MyH2O', fuelName='AlMg_Fuel')
            mr_range = np.linspace(1.5, 10.0, 40)
            data = []

            for mr in mr_range:
                isp_vac = cea.get_Isp(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
                temps = cea.get_Temperatures(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
                mw_gam = cea.get_Chamber_MolWt_gamma(Pc=pc, MR=mr, eps=eps)

                data.append({
                    'MR': mr,
                    'Isp_vac_s': isp_vac,
                    'T_Chamber_K': temps[0] if len(temps)>0 else np.nan,
                    'T_Exit_K': temps[2] if len(temps)>2 else np.nan,
                    'MolecularWeight': mw_gam[0],
                    'Gamma': mw_gam[1]
                })

            df = pd.DataFrame(data)
            self.global_df = df

            self.ax1_1.clear()
            self.ax1_2.clear()

            self.ax1_1.plot(df['MR'], df['Isp_vac_s'], 'b-', label='Isp (s)')
            self.ax1_1.set_ylabel('Vacuum Isp (s)', color='b')
            ax_temp = self.ax1_1.twinx()
            ax_temp.plot(df['MR'], df['T_Chamber_K'], 'r--', label='Chamber T (K)')
            ax_temp.set_ylabel('Temp (K)', color='r')
            self.ax1_1.set_title("Water-Ramjet Performance & Thermal Profile")
            self.ax1_1.grid(True, linestyle=':', alpha=0.6)

            self.ax1_2.plot(df['MR'], df['MolecularWeight'], 'g-', label='MW (g/mol)')
            self.ax1_2.set_ylabel('Molecular Weight', color='g')
            ax_gam = self.ax1_2.twinx()
            ax_gam.plot(df['MR'], df['Gamma'], 'm:', label='Gamma')
            ax_gam.set_ylabel('Gamma', color='m')
            self.ax1_2.set_xlabel('Water-to-Fuel Ratio (O/F)')
            self.ax1_2.set_title("Gas Thermochemical Properties for CFD")
            self.ax1_2.grid(True, linestyle=':', alpha=0.6)

            self.fig1.tight_layout()
            self.canvas1.draw()

            self._update_text_preview(df)
            self.status_label.config(text=f"CEA Water Sweep Complete ({len(df)} points).", foreground="green")

        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text=f"Error: {e}", foreground="red")

    def run_primary_expulsion_model(self):
        try:
            al_pct = self._get_float(self.al_entry, "Al Mass %", 0.0, 100.0, 75.0)
            htpb_pct = self._get_float(self.htpb_entry, "HTPB Binder %", 1.0, 50.0, 15.0)
            pc = self._get_float(self.pc_entry, "Chamber Pressure", 10.0, 5000.0, 200.0)
            dp_um = self._get_float(self.dp_entry, "Particle Size d_p0", 1.0, 500.0, 50.0)
            lc_m = self._get_float(self.lc_entry, "Chamber Length L_c", 0.01, 2.0, 0.20)
            primary_of = self._get_float(self.of_entry, "Primary O/F Ratio", 0.01, 2.0, 0.25)

            df_traj, metrics = simulate_1d_primary_combustor(
                al_pct=al_pct,
                htpb_pct=htpb_pct,
                d_p0_um=dp_um,
                pc_psia=pc,
                primary_of=primary_of,
                chamber_length_m=lc_m
            )

            self.global_df = df_traj

            self.ax2_1.clear()
            self.ax2_2.clear()

            self.ax2_1.plot(df_traj['x_m'], df_traj['d_p_um'], 'b-', label='d_p (µm)')
            self.ax2_1.set_ylabel('Droplet Diameter (µm)', color='b')
            ax_eff = self.ax2_1.twinx()
            ax_eff.plot(df_traj['x_m'], df_traj['expulsion_eff_pct'], 'g--', label='Expulsion Eff (%)')
            ax_eff.axhline(90.0, color='r', linestyle=':', linewidth=2, label='Target Threshold (90%)')
            ax_eff.set_ylabel('Expulsion Efficiency (%)', color='g')
            self.ax2_1.set_title('Primary Combustion: Droplet Size & Metal Expulsion Efficiency')
            self.ax2_1.grid(True, linestyle=':', alpha=0.6)

            self.ax2_2.plot(df_traj['x_m'], df_traj['u_p_ms'], 'm-', label='Particle Vel u_p (m/s)')
            self.ax2_2.plot(df_traj['x_m'], df_traj['u_g_ms'], 'c--', label='Gas Vel u_g (m/s)')
            self.ax2_2.set_xlabel('Combustor Position x (m)')
            self.ax2_2.set_ylabel('Velocity (m/s)')
            self.ax2_2.set_title('Particle Drag Acceleration & Gas Expansion Profile')
            self.ax2_2.grid(True, linestyle=':', alpha=0.6)
            self.ax2_2.legend()

            self.fig2.tight_layout()
            self.canvas2.draw()

            self._update_text_preview(df_traj)

            eff_val = metrics['exit_expulsion_eff_pct']
            color = "green" if eff_val >= 90.0 else "orange"
            self.status_label.config(
                text=f"Expulsion Eff @ Exit: {eff_val:.1f}% | Slag Accumulation: {metrics['total_slag_mass_pct']:.2f}%",
                foreground=color
            )

        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text=f"Error: {e}", foreground="red")

    def run_parametric_optimization_sweep(self):
        try:
            al_list = [25.0, 50.0, 75.0, 90.0]
            dp_list = [15.0, 30.0, 50.0, 75.0, 100.0]
            length_list = np.linspace(0.05, 0.50, 8)
            pc = float(self.pc_entry.get())
            primary_of = float(self.of_entry.get())

            results = []
            for al in al_list:
                for dp in dp_list:
                    for Lc in length_list:
                        _, m = simulate_1d_primary_combustor(
                            al_pct=al,
                            d_p0_um=dp,
                            pc_psia=pc,
                            primary_of=primary_of,
                            chamber_length_m=Lc
                        )
                        results.append({
                            'Al_pct': al,
                            'd_p0_um': dp,
                            'L_c_m': Lc,
                            'Expulsion_Eff_pct': m['exit_expulsion_eff_pct'],
                            'Slag_Pct': m['total_slag_mass_pct'],
                            'Target_>=90_Met': m['target_met']
                        })

            df_sweep = pd.DataFrame(results)
            self.sweep_df = df_sweep

            self.ax3_1.clear()
            self.ax3_2.clear()

            df_75 = df_sweep[df_sweep['Al_pct'] == 75.0]
            for dp in sorted(df_75['d_p0_um'].unique()):
                sub = df_75[df_75['d_p0_um'] == dp]
                self.ax3_1.plot(sub['L_c_m'], sub['Expulsion_Eff_pct'], '-o', label=f'd_p0={dp}µm')
                self.ax3_2.plot(sub['L_c_m'], sub['Slag_Pct'], '-s', label=f'd_p0={dp}µm')

            self.ax3_1.axhline(90.0, color='r', linestyle='--', label='Target 90%')
            self.ax3_1.set_xlabel('Chamber Length L_c (m)')
            self.ax3_1.set_ylabel('Expulsion Efficiency (%)')
            self.ax3_1.set_title('Expulsion Efficiency vs L_c (75% Al)')
            self.ax3_1.grid(True, linestyle=':', alpha=0.6)
            self.ax3_1.legend()

            self.ax3_2.set_xlabel('Chamber Length L_c (m)')
            self.ax3_2.set_ylabel('Slag Mass Fraction (%)')
            self.ax3_2.set_title('Slag Accumulation vs L_c')
            self.ax3_2.grid(True, linestyle=':', alpha=0.6)
            self.ax3_2.legend()

            self.fig3.tight_layout()
            self.canvas3.draw()

            self._update_text_preview(df_sweep)
            valid_cnt = len(df_sweep[df_sweep['Target_>=90_Met'] == True])
            self.status_label.config(text=f"Optimization Sweep Complete. Valid Configurations (>=90% Expulsion): {valid_cnt}/{len(df_sweep)}", foreground="green")

        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text=f"Error: {e}", foreground="red")

    def run_bates_simulation(self):
        try:
            df_bates = simulate_bates_grain_regression()
            self.bates_df = df_bates
            self.global_df = df_bates

            self.ax5_1.clear()
            self.ax5_2.clear()

            self.ax5_1.plot(df_bates['time_s'], df_bates['P_c_MPa'], 'r-', label='P_c (MPa)')
            self.ax5_1.set_ylabel('Chamber Pressure (MPa)', color='r')
            ax_rb = self.ax5_1.twinx()
            ax_rb.plot(df_bates['time_s'], df_bates['r_b_mm_s'], 'b--', label='r_b (mm/s)')
            ax_rb.set_ylabel('Burn Rate (mm/s)', color='b')
            self.ax5_1.set_title('BATES Grain Internal Ballistics: Pressure & Burn Rate History')
            self.ax5_1.grid(True, linestyle=':', alpha=0.6)

            self.ax5_2.plot(df_bates['time_s'], df_bates['A_b_m2'], 'g-', label='A_b (m²)')
            self.ax5_2.set_ylabel('Burning Area A_b (m²)', color='g')
            ax_flux = self.ax5_2.twinx()
            ax_flux.plot(df_bates['time_s'], df_bates['mass_flux_wall_kg_m2s'], 'm--', label='Mass Flux (kg/m²s)')
            ax_flux.set_ylabel('Wall Mass Flux [kg/(m² s)]', color='m')
            self.ax5_2.set_xlabel('Time t (s)')
            self.ax5_2.set_title('BATES Grain Regression & CONVERGE Wall Inflow Mass Flux')
            self.ax5_2.grid(True, linestyle=':', alpha=0.6)

            self.fig5.tight_layout()
            self.canvas5.draw()

            self._update_text_preview(df_bates)
            self.status_label.config(
                text=f"BATES Simulation Complete. Burn Time: {df_bates['time_s'].iloc[-1]:.2f}s | Peak P_c: {df_bates['P_c_MPa'].max():.2f} MPa",
                foreground="green"
            )
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text=f"Error: {e}", foreground="red")

    def export_converge_files(self):
        if self.bates_df is None:
            messagebox.showwarning("No Data", "Run BATES Grain Regression simulation first!")
            return
        try:
            export_converge_cfd_inputs(self.bates_df)
            messagebox.showinfo("Export Success", "CONVERGE CFD boundary files successfully created in project directory!\n\nFiles Generated:\n- converge_bates_boundary.in\n- converge_inflow_mass_flux.dat\n- converge_thermo.dat\n- bates_internal_ballistics.csv")
        except Exception as e:
            messagebox.showerror("Export Error", str(e))

    def generate_cfd_header(self):
        if self.global_df is None:
            messagebox.showwarning("No Data", "Run a simulation sweep first!")
            return

        header = f"# ==========================================================================\n"
        header += f"# CFD-Ready Propellant Table (NASA CEA & 1D Droplet Dynamics)\n"
        header += f"# Generated by Unified Al-Mg Hydro-Reactive WBR Research Suite\n"
        header += f"# ==========================================================================\n"
        header += self.global_df.head(20).to_string(index=False)
        header += f"\n\n# Total Rows: {len(self.global_df)}\n"

        self.text_preview.delete("1.0", tk.END)
        self.text_preview.insert(tk.END, header)

    def _update_text_preview(self, df):
        self.text_preview.delete("1.0", tk.END)
        self.text_preview.insert(tk.END, df.to_string(index=False))

    def export_csv(self):
        if self.global_df is None:
            messagebox.showwarning("No Data", "Run a simulation sweep first!")
            return
        fname = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if fname:
            self.global_df.to_csv(fname, index=False)
            messagebox.showinfo("Exported", f"Data exported successfully to: {fname}")

    def export_sweep_csv(self):
        if self.sweep_df is None:
            messagebox.showwarning("No Data", "Run optimization sweep first!")
            return
        fname = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if fname:
            self.sweep_df.to_csv(fname, index=False)
            messagebox.showinfo("Exported", f"Optimization sweep exported to: {fname}")

if __name__ == '__main__':
    app = UnifiedWBRResearchLab()
    app.mainloop()

