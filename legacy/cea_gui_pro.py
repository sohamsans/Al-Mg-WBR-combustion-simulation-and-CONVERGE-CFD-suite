import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rocketcea.cea_obj import CEA_Obj
import rocketcea.cea_obj as cea_obj
import numpy as np

# Propellant list setup
try:
    all_props = sorted(cea_obj.propDict.keys())
except:
    all_props = ["LOX", "LH2", "RP1", "CH4", "NTO", "MMH", "UDMH", "A50"]

def run_sweep():
    try:
        ox = ox_var.get()
        fuel = fuel_var.get()
        pc = float(pc_entry.get())
        eps = float(eps_entry.get())

        # Create range of Mixture Ratios (from 1 to 10)
        mr_range = np.linspace(1.0, 10.0, 50)
        isp_list = []
        temp_list = []

        isp_obj = CEA_Obj(oxName=ox, fuelName=fuel)

        for mr in mr_range:
            isp_list.append(isp_obj.get_Isp(Pc=pc, MR=mr, eps=eps))
            temp_list.append(isp_obj.get_Temperatures(Pc=pc, MR=mr, eps=eps)[0])

        # Create Plotting Window
        plot_window = tk.Toplevel(root)
        plot_window.title(f"Performance Sweep: {ox}/{fuel}")

        fig, ax1 = plt.subplots(figsize=(6, 4))

        # Plot Isp
        ax1.set_xlabel('Mixture Ratio (O/F)')
        ax1.set_ylabel('Vac Isp (s)', color='tab:blue')
        ax1.plot(mr_range, isp_list, color='tab:blue', lw=2, label='Isp')
        ax1.tick_params(axis='y', labelcolor='tab:blue')

        # Create second axis for Temperature
        ax2 = ax1.twinx()
        ax2.set_ylabel('Chamber Temp (K)', color='tab:red')
        ax2.plot(mr_range, temp_list, color='tab:red', linestyle='--', label='Temp')
        ax2.tick_params(axis='y', labelcolor='tab:red')

        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    except Exception as e:
        messagebox.showerror("Error", f"Sweep failed: {e}")

# --- Keep your previous run_simulation function here ---
def run_simulation():
    try:
        isp_obj = CEA_Obj(oxName=ox_var.get(), fuelName=fuel_var.get())
        pc, mr, eps = float(pc_entry.get()), float(mr_entry.get()), float(eps_entry.get())
        
        # Positional arguments for robustness
        vac = isp_obj.get_Isp(pc, mr, eps)
        sl = isp_obj.get_Isp(pc, mr, eps, 14.7)
        temp = isp_obj.get_Temperatures(pc, mr, eps)[0]

        res_isp_vac.config(text=f"{vac:.2f} s", foreground="green")
        res_isp_sl.config(text=f"{sl:.2f} s", foreground="blue")
        res_temp.config(text=f"{temp:.2f} K", foreground="red")
    except Exception as e:
        messagebox.showerror("Error", f"Calc failed: {e}")

# --- GUI Layout (Main) ---
root = tk.Tk()
root.title("NASA CEA Pro Tool")
root.geometry("450x600")

frame = ttk.Frame(root, padding="20")
frame.pack(fill="both", expand=True)

# (Re-use your previous UI labels/entries here)
ttk.Label(frame, text="Oxidizer:").grid(row=0, column=0, sticky="w")
ox_var = tk.StringVar(value="LOX")
ttk.Combobox(frame, textvariable=ox_var, values=all_props).grid(row=0, column=1, pady=5)

ttk.Label(frame, text="Fuel:").grid(row=1, column=0, sticky="w")
fuel_var = tk.StringVar(value="RP1")
ttk.Combobox(frame, textvariable=fuel_var, values=all_props).grid(row=1, column=1, pady=5)

ttk.Label(frame, text="Chamber Pressure:").grid(row=2, column=0, sticky="w")
pc_entry = ttk.Entry(frame); pc_entry.insert(0, "1000"); pc_entry.grid(row=2, column=1, pady=5)

ttk.Label(frame, text="Mixture Ratio:").grid(row=3, column=0, sticky="w")
mr_entry = ttk.Entry(frame); mr_entry.insert(0, "2.56"); mr_entry.grid(row=3, column=1, pady=5)

ttk.Label(frame, text="Expansion Ratio:").grid(row=4, column=0, sticky="w")
eps_entry = ttk.Entry(frame); eps_entry.insert(0, "15.0"); eps_entry.grid(row=4, column=1, pady=5)

# Buttons
ttk.Button(frame, text="Single Point Calc", command=run_simulation).grid(row=5, column=0, columnspan=2, pady=10, sticky="ew")
ttk.Button(frame, text="Generate Performance Graph", command=run_sweep).grid(row=6, column=0, columnspan=2, pady=10, sticky="ew")

# Result Display
res_isp_vac = ttk.Label(frame, text="--")
res_isp_vac.grid(row=7, column=1)
ttk.Label(frame, text="Vac Isp:").grid(row=7, column=0)

res_isp_sl = ttk.Label(frame, text="--")
res_isp_sl.grid(row=8, column=1)
ttk.Label(frame, text="SL Isp:").grid(row=8, column=0)

res_temp = ttk.Label(frame, text="--")
res_temp.grid(row=9, column=1)
ttk.Label(frame, text="Chamber Temp:").grid(row=9, column=0)

root.mainloop()