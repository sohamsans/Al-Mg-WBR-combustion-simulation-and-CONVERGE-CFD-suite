import tkinter as tk
from tkinter import ttk, messagebox
from rocketcea.cea_obj import CEA_Obj
import rocketcea.cea_obj as cea_obj

# Robust way to get propellants across different versions
try:
    # Try to get from the internal data dictionary
    all_props = sorted(cea_obj.propDict.keys())
except AttributeError:
    # Fallback common presets if dictionary access fails
    all_props = ["LOX", "LH2", "RP1", "CH4", "NTO", "MMH", "UDMH", "A50", "GH2", "GO2"]

def run_simulation():
    try:
        ox = ox_var.get()
        fuel = fuel_var.get()
        pc = float(pc_entry.get())
        mr = float(mr_entry.get())
        eps = float(eps_entry.get())

        # Initialize
        isp_obj = CEA_Obj(oxName=ox, fuelName=fuel)
        
        # Calculate Vacuum Isp
        isp_vac = isp_obj.get_Isp(Pc=pc, MR=mr, eps=eps)
        
        # Calculate Sea Level Isp 
        # Using 14.7 as the 4th positional argument if 'pamb' failed
        isp_sl = isp_obj.get_Isp(pc, mr, eps, 14.7) 
        
        temps = isp_obj.get_Temperatures(Pc=pc, MR=mr, eps=eps)
        t_cham = temps[0]

        # Update Results
        res_isp_vac.config(text=f"{isp_vac:.2f} s", foreground="green")
        res_isp_sl.config(text=f"{isp_sl:.2f} s", foreground="blue")
        res_temp.config(text=f"{t_cham:.2f} K", foreground="red")
        
    except Exception as e:
        messagebox.showerror("Error", f"Simulation failed: {e}")

# --- GUI Setup ---
root = tk.Tk()
root.title("NASA CEA Simulation Tool")
root.geometry("450x500")

frame = ttk.Frame(root, padding="20")
frame.pack(fill="both", expand=True)

# Selection Menus
ttk.Label(frame, text="Oxidizer:").grid(row=0, column=0, sticky="w")
ox_var = tk.StringVar(value="LOX")
ox_menu = ttk.Combobox(frame, textvariable=ox_var, values=all_props)
ox_menu.grid(row=0, column=1, pady=5, sticky="ew")

ttk.Label(frame, text="Fuel:").grid(row=1, column=0, sticky="w")
fuel_var = tk.StringVar(value="RP1")
fuel_menu = ttk.Combobox(frame, textvariable=fuel_var, values=all_props)
fuel_menu.grid(row=1, column=1, pady=5, sticky="ew")

# Inputs
ttk.Label(frame, text="Chamber Pressure (psia):").grid(row=2, column=0, sticky="w")
pc_entry = ttk.Entry(frame)
pc_entry.insert(0, "1000")
pc_entry.grid(row=2, column=1, pady=5, sticky="ew")

ttk.Label(frame, text="Mixture Ratio (O/F):").grid(row=3, column=0, sticky="w")
mr_entry = ttk.Entry(frame)
mr_entry.insert(0, "2.56")
mr_entry.grid(row=3, column=1, pady=5, sticky="ew")

ttk.Label(frame, text="Expansion Ratio (ε):").grid(row=4, column=0, sticky="w")
eps_entry = ttk.Entry(frame)
eps_entry.insert(0, "15.0")
eps_entry.grid(row=4, column=1, pady=5, sticky="ew")

# Action
run_btn = ttk.Button(frame, text="Calculate Performance", command=run_simulation)
run_btn.grid(row=5, column=0, columnspan=2, pady=20)

# Detailed Results
ttk.Label(frame, text="Vacuum Isp:").grid(row=6, column=0, sticky="w")
res_isp_vac = ttk.Label(frame, text="--", font=("Arial", 10, "bold"))
res_isp_vac.grid(row=6, column=1, sticky="w")

ttk.Label(frame, text="Sea Level Isp:").grid(row=7, column=0, sticky="w")
res_isp_sl = ttk.Label(frame, text="--", font=("Arial", 10, "bold"))
res_isp_sl.grid(row=7, column=1, sticky="w")

ttk.Label(frame, text="Chamber Temp:").grid(row=8, column=0, sticky="w")
res_temp = ttk.Label(frame, text="--", font=("Arial", 10, "bold"))
res_temp.grid(row=8, column=1, sticky="w")

root.mainloop()