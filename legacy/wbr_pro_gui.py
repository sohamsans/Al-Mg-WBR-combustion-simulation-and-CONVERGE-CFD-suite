import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
import numpy as np

def run_simulation():
    try:
        # 1. Get Parameters from UI
        al_frac = float(al_entry.get()) / 100.0
        mg_frac = 1.0 - al_frac
        pc = float(pc_entry.get())
        eps = float(eps_entry.get())
        
        # 2. Define Custom Propellants (Manual Card method for 100% reliability)
        # Water Oxidizer
        ox_card = "ox MyH2O H 2 O 1 wt%=100.0 h,cal=-3788.5 t(k)=298.15"
        add_new_oxidizer('MyH2O', ox_card)
        
        # Al-Mg Composite Fuel
        fuel_card = f"fuel AlMg_Fuel Al {al_frac:.4f} Mg {mg_frac:.4f} wt%=100.0 h,cal=0.0 t(k)=298.15"
        add_new_fuel('AlMg_Fuel', fuel_card)
        
        # 3. Initialize CEA
        isp_obj = CEA_Obj(oxName='MyH2O', fuelName='AlMg_Fuel')
        
        # 4. Sweep MR (Water-to-Fuel ratio)
        mr_range = np.linspace(1.5, 8.0, 50)
        isp_list = []
        temp_list = []
        
        for mr in mr_range:
            isp_list.append(isp_obj.get_Isp(Pc=pc, MR=mr, eps=eps))
            temp_list.append(isp_obj.get_Temperatures(Pc=pc, MR=mr, eps=eps)[0])
            
        # 5. Plotting
        for widget in plot_frame.winfo_children():
            widget.destroy()
            
        fig, ax1 = plt.subplots(figsize=(5, 4), dpi=100)
        ax1.plot(mr_range, isp_list, 'b-', label='Isp')
        ax1.set_xlabel('Water/Fuel Ratio (Mass)')
        ax1.set_ylabel('Vacuum Isp (s)', color='b')
        
        ax2 = ax1.twinx()
        ax2.plot(mr_range, temp_vals := temp_list, 'r--', label='Temp')
        ax2.set_ylabel('Temp (K)', color='r')
        
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()
        
        # Update Single Point Labels (at 3.0 MR for reference)
        ref_isp = isp_obj.get_Isp(pc, 3.0, eps)
        res_label.config(text=f"Ref Isp @ 3.0 MR: {ref_isp:.2f} s")

    except Exception as e:
        messagebox.showerror("Simulation Error", str(e))

# --- UI Setup ---
root = tk.Tk()
root.title("WBR Analysis Tool: Al-Mg Composite")
root.geometry("900x600")

# Input Panel
input_frame = ttk.LabelFrame(root, text="Engine Parameters", padding=10)
input_frame.pack(side="left", fill="y", padx=10, pady=10)

ttk.Label(input_frame, text="Aluminum % (Rest is Mg):").pack(anchor="w")
al_entry = ttk.Entry(input_frame); al_entry.insert(0, "70"); al_entry.pack(fill="x", pady=5)

ttk.Label(input_frame, text="Combustor Pressure (psia):").pack(anchor="w")
pc_entry = ttk.Entry(input_frame); pc_entry.insert(0, "150"); pc_entry.pack(fill="x", pady=5)

ttk.Label(input_frame, text="Nozzle Expansion Ratio (ε):").pack(anchor="w")
eps_entry = ttk.Entry(input_frame); eps_entry.insert(0, "10"); eps_entry.pack(fill="x", pady=5)

run_btn = ttk.Button(input_frame, text="Run Analysis", command=run_simulation)
run_btn.pack(fill="x", pady=20)

res_label = ttk.Label(input_frame, text="Results will appear here", font=("Arial", 9, "italic"))
res_label.pack(pady=10)

# Plotting Panel
plot_frame = ttk.Frame(root, padding=10)
plot_frame.pack(side="right", fill="both", expand=True)

root.mainloop()