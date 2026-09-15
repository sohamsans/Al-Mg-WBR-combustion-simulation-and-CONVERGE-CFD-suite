import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

# 1. DEFINE WATER MANUALLY (Oxidizer)
# H2O(L) at 298.15K. Heat of formation is -285.83 kJ/mol -> -3788.5 cal/g
ox_card = "ox H2O(L) H 2 O 1  wt%=100.0  h,cal=-3788.5  t(k)=298.15"
add_new_oxidizer('MyWater', ox_card)

# 2. DEFINE AL-MG COMPOSITE (Fuel)
# Al: h=0, Mg: h=0 (pure elements)
fuel_card = "fuel AlMg_Fuel  Al 0.7  Mg 0.3  wt%=100.0  h,cal=0.0  t(k)=298.15"
add_new_fuel('AlMg_Fuel', fuel_card)

# 3. INITIALIZE
try:
    isp_obj = CEA_Obj(oxName='MyWater', fuelName='AlMg_Fuel')
    
    # Range of Water-to-Fuel ratios
    mr_range = np.linspace(1.5, 6.0, 40)
    isp_vals = []
    temp_vals = []

    print("Running Al-Mg + Water Ramjet Simulation...")
    for mr in mr_range:
        # 150 psia is a common combustor pressure for ramjets
        isp_vals.append(isp_obj.get_Isp(Pc=150.0, MR=mr, eps=10.0))
        temp_vals.append(isp_obj.get_Temperatures(Pc=150.0, MR=mr, eps=10.0)[0])

    # 4. PLOT
    fig, ax1 = plt.subplots(figsize=(9, 5))
    ax1.plot(mr_range, isp_vals, 'b-', label='Isp')
    ax1.set_ylabel('Vacuum Isp (s)', color='b')
    ax1.set_xlabel('Water-to-Fuel Ratio')
    
    ax2 = ax1.twinx()
    ax2.plot(mr_range, temp_vals, 'r--', label='Temp')
    ax2.set_ylabel('Combustion Temp (K)', color='r')
    
    plt.title('Propulsion Analysis: Al-Mg Composite with Ingested Water')
    plt.grid(True, alpha=0.3)
    plt.show()

except Exception as e:
    print(f"Error: {e}")