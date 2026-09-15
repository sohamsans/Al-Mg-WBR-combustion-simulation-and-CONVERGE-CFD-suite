from rocketcea.cea_obj import CEA_Obj

# Define the engine object: (Oxidizer, Fuel)
# RocketCEA understands standard names like LOX, LH2, RP1, NTO, etc.
isp_obj = CEA_Obj(oxName='LOX', fuelName='LH2')

# Calculate Vacuum Isp
# Pc = Pressure (psia), MR = Mixture Ratio (O/F)
# eps = Area Ratio (Exit Area / Throat Area)
isp_vac = isp_obj.get_Isp(Pc=500.0, MR=6.0, eps=40.0)

print(f"\n--- Simulation Results ---")
print(f"Vacuum Isp: {isp_vac:.2f} seconds")
print(f"Chamber Temp: {isp_obj.get_Temperatures(Pc=500.0, MR=6.0, eps=40.0)[0]:.2f} K")