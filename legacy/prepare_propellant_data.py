import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Optional: Uncomment if you're using rocketcea and need temperatures
try:
    from rocketcea.cea_obj import CEA_Obj
    has_rocketcea = True
except ImportError:
    has_rocketcea = False

INPUT_CSV = "result.csv"
OUTPUT_CSV = "WBR_Research_Data_labeled.csv"
OUTPUT_DAT = "Propellant_Data.dat"
PLOTS_DIR = "plots"

os.makedirs(PLOTS_DIR, exist_ok=True)

# 1) Load CSV (no header expected) and assign column names
cols = ['MR', 'Isp_s', 'Pc_psia', 'P_out_psia', 'MW_gmol', 'Gamma', 'SlagFraction']
df = pd.read_csv(INPUT_CSV, header=None, names=cols)

# 2) Convert columns to appropriate types while ignoring errors
for col in ['MR', 'Isp_s', 'Pc_psia', 'P_out_psia', 'MW_gmol', 'Gamma', 'SlagFraction']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Check for NaNs and warn if any are found
if df.isnull().values.any():
    print("Warning: Non-numeric values were found and converted to NaN.")

# 3) Basic unit conversions and derived fields
df['MW_kg_per_mol'] = df['MW_gmol'] / 1000.0          # Convert MW to kg/mol
R_u = 8.31446261815324  # J/(mol K)
df['R_specific_J_per_kgK'] = R_u / df['MW_kg_per_mol']  # Specific gas constant

# Compute cp, cv from gamma
df['cp_J_per_kgK'] = df['Gamma'] * df['R_specific_J_per_kgK'] / (df['Gamma'] - 1.0)
df['cv_J_per_kgK'] = df['R_specific_J_per_kgK'] / (df['Gamma'] - 1.0)

# 4) Optional temperature fetch from rocketcea if installed
df['T_chamber_K'] = np.nan
df['T_exit_K'] = np.nan

if has_rocketcea:
    try:
        cea = CEA_Obj(oxName='H2O(L)', fuelName='AL')
        for idx, row in df.iterrows():
            mr = row['MR']
            temps = cea.get_Temperatures(Pc=row['Pc_psia'], MR=mr, eps=1.0)
            df.at[idx, 'T_chamber_K'] = temps[0]
            df.at[idx, 'T_exit_K'] = temps[2]
    except Exception as e:
        print("rocketcea attempt failed; proceeding without temps. Error:", e)

# 5) Save a labeled CSV
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved labeled CSV: {OUTPUT_CSV}")

# 6) Create Propellant_Data.dat
df[['MR', 'Isp_s', 'MW_gmol', 'Gamma', 'T_chamber_K', 'T_exit_K', 'SlagFraction', 'R_specific_J_per_kgK', 'cp_J_per_kgK', 'cv_J_per_kgK']].to_csv(OUTPUT_DAT, sep='\t', index=False)
print(f"Saved CFD-ready table: {OUTPUT_DAT}")

# 7) Plotting
plt.figure(figsize=(8,5))
plt.plot(df['MR'], df['Isp_s'], '-o')
plt.xlabel('MR')
plt.ylabel('Vacuum Isp (s)')
plt.title('MR vs Isp')
plt.grid(True)
plt.savefig(os.path.join(PLOTS_DIR, 'MR_vs_Isp.png'))
plt.close()

# Continue similarly for MW and Gamma plots...
print(f"Plots saved in {PLOTS_DIR}/")