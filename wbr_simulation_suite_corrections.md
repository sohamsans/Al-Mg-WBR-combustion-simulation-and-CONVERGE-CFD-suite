# WBR Simulation Suite: Code Corrections and Physical Fault Analysis

This document outlines critical logical, physical, and syntactical faults identified in the Unified Al-Mg Hydro-Reactive WBR Research Suite. The corrections are divided by their respective modules.

---

## 1. Thermochemical Setup & GUI (`cea_wbr_final_lab.py`)

### 1.1. Omission of HTPB Binder Mass in Fuel String
*   **The Fault:** The GUI collects an "HTPB Binder %" value, but the `run_cea_water_sweep` function never integrates it into the `AlMg_Fuel` string. The fuel is defined purely as `Al` and `Mg` summing to 100%.
*   **Reasoning:** Solid propellants require a polymeric binder to hold the metal powders together. Ignoring the 10%–15% binder mass drastically skews the calculated chamber temperature, molecular weight, and $\gamma$, resulting in physically invalid data for the CFD boundary conditions.
*   **Correction:** Normalize the metal mass fractions against the binder fraction.
    ```python
    htpb_pct = self._get_float(self.htpb_entry, "HTPB Binder %", 0.0, 100.0, 15.0)
    al_raw = self._get_float(self.al_entry, "Al Mass %", 0.0, 100.0, 75.0)
    mg_raw = 100.0 - al_raw
    
    # Calculate true mass fractions in the total mixture
    metal_fraction = (100.0 - htpb_pct) / 100.0
    al_true_pct = al_raw * metal_fraction
    mg_true_pct = mg_raw * metal_fraction
    ```

### 1.2. NASA CEA String Syntax Error (Mole vs. Mass Fractions)
*   **The Fault:** The fuel string is formatted as `Al {al_pct/100:.4f} Mg {mg_pct/100:.4f} wt%=100.0`.
*   **Reasoning:** In standard RocketCEA/NASA CEA string formatting, a decimal placed immediately after a chemical symbol denotes an *atomic or mole ratio*, not a mass fraction. The `wt%=100.0` at the end only defines the weight of the total mixture, failing to fix the internal mass distribution.
*   **Correction:** Use the `wt%=` identifier directly for each constituent to force mass-based evaluation. Also, assign a realistic heat of formation for HTPB (e.g., $-12.5 \text{ kcal/mol}$).
    ```python
    fuel_str = (f"fuel AlMg_Fuel Al 1 wt%={al_true_pct:.2f} "
                f"Mg 1 wt%={mg_true_pct:.2f} "
                f"HTPB 1 wt%={htpb_pct:.2f} h,cal=-12.5 t(k)=298.15")
    add_new_fuel('AlMg_Fuel', fuel_str)
    ```

### 1.3. Vacuum vs. Optimum Specific Impulse Mislabeling
*   **The Fault:** The script calls `cea.get_Isp(...)` but labels and plots the result as `Isp_vac_s`.
*   **Reasoning:** The `get_Isp` method returns the *optimum* specific impulse assuming the nozzle perfectly expands the gas to match the ambient atmospheric pressure ($p_e = p_a$). True vacuum specific impulse ($I_{vac}$) must account for the additional pressure thrust term in a vacuum environment ($p_a = 0$).
*   **Correction:** Switch to the native vacuum function.
    ```python
    # Change this:
    isp_vac = cea.get_Isp(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
    # To this:
    isp_vac = cea.get_Ivac(Pc=pc, MR=mr, eps=eps, frozen=sub_frozen)
    ```

---

## 2. Primary Combustor Physics (`primary_combustion_model.py`)

### 2.1. Inverted Expulsion/Combustion Efficiency Logic
*   **The Fault:** `expulsion_eff = float(np.clip((m_p / m_p0) * 100.0, 0.0, 100.0))`. The optimization sweep aims for this metric to be $\ge 90\%$.
*   **Reasoning:** `m_p` is the *current* mass of the droplet. `m_p / m_p0` calculates the percentage of the original droplet mass that remains **unburned**. By targeting $\ge 90\%$, the code actively optimizes for a combustor that fails to burn the metal, ejecting it raw.
*   **Correction:** Invert the logic to track *combustion* efficiency (how much mass was converted to gas/slag).
    ```python
    # Calculate how much of the particle has successfully burned
    combustion_eff_pct = float(np.clip((1.0 - (m_p / m_p0)) * 100.0, 0.0, 100.0))
    ```
    *(Note: You will also need to update the variable names and dictionary keys from `expulsion_eff_pct` to `combustion_eff_pct` throughout both files).*

### 2.2. Hardcoded Stoichiometric Oxide Ratio
*   **The Fault:** Slag mass is calculated using a static multiplier: `oxide_mass = dm_dt * dt * 1.85`.
*   **Reasoning:** The molecular weight ratio of Aluminum to Alumina ($Al \to Al_2O_3$) is $\approx 1.889$. The ratio of Magnesium to Magnesia ($Mg \to MgO$) is $\approx 1.658$. Using a hardcoded factor of $1.85$ violates mass conservation whenever the parametric sweep changes the $Al:Mg$ ratio.
*   **Correction:** Dynamically weight the oxide multiplier based on the current alloy fraction.
    ```python
    # Dynamic oxidation mass multiplier
    oxide_ratio = al_frac * 1.889 + (1.0 - al_frac) * 1.658
    oxide_mass = dm_dt * dt * oxide_ratio
    ```

### 2.3. Kinematic Gas Velocity Profile Violates Mass Conservation
*   **The Fault:** The local gas velocity is linearly interpolated: `u_g_local = u_g_inlet + (u_g_exit - u_g_inlet) * (x / L)`. Furthermore, `u_g_exit` is hardcoded to `u_g_inlet * 1.8`.
*   **Reasoning:** In a solid rocket combustor, the gas velocity accelerates dynamically due to mass addition from the regressing grain wall along the length of the chamber ($x$). A hardcoded $1.8$ expansion factor decouples the particle drag physics from the actual internal ballistics and local surface mass flux.
*   **Correction:** Calculate velocity strictly using mass conservation ($\dot{m} = \rho \cdot u \cdot A$).
    ```python
    # Inside the integration loop:
    # Calculate total mass flow rate at position x
    # Assuming mdot_inlet is the head-end injection, plus the integral of wall mass flux
    mdot_local = mass_flow_rate_kg_s + (mass_flux_wall_kg_m2s * np.pi * chamber_diameter_m * x)
    u_g_local = mdot_local / (gas['rho_g'] * area_combustor)
    ```

---

## 3. CFD Boundary Export (`bates_converge_cfd_exporter.py`)

### 3.1. Hardcoded CFD Boundary Configurations
*   **The Fault:** `f.write(f"PROPELLANT_DENSITY: {1750.0} kg/m^3\n")` and the entire `SPECIES_MASS_FRACTIONS` block are hardcoded strings.
*   **Reasoning:** The module is designed to export dynamic simulation parameters, but it forces a static density of $1750.0 \text{ kg/m}^3$ regardless of what is passed via `rho_prop`. Furthermore, hardcoding the exhaust species (like $CO$, $CO_2$, $Al_2O_3$) completely overrides and ignores the accurate thermochemical equilibrium calculations performed by RocketCEA, making the parametric sweeps useless for the CFD export.
*   **Correction:** Use the dynamic `rho_prop` variable and extract species fractions directly from the CEA object.
    ```python
    # For Density:
    f.write(f"PROPELLANT_DENSITY: {rho_prop} kg/m^3\n")
    
    # For Species (requires passing the species dictionary from CEA output):
    # (Assuming 'species_dict' is retrieved from RocketCEA and passed to this function)
    f.write("SPECIES_MASS_FRACTIONS:\n")
    for species, fraction in species_dict.items():
        if fraction > 1e-4:  # Filter out trace elements
            f.write(f"  {species:<5} : {fraction:.4f}\n")
    ```