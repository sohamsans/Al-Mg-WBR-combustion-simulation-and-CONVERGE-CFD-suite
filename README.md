# Hydro-Reactive Combustion and BATES Grain Internal Ballistics Simulation Suite

An integrated physics-based simulation framework and analytical toolset for modeling fuel-rich primary combustion, metallic fuel expulsion efficiency, oxide slag accumulation, and BATES grain internal ballistics for water-breathing ramjet (WBR) propulsion systems. The suite includes complete boundary condition generators for 3D CONVERGE CFD simulations.

## Table of Contents
1. Introduction and Propulsion Theory
2. Physics and Mathematical Formulations
   - Thermochemical Equilibrium and Elemental Atom Balancing
   - Al-Mg Alloy Properties and Ignition Kinetics
   - 1D Lagrangian Droplet Dynamics and Drag Acceleration
   - Convective Heating and Pressure-Scaled d^1.8 Combustion Kinetics
   - Expulsion Efficiency and Slag Deposition Formulations
   - BATES Grain Geometry and Internal Ballistics
   - CONVERGE CFD Boundary Formulations
3. Repository Architecture
4. Installation and Dependencies
5. Usage Instructions
   - Graphical User Interface
   - Automated Test Suite
   - Command-Line Simulation and Optimization
   - CONVERGE CFD Export
6. Standalone Executable Build

---

## 1. Introduction and Propulsion Theory

Water-breathing ramjet (WBR) engines utilize surrounding liquid water as an oxidizer, providing significantly higher theoretical specific impulse than conventional closed-system underwater rockets. Operating a WBR involves two sequential combustion stages:

1. Primary Combustor (Fuel-Rich Gas Generator):
   A composite solid propellant grain containing metallic fuel (Aluminum and Magnesium powders) embedded within a hydroxyl-terminated polybutadiene (HTPB) binder and ammonium perchlorate (AP) oxidizer burns under fuel-rich conditions. This stage heats the propellant, gasifies the binder, and expels hot, partially reacted metallic droplets and fuel gases into the combustor duct.

2. Secondary Combustor (Hydro-Reactive Stage):
   Downstream of the primary nozzle, the expelled molten metal droplets react vigorously with ingested ambient water, generating superheated steam, hydrogen gas, and substantial secondary thrust.

### The Metal Combustion Dilemma: Expulsion vs. Slag Accumulation
A critical design challenge in primary chamber engineering is maintaining metallic fuel expulsion efficiency. If the metal particles burn excessively in the primary combustor:
- The energy available for secondary reaction with ingested water is depleted, diminishing the ramjet benefit.
- Molten droplets coalesce and oxidize, forming refractory oxides (Al2O3 and MgO) that deposit on the combustor walls and nozzle throat as slag, degrading nozzle throat geometry and causing thermal damage.

Conversely, if the propellant does not ignite sufficiently, primary flame stability is lost. Magnesium provides a low ignition threshold (approximately 1100 K) and fast reaction kinetics, whereas Aluminum provides high energy density but requires higher temperatures (approximately 2030 K) to dissolve its protective oxide film.

This software models the physical trade-off between particle size distribution, chamber length, pressure, and metal composition, providing design boundaries that satisfy the research requirement of at least 90% metal expulsion efficiency while minimizing slag deposition.

---

## 2. Physics and Mathematical Formulations

### 2.1 Thermochemical Equilibrium and Elemental Atom Balancing
NASA Chemical Equilibrium with Applications (CEA) calculates adiabatic flame temperature, gas molecular weight, and thermodynamic properties. To accurately represent blended composite propellants in RocketCEA, the mass percentages of HTPB binder, AP oxidizer, and Al-Mg alloys are converted into stoichiometric elemental mole numbers per 100 grams of fuel mixture:

HTPB binder is modeled as C7.07 H10.12 O0.20 with an average molecular weight of 98.32 g/mol and a heat of formation of -12.5 cal/g.
Ammonium perchlorate (AP) is modeled as NH4ClO4 with a molecular weight of 117.49 g/mol.

For a fuel blend with HTPB mass fraction w_HTPB and metal mass fraction (100 - w_HTPB):
- Moles of HTPB = w_HTPB / 98.32
- n_C = 7.07 * Moles_HTPB
- n_H = 10.12 * Moles_HTPB
- n_O = 0.20 * Moles_HTPB
- n_Al = w_Al / 26.9815
- n_Mg = w_Mg / 24.305
- Bulk Enthalpy h_cal = (w_HTPB * -12.5) / 100.0 [cal/g]

The primary gas density is evaluated using the ideal gas equation of state:
rho_g = (P_c * MW_g) / (R_univ * T_c)

Specific heat capacity at constant pressure is derived from the specific heat ratio gamma:
c_p = (gamma * R_univ) / ((gamma - 1.0) * MW_g)

### 2.2 Al-Mg Alloy Properties and Ignition Kinetics
Alloy droplets exhibit composition-dependent density and specific heat:
1 / rho_p = (w_Al / rho_Al) + ((1.0 - w_Al) / rho_Mg)
c_p,p = w_Al * c_p,Al + (1.0 - w_Al) * c_p,Mg

where rho_Al = 2700 kg/m^3, rho_Mg = 1738 kg/m^3, c_p,Al = 900 J/(kg K), and c_p,Mg = 1020 J/(kg K).

Pure Aluminum requires reaching the oxide shell melting point (2030 K) for ignition, while Magnesium vaporizes and ignites at 1100 K. In Al-Mg alloys, Magnesium vapor breakout disrupts the continuous Al2O3 layer. This non-linear transition is modeled as:
T_ign(w_Al) = T_ign,Mg + (T_ign,Al - T_ign,Mg) * (w_Al)^1.5

### 2.3 1D Lagrangian Droplet Dynamics and Drag Acceleration
The equation of motion for a spherical droplet moving through primary combustion gases is governed by standard aerodynamic drag:
d(u_p) / dt = (3 / 4) * (rho_g * C_d / (rho_p * d_p)) * (u_g - u_p) * |u_g - u_p|

The particle Reynolds number is defined by relative velocity:
Re_p = (rho_g * |u_g - u_p| * d_p) / mu_g

The drag coefficient C_d is calculated using the Schiller-Naumann correlation:
- For Re_p < 1000: C_d = (24 / Re_p) * (1.0 + 0.15 * Re_p^0.687)
- For Re_p >= 1000: C_d = 0.44

Local gas velocity u_g(x) accounts for thermal expansion along the combustor length x:
u_g(x) = u_g,inlet + (u_g,exit - u_g,inlet) * (x / L_c)

### 2.4 Convective Heating and Pressure-Scaled d^1.8 Combustion Kinetics
Prior to ignition (T_p < T_ign), droplet heating is driven by forced convection:
Nu_p = 2.0 + 0.6 * (Re_p^0.5) * (Pr_g^(1/3))
h_p = (Nu_p * k_g) / d_p
d(T_p) / dt = (h_p * pi * d_p^2 * (T_c - T_p)) / (m_p * c_p,p)

Once T_p reaches T_ign, droplet combustion initiates. Metal droplet diameter decay follows a pressure-scaled power law (d^1.8 law):
d(d_p^1.8) / dt = -K_b
K_b = K_b0 * (P_c / P_ref)^0.27
K_b0 = 0.75e-6 * w_Al + 1.85e-6 * (1.0 - w_Al) [m^1.8 / s]

Instantaneous droplet mass is updated from diameter:
m_p = (4 / 3) * pi * (d_p / 2)^3 * rho_p

### 2.5 Expulsion Efficiency and Slag Deposition Formulations
Metal expulsion efficiency represents the fraction of metallic fuel that leaves the primary combustor without burning:
eta_expulsion(x) = (m_p(x) / m_p0) * 100 [%]

Burned metal produces condensed oxide species (Al2O3 and MgO). The mass of oxide generated during time step dt is:
m_oxide = (m_p(t) - m_p(t + dt)) * 1.85

A fraction of this oxide deposits onto combustor surfaces based on droplet inertia and residence time:
f_dep = clamp(0.10 * (1.0 + 0.04 * (d_p0_um / 20.0)), 0.0, 0.50)
Slag_accumulated = sum(m_oxide * f_dep)

### 2.6 BATES Grain Geometry and Internal Ballistics
A Ballistic Test and Evaluation System (BATES) grain consists of N_seg cylindrical segments with outer diameter D_o, initial inner port diameter d_i0, and segment length L_g0.

As combustion proceeds, the web distance burned is y(t) = integral(r_b dt).
- Instantaneous inner diameter: d_i(t) = d_i0 + 2 * y(t)
- Instantaneous segment length: L_g(t) = L_g0 - 2 * y(t) (for uninhibited ends)
- Burning surface area:
  A_b(t) = N_seg * [ pi * d_i(t) * L_g(t) + 0.5 * pi * (D_o^2 - d_i(t)^2) ]
- Klemmung ratio: K_b(t) = A_b(t) / A_t, where A_t = (pi / 4) * d_t^2

Internal chamber pressure equilibrium follows Saint-Robert's burn rate law r_b = a * P_c^n:
P_c(t) [MPa] = [ (a * rho_prop * c* * K_b(t)) / 1.0e6 ]^(1 / (1 - n))
r_b(t) [m/s] = a * (P_c(t))^n

### 2.7 CONVERGE CFD Boundary Formulations
For 3D solid rocket motor combustion and regression CFD:
- Surface Mass Flux: m''_wall(t) = rho_prop * r_b(t) [kg / (m^2 * s)]
- Gas Inflow Velocity: v_inflow(t) = (rho_prop / rho_g(t)) * r_b(t) [m / s]
- Regressing Boundary Displacement: dy / dt = r_b(t)

---

## 3. Repository Architecture

```text
d:/CEA/
├── cea_wbr_final_lab.py           Master graphical application (Tkinter, 5 specialized tabs)
├── primary_combustion_model.py    1D Lagrangian droplet dynamics and expulsion physics core
├── bates_converge_cfd_exporter.py BATES grain regression and CONVERGE CFD boundary generator
├── optimize_expulsion.py          Parametric sweep and design envelope optimization engine
├── test_physics_and_math.py       Automated unit and validation test suite
├── MATH_AND_TEST_AUDIT.md         Technical mathematical audit and error catalog
├── cea_wbr_final_lab.spec         PyInstaller executable build configuration
├── requirements.txt               Python dependency specifications
├── .gitignore                     Git tracking configuration
├── dist/
│   └── Unified_WBR_CEA_Lab.exe    Standalone Windows executable (447 MB)
└── legacy/                        Archived early prototypes and preliminary scripts
```

---

## 4. Installation and Dependencies

### Prerequisites
- Python 3.10, 3.11, or 3.12 (64-bit)
- Microsoft C++ Build Tools (required by RocketCEA for Fortran/C wrappers)

### Installation
Clone the repository and install the required dependencies:
```bash
git clone <repository-url>
cd <repository-directory>
pip install -r requirements.txt
```

---

## 5. Usage Instructions

### 5.1 Graphical User Interface
Launch the master simulation dashboard:
```bash
python cea_wbr_final_lab.py
```
The interface contains five operational tabs:
1. NASA CEA Equilibrium and Water Sweep: Evaluates specific impulse, combustion temperature, molecular weight, and gamma across water-to-fuel ratios.
2. 1D Droplet Trajectory and Expulsion Model: Tracks droplet diameter decay, drag velocity, gas acceleration, and expulsion efficiency along the combustor.
3. Parametric Sweeps and Optimization: Executes multi-variable sweeps over particle diameter, combustor length, and pressure to identify configurations where expulsion efficiency exceeds 90%.
4. CFD Exporter and Data Viewer: Previews numerical tables and formats CFD-ready data.
5. BATES Grain and CONVERGE CFD Setup: Simulates transient internal ballistics and exports CONVERGE CFD boundary files.

### 5.2 Automated Test Suite
Run the verification suite to check physics limits, boundary conditions, and mathematical sanity:
```bash
python test_physics_and_math.py
```
Expected output:
```text
test_1d_droplet_combustion_bounds ... ok
test_alloy_properties_bounds ... ok
test_bates_grain_regression_bounds ... ok
test_converge_cfd_export ... ok
test_primary_gas_state_robustness ... ok

----------------------------------------------------------------------
Ran 5 tests in 4.170s

OK
```

### 5.3 Command-Line Simulation and Optimization
Execute standalone 1D droplet trajectory physics:
```bash
python primary_combustion_model.py
```

Execute multi-variable parametric optimization sweep:
```bash
python optimize_expulsion.py
```

### 5.4 CONVERGE CFD Export
Generate input boundary condition files for CONVERGE CFD:
```bash
python bates_converge_cfd_exporter.py
```
Generated artifacts:
- `converge_bates_boundary.in`: Inflow and surface regression boundary setup.
- `converge_inflow_mass_flux.dat`: Transient wall mass flux and inflow injection velocity table.
- `converge_thermo.dat`: Gas molecular weight, gamma, and species mass fractions.
- `bates_internal_ballistics.csv`: Time-resolved grain geometry and chamber pressure dataset.

---

## 6. Standalone Executable Build

To rebuild the single-file executable using PyInstaller:
```bash
python -m PyInstaller cea_wbr_final_lab.spec --noconfirm
```
The compiled application will be generated in `dist/Unified_WBR_CEA_Lab.exe`.
