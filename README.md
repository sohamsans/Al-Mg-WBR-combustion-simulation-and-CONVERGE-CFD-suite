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
A critical design challenge in primary chamber engineering is maintaining metallic fuel expulsion efficiency:
* Excess combustion of metal particles in the primary combustor depletes the energy available for secondary reaction with ingested water, diminishing the overall ramjet benefit.
* Molten droplets coalesce and oxidize, forming refractory oxides ($\text{Al}_2\text{O}_3$ and $\text{MgO}$) that deposit on the combustor walls and nozzle throat as slag, degrading nozzle throat geometry and causing thermal damage.

Conversely, if the propellant does not ignite sufficiently, primary flame stability is lost. Magnesium provides a low ignition threshold (approximately $1100\text{ K}$) and fast reaction kinetics, whereas Aluminum provides high energy density but requires higher temperatures (approximately $2030\text{ K}$) to dissolve its protective oxide film.

This software models the physical trade-off between particle size distribution, chamber length, pressure, and metal composition, providing design boundaries that satisfy the research requirement of at least $90\%$ metal expulsion efficiency while minimizing slag deposition.

---

## 2. Physics and Mathematical Formulations

### 2.1 Thermochemical Equilibrium and Elemental Atom Balancing
NASA Chemical Equilibrium with Applications (CEA) calculates adiabatic flame temperature, gas molecular weight, and thermodynamic properties. To accurately represent blended composite propellants in RocketCEA, the mass percentages of HTPB binder, AP oxidizer, and Al-Mg alloys are converted into stoichiometric elemental mole numbers per $100\text{ g}$ of fuel mixture:

HTPB binder is modeled as $\text{C}_{7.07}\text{H}_{10.12}\text{O}_{0.20}$ with an average molecular weight of $98.32\text{ g/mol}$ and a heat of formation of $-12.5\text{ cal/g}$.
Ammonium perchlorate (AP) is modeled as $\text{NH}_4\text{ClO}_4$ with a molecular weight of $117.49\text{ g/mol}$.

For a fuel blend with HTPB mass fraction $w_{\text{HTPB}}$ and metal mass fraction $(100 - w_{\text{HTPB}})$:

$$\text{Moles of HTPB} = \frac{w_{\text{HTPB}}}{98.32 \text{ g/mol}}$$
$$n_{\text{C}} = 7.07 \times \text{Moles}_{\text{HTPB}}$$
$$n_{\text{H}} = 10.12 \times \text{Moles}_{\text{HTPB}}$$
$$n_{\text{O}} = 0.20 \times \text{Moles}_{\text{HTPB}}$$
$$n_{\text{Al}} = \frac{w_{\text{Al}}}{26.9815 \text{ g/mol}}$$
$$n_{\text{Mg}} = \frac{w_{\text{Mg}}}{24.305 \text{ g/mol}}$$
$$\text{Bulk Enthalpy } h_{\text{cal}} = \frac{w_{\text{HTPB}} \times (-12.5 \text{ cal/g})}{100.0}$$

The primary gas density is evaluated using the ideal gas equation of state:

$$\rho_g = \frac{P_c \cdot MW_g}{R_{\text{univ}} \cdot T_c}$$

Specific heat capacity at constant pressure is derived from the specific heat ratio $\gamma$:

$$c_{p,g} = \frac{\gamma R_{\text{univ}}}{(\gamma - 1.0) MW_g}$$

### 2.2 Al-Mg Alloy Properties and Ignition Kinetics
Alloy droplets exhibit composition-dependent density and specific heat:

$$\frac{1}{\rho_p} = \frac{w_{\text{Al}}}{\rho_{\text{Al}}} + \frac{1.0 - w_{\text{Al}}}{\rho_{\text{Mg}}}$$
$$c_{p,p} = w_{\text{Al}} c_{p,\text{Al}} + (1.0 - w_{\text{Al}}) c_{p,\text{Mg}}$$

where $\rho_{\text{Al}} = 2700\text{ kg/m}^3$, $\rho_{\text{Mg}} = 1738\text{ kg/m}^3$, $c_{p,\text{Al}} = 900\text{ J/(kg K)}$, and $c_{p,\text{Mg}} = 1020\text{ J/(kg K)}$.

Pure Aluminum requires reaching the oxide shell melting point ($2030\text{ K}$) for ignition, while Magnesium vaporizes and ignites at $1100\text{ K}$. In Al-Mg alloys, Magnesium vapor breakout disrupts the continuous $\text{Al}_2\text{O}_3$ layer. This non-linear transition is modeled as:

$$T_{\text{ign}}(w_{\text{Al}}) = T_{\text{ign,Mg}} + (T_{\text{ign,Al}} - T_{\text{ign,Mg}}) \cdot (w_{\text{Al}})^{1.5}$$

### 2.3 1D Lagrangian Droplet Dynamics and Drag Acceleration
The equation of motion for a spherical droplet moving through primary combustion gases is governed by standard aerodynamic drag:

$$\frac{du_p}{dt} = \frac{3}{4} \cdot \left( \frac{\rho_g C_d}{\rho_p d_p} \right) \cdot (u_g - u_p) \cdot |u_g - u_p|$$

The particle Reynolds number is defined by relative velocity:

$$Re_p = \frac{\rho_g \cdot |u_g - u_p| \cdot d_p}{\mu_g}$$

The drag coefficient $C_d$ is calculated using the Schiller-Naumann correlation:
* For $Re_p < 1000$: $C_d = \frac{24}{Re_p} \cdot \left(1.0 + 0.15 \cdot Re_p^{0.687}\right)$
* For $Re_p \ge 1000$: $C_d = 0.44$

Local gas velocity $u_g(x)$ accounts for thermal expansion along the combustor length $x$:

$$u_g(x) = u_{g,\text{inlet}} + (u_{g,\text{exit}} - u_{g,\text{inlet}}) \cdot \left(\frac{x}{L_c}\right)$$

### 2.4 Convective Heating and Pressure-Scaled d^1.8 Combustion Kinetics
Prior to ignition ($T_p < T_{\text{ign}}$), droplet heating is driven by forced convection:

$$Nu_p = 2.0 + 0.6 \cdot Re_p^{0.5} \cdot Pr_g^{1/3}$$
$$h_p = \frac{Nu_p \cdot k_g}{d_p}$$
$$\frac{dT_p}{dt} = \frac{h_p \cdot \pi d_p^2 \cdot (T_c - T_p)}{m_p \cdot c_{p,p}}$$

Once $T_p$ reaches $T_{\text{ign}}$, droplet combustion initiates. Metal droplet diameter decay follows a pressure-scaled power law ($d^{1.8}$ law):

$$\frac{d(d_p^{1.8})}{dt} = -K_b$$
$$K_b = K_{b0} \cdot \left(\frac{P_c}{P_{\text{ref}}}\right)^{0.27}$$
$$K_{b0} = 0.75 \times 10^{-6} \cdot w_{\text{Al}} + 1.85 \times 10^{-6} \cdot (1.0 - w_{\text{Al}}) \quad [\text{m}^{1.8}/\text{s}]$$

Instantaneous droplet mass is updated from diameter:

$$m_p = \frac{4}{3} \cdot \pi \cdot \left(\frac{d_p}{2}\right)^3 \cdot \rho_p$$

### 2.5 Expulsion Efficiency and Slag Deposition Formulations
Metal expulsion efficiency represents the fraction of metallic fuel that leaves the primary combustor without burning:

$$\eta_{\text{expulsion}}(x) = \text{clamp}\left( \frac{m_p(x)}{m_{p0}} \times 100\%, \, 0.0\%, \, 100.0\% \right)$$

Burned metal produces condensed oxide species ($\text{Al}_2\text{O}_3$ and $\text{MgO}$). The mass of oxide generated during time step $dt$ is:

$$m_{\text{oxide}} = (m_p(t) - m_p(t + dt)) \times 1.85$$

A fraction of this oxide deposits onto combustor surfaces based on droplet inertia and residence time:

$$f_{\text{dep}} = \text{clamp}\left( 0.10 \cdot \left(1.0 + 0.04 \cdot \frac{d_{p0,\mu\text{m}}}{20.0}\right), \, 0.0, \, 0.50 \right)$$
$$\text{Slag}_{\text{accumulated}} = \sum \left( m_{\text{oxide}} \cdot f_{\text{dep}} \right)$$

### 2.6 BATES Grain Geometry and Internal Ballistics
A Ballistic Test and Evaluation System (BATES) grain consists of $N_{\text{seg}}$ cylindrical segments with outer diameter $D_o$, initial inner port diameter $d_{i0}$, and segment length $L_{g0}$.

As combustion proceeds, the web distance burned is $y(t) = \int r_b \, dt$.
* Instantaneous inner diameter: $d_i(t) = d_{i0} + 2y(t)$
* Instantaneous segment length: $L_g(t) = L_{g0} - 2y(t)$ (for uninhibited ends)
* Burning surface area:

$$A_b(t) = N_{\text{seg}} \cdot \left[ \pi \cdot d_i(t) \cdot L_g(t) + \frac{\pi}{2} \cdot \left(D_o^2 - d_i(t)^2\right) \right]$$

* Klemmung ratio: $K_b(t) = \frac{A_b(t)}{A_t}$, where $A_t = \frac{\pi}{4} \cdot d_t^2$

Internal chamber pressure equilibrium follows Saint-Robert's burn rate law $r_b = a \cdot P_c^n$:

$$P_c(t) \, [\text{MPa}] = \left[ \frac{a \cdot \rho_{\text{prop}} \cdot c^* \cdot K_b(t)}{10^6} \right]^{\frac{1}{1 - n}}$$
$$r_b(t) \, [\text{m/s}] = a \cdot (P_c(t))^n$$

### 2.7 CONVERGE CFD Boundary Formulations
For 3D solid rocket motor combustion and regression CFD:
* Surface Mass Flux: $m''_{\text{wall}}(t) = \rho_{\text{prop}} \cdot r_b(t) \quad [\text{kg}/(\text{m}^2 \cdot \text{s})]$
* Gas Inflow Velocity: $v_{\text{inflow}}(t) = \frac{\rho_{\text{prop}}}{\rho_g(t)} \cdot r_b(t) \quad [\text{m/s}]$
* Regressing Boundary Displacement: $\frac{dy}{dt} = r_b(t)$

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
├── setup.py                       Python package configuration script
├── cea_wbr_final_lab.spec         PyInstaller executable build configuration
├── requirements.txt               Python dependency specifications
├── .gitignore                     Git tracking configuration
├── .gitattributes                 Git LFS tracking configuration
├── dist/
│   ├── Unified_WBR_CEA_Lab.exe    Standalone Windows executable (447 MB)
│   └── Unified_WBR_CEA_Lab_v1.0.0_win64.zip Archive build package
└── legacy/                        Archived early prototypes and preliminary scripts
```

---

## 4. Installation and Dependencies

### Prerequisites
* Python 3.10, 3.11, or 3.12 (64-bit)
* Microsoft C++ Build Tools (required by RocketCEA for Fortran/C wrappers)

### Installation Options

#### Standard pip install (from directory)
```bash
git clone https://github.com/sohamsans/Al-Mg-WBR-combustion-simulation-and-CONVERGE-CFD-suite.git
cd Al-Mg-WBR-combustion-simulation-and-CONVERGE-CFD-suite
pip install -e .
```

#### Installing from requirements file
```bash
pip install -r requirements.txt
```

---

## 5. Usage Instructions

### 5.1 Graphical User Interface
Launch the master simulation dashboard:
```bash
wbr-lab-gui
```
or via Python:
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
wbr-physics-test
```
or via Python:
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
Ran 5 tests in 4.358s

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
wbr-cfd-export
```
or via Python:
```bash
python bates_converge_cfd_exporter.py
```
Generated artifacts:
* `converge_bates_boundary.in`: Inflow and surface regression boundary setup.
* `converge_inflow_mass_flux.dat`: Transient wall mass flux and inflow injection velocity table.
* `converge_thermo.dat`: Gas molecular weight, gamma, and species mass fractions.
* `bates_internal_ballistics.csv`: Time-resolved grain geometry and chamber pressure dataset.

---

## 6. Standalone Executable Build

To rebuild the single-file executable using PyInstaller:
```bash
python -m PyInstaller cea_wbr_final_lab.spec --noconfirm
```
The compiled application will be generated in `dist/Unified_WBR_CEA_Lab.exe` and can be compressed into `dist/Unified_WBR_CEA_Lab_v1.0.0_win64.zip`.
