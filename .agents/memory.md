# Project Memory: CEA (Chemical Equilibrium with Applications) Water-Ramjet Lab

## Executive Summary
This repository implements thermodynamic simulation tools, parametric sweep engines, and GUI dashboards for **Al-Mg Metal Composites with HTPB Binder** reacting with **Ingested Water (H2O)** or **Ammonium Perchlorate (AP)** in **Water-Ramjet (WBR)** and solid rocket motor configurations. It relies on `rocketcea` (Python wrapper around NASA CEA) to analyze rocket performance ($I_{sp}$, $c^*$, $C_f$, Thrust), gas properties ($MW$, $\gamma$, $c_p$, $c_v$, $R$), and condensed phase slag fractions ($Al_2O_3$, $MgO$, unburned metal carryover) for CFD simulation input.

---

## File Map & Module Roles

| File | Type | Primary Purpose | Key Features / Dependencies |
| :--- | :--- | :--- | :--- |
| `test_physics_and_math.py` | Python Automated Test Suite | Unit & boundary tests for physics, numerical limits, and CONVERGE export | 5 suites, zero NaNs, zero Infs, mass/diameter monotonicity. |
| `MATH_AND_TEST_AUDIT.md` | Technical Documentation | Complete math audit, error catalog, equations, and test result report | Formulations, discovered bugs, fixes, and benchmark tables. |
| `cea_wbr_final_lab.spec` | PyInstaller Spec | Specification for building standalone Windows executable | Bundles RocketCEA, Fortran libs, Matplotlib, Tkinter into `dist/`. |
| `bates_converge_cfd_exporter.py` | Python Regression & CFD Exporter | BATES grain geometry regression ($A_b(t)$) & CONVERGE CFD boundary setup | Saint-Robert burn law, $P_c(t)$, wall mass flux $m''_{wall}(t)$, inflow velocity $v_{inflow}(t)$. |
| `cea_wbr_final_lab.py` | Tkinter GUI Suite | Unified Master Research Suite (5 Tabs) | NASA CEA equilibrium, 1D expulsion trajectory, parametric sweeps, BATES CFD setup. |
| `primary_combustion_model.py` | Python Physics Model | High-Fidelity 1D Lagrangian droplet dynamics & metal expulsion | Drag, non-linear alloy $T_{ign}$, $d^{1.8}$-law kinetics, gas acceleration, $\ge 90\%$ target. |
| `optimize_expulsion.py` | Python Sweep Engine | Multi-variable optimization for Al:Mg ratio, $d_{p0}$, $L_c$, and $P_c$ | Maps safe operating envelope contours for $\ge 90\%$ expulsion & low slag. |
| `newcea.py` | Tkinter GUI | DRDO baseline fuel sweep (15% HTPB + 85% Al-Mg) vs AP Oxidizer | Slag & unburned metal tracking, thrust computation, 3-plot dashboard. |
| `cea_gui_tool.py` | Tkinter GUI | Main monolithic dashboard (~44KB) | Advanced multi-sweep, fuel/ox custom definitions, Matplotlib embed, CSV export. |
| `unified_propellant_analysis.py` | Tkinter GUI | Al-Mg composite with H2O(L) oxidizer sweep | Robust slag fraction parsing, temperature/MW/Gamma plots. |
| `cea_almg_sweep.py` | CLI Script | Al (65%) + Mg (25%) + HTPB (10%) ternary blend vs H2O | Uses `rocketcea.blends.newFuelBlend`, standard sweep (MR 3-14 at 435.11 psia). |
| `prepare_propellant_data.py` | CLI Script | Post-processing script for CFD input generation | Computes $R_{specific}, c_p, c_v$, converts CSV to `Propellant_Data.dat`. |
| `wbr_research_lab.py` | Tkinter GUI | WBR Slag & Primary Expulsion Physics Lab | Al-Mg fuel + H2O/AP oxidizer, 1D expulsion trajectory plots & CFD export. |
| `WBR_CEA_gui.py` | Tkinter GUI | Defensive GUI variant | Fallback logic for `rocketcea` API variations. |
| `almg.py` | Script | Minimal Al-Mg + H2O ramjet demo | Quick check script for custom ox/fuel registration and basic plot. |
| `test_cea.py` | Script | Environment check | Standard LOX/LH2 baseline check. |
| `newcea.spec` | PyInstaller | Executable build spec | Packaging spec file for PyInstaller executable. |

---

## CONVERGE CFD Export Artifacts Generated
- `converge_bates_boundary.in` (CONVERGE boundary configuration file)
- `converge_inflow_mass_flux.dat` (Time-dependent wall mass flux $m''_{wall}$ & velocity $v_{inflow}$ table)
- `converge_thermo.dat` (CONVERGE thermodynamic properties & species mass fractions)
- `bates_internal_ballistics.csv` (Full grain regression time history dataset)

---

## Chemistry & Thermochemical Formulations

1. **Oxidizer Formulations:**
   - **H2O (Liquid):** `ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-68.32 t(k)=298.15` (or built-in `'H2O'`)
   - **AP (Ammonium Perchlorate):** `ox NH4CLO4(I) N 1.0 H 4.0 CL 1.0 O 4.0 wt%=100.0 h,cal=-70690.0 t(k)=298.15`

2. **Fuel Formulations:**
   - **Al-Mg Metals:** $Al$ ($h=0 \text{ cal/g}$), $Mg$ ($h=0 \text{ cal/g}$, $\rho=1.738 \text{ g/cc}$)
   - **HTPB Binder:** $C_{7.07} H_{10.12} O_{0.20}$ ($MW \approx 98.32 \text{ g/mol}$, $h = -12.5 \text{ cal/g}$)

3. **Slag / Condensed Phase Tracking:**
   - Target species: `AL2O3`, `MGO`, `AL2O3(L)`, `AL2O3(S)`, `MGO(S)`, `MG(S)`, `AL(L)`
   - Extracted via `cea.get_SpeciesMassFractions(Pc, MR, eps)`

4. **1D Droplet Kinetics & Expulsion Target:**
   - Drag force: $\frac{du_p}{dt} = \frac{3}{4} \frac{\rho_g C_d}{\rho_p d_p} (u_g - u_p) |u_g - u_p|$
   - Non-linear alloy ignition: $T_{ign}(w_{Al}) = T_{ign,Mg} + (T_{ign,Al} - T_{ign,Mg}) w_{Al}^{1.5}$
   - Droplet combustion: $\frac{d(d_p^{1.8})}{dt} = -K_b (P_c / P_{ref})^{0.27} Y_{O2}^{0.9}$
   - Target expulsion efficiency: $\eta_{expulsion} = \frac{m_{metal}(L_c)}{m_{metal, 0}} \ge 90\%$

---

## Verified Commands & Environment Notes
- **Python Dependencies:** `rocketcea`, `numpy`, `pandas`, `matplotlib`, `tkinter`
- **Execution Commands:**
  - `python test_physics_and_math.py` (Run full automated test suite)
  - `python bates_converge_cfd_exporter.py` (Run BATES grain regression & generate CONVERGE CFD files)
  - `python cea_wbr_final_lab.py` (Launch Unified Master Research & Simulation Suite)
  - `python primary_combustion_model.py` (Run 1D primary combustor trajectory physics benchmark)
  - `python optimize_expulsion.py` (Run multi-variable parametric optimization sweep)
  - `python -m PyInstaller cea_wbr_final_lab.spec --noconfirm` (Build standalone Windows executable)



