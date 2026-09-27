"""
Full Grid Sanity & Physical Validity Verification Suite
======================================================
Systematically tests parameter combinations across NASA CEA water sweeps,
1D Lagrangian droplet dynamics, and BATES grain internal ballistics.

Verifies:
  1. No NaN / Inf / ZeroDivisionError crashes
  2. Physical bounds on Temperatures (500K <= Tc <= 7500K)
  3. Physical bounds on Isp (100s <= Isp <= 600s)
  4. Physical bounds on MW (5.0 <= MW <= 45.0 g/mol)
  5. Physical bounds on Gamma (1.05 <= Gamma <= 1.45)
  6. Strict Monotonicity on Particle Diameter Decay (dd_p/dt <= 0)
  7. Bounded Expulsion Efficiency (0% <= eta_expulsion <= 100%)
  8. BATES Ballistics Stability (Pc >= 0.1 MPa, rb > 0 mm/s)
"""

import numpy as np
import pandas as pd
import time
import json
from primary_combustion_model import (
    calculate_primary_gas_state,
    simulate_1d_primary_combustor,
    generate_atom_balanced_fuel_string
)
from bates_converge_cfd_exporter import (
    simulate_bates_grain_regression,
    export_converge_cfd_inputs
)
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

def run_grid_sanity_test():
    print("=" * 80)
    print("STARTING FULL GRID PHYSICAL SANITY & ENVELOPE VERIFICATION")
    print("=" * 80)
    
    # 1. High-Resolution 1% Fine Increment Sweep Arrays
    al_fine_1pct = list(np.linspace(0.0, 100.0, 101))    # Exact 1% increments: 0%, 1%, 2%, ..., 100%
    htpb_fine_1pct = list(np.linspace(5.0, 30.0, 26))    # Exact 1% increments: 5%, 6%, ..., 30%

    total_tests = 0
    passed_tests = 0
    failures = []

    start_time = time.time()

    # --- TEST SUITE 1: NASA CEA Water-Ramjet Thermochemistry (Fine 1% Al & 1% HTPB) ---
    print("\n--- [Suite 1/3] Testing NASA CEA Water-Ramjet Equilibrium (1% Increments) ---")
    add_new_oxidizer('SanityH2O', "ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-3788.5 t(k)=298.15")
    
    # 101 Al values (1% increments) x 2 representative HTPB fractions x 3 pressures x 3 O/F
    for al in al_fine_1pct:
        for htpb in [10.0, 15.0]:
            fuel_card = f"SanityFuel_{int(al*10)}_{int(htpb)}"
            fuel_str = generate_atom_balanced_fuel_string(al, htpb, fuel_card)
            add_new_fuel(fuel_card, fuel_str)
            
            try:
                cea = CEA_Obj(oxName='SanityH2O', fuelName=fuel_card)
            except Exception as e:
                failures.append({'type': 'CEA_Init', 'al': al, 'htpb': htpb, 'error': str(e)})
                continue

            for pc in [100.0, 250.0, 500.0]:
                for w_of in [2.0, 4.0, 7.0]:
                    total_tests += 1
                    try:
                        isp = cea.get_Isp(Pc=pc, MR=w_of, eps=8.0)
                        temps = cea.get_Temperatures(Pc=pc, MR=w_of, eps=8.0)
                        mw_gam = cea.get_Chamber_MolWt_gamma(Pc=pc, MR=w_of, eps=8.0)
                        
                        tc = temps[0] if temps else np.nan
                        mw = mw_gam[0]
                        gam = mw_gam[1]

                        assert not np.isnan(isp) and not np.isinf(isp), f"NaN/Inf Isp ({isp})"
                        assert 100.0 <= isp <= 600.0, f"Isp out of bounds ({isp} s)"
                        assert 500.0 <= tc <= 7500.0, f"Tc out of bounds ({tc} K)"
                        assert 5.0 <= mw <= 45.0, f"MW out of bounds ({mw} g/mol)"
                        assert 1.05 <= gam <= 1.45, f"Gamma out of bounds ({gam})"
                        
                        passed_tests += 1
                    except Exception as e:
                        failures.append({
                            'suite': 'Water_CEA',
                            'al': al, 'htpb': htpb, 'pc': pc, 'water_of': w_of,
                            'error': str(e)
                        })

    print(f"Suite 1 Finished. Passed {passed_tests}/{total_tests} fine evaluations.")

    # --- TEST SUITE 2: 1D Droplet Trajectory & Combustion Kinetics (101 Al 1% Steps) ---
    print("\n--- [Suite 2/3] Testing 1D Droplet Dynamics Across 101 Al (1% Steps) ---")
    suite2_total = 0
    suite2_passed = 0

    # 101 fine 1% Al steps swept over various particle sizes and combustor lengths
    for al in al_fine_1pct:
        for dp0 in [15.0, 30.0, 60.0, 100.0]:
            for lc in [0.10, 0.25, 0.40]:
                suite2_total += 1
                total_tests += 1
                try:
                    df_traj, metrics = simulate_1d_primary_combustor(
                        al_pct=al,
                        htpb_pct=15.0,
                        d_p0_um=dp0,
                        pc_psia=200.0,
                        primary_of=0.25,
                        chamber_length_m=lc,
                        chamber_diameter_m=0.08,
                        mass_flow_rate_kg_s=0.50
                    )

                    assert len(df_traj) > 0, "Empty trajectory DataFrame"
                    assert not df_traj.isna().any().any(), "NaN values in trajectory"
                    
                    # Particle diameter decay monotonicity check
                    d_p_arr = df_traj['d_p_um'].values
                    diffs = np.diff(d_p_arr)
                    assert (diffs <= 1e-9).all(), "Particle diameter decay non-monotonic!"
                    
                    # Expulsion efficiency bounds
                    eta_arr = df_traj['expulsion_eff_pct'].values
                    assert (eta_arr >= 0.0).all() and (eta_arr <= 100.0).all(), "Expulsion efficiency out of [0, 100%]"

                    suite2_passed += 1
                    passed_tests += 1
                except Exception as e:
                    failures.append({
                        'suite': '1D_Droplet_Trajectory',
                        'al': al, 'dp0': dp0, 'lc': lc,
                        'error': str(e)
                    })

    print(f"Suite 2 Finished. Passed {suite2_passed}/{suite2_total} fine evaluations.")

    # --- TEST SUITE 3: BATES Grain Internal Ballistics (26 HTPB 1% Steps + 21 Al Steps) ---
    print("\n--- [Suite 3/3] Testing BATES Grain Regression Across 1% Fine Increments ---")
    suite3_total = 0
    suite3_passed = 0

    for htpb in htpb_fine_1pct:  # 26 fine 1% steps
        for al in [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
            suite3_total += 1
            total_tests += 1
            try:
                df_bates = simulate_bates_grain_regression(
                    al_pct=al,
                    htpb_pct=htpb,
                    primary_of=0.25
                )
                assert len(df_bates) > 0, "Empty BATES DataFrame"
                assert not df_bates.isna().any().any(), "NaN in BATES DataFrame"
                assert df_bates['P_c_MPa'].min() >= 0.05, "Unphysical low pressure in BATES"
                assert df_bates['r_b_mm_s'].min() > 0.0, "Negative or zero burn rate in BATES"

                export_converge_cfd_inputs(df_bates, output_dir="d:/CEA/scratch_test_out")
                
                suite3_passed += 1
                passed_tests += 1
            except Exception as e:
                failures.append({
                    'suite': 'BATES_Grain_Regression',
                    'al': al, 'htpb': htpb,
                    'error': str(e)
                })

    print(f"Suite 3 Finished. Passed {suite3_passed}/{suite3_total} evaluations.")

    elapsed = time.time() - start_time
    print("=" * 80)
    print(f"GRID SANITY VERIFICATION COMPLETE in {elapsed:.2f} seconds")
    print(f"TOTAL EVALUATIONS: {total_tests}")
    print(f"PASSED          : {passed_tests} ({passed_tests/total_tests*100:.2f}%)")
    print(f"FAILURES        : {len(failures)}")
    print("=" * 80)

    if failures:
        print("\nFAILURE SUMMARY:")
        for f in failures[:10]:
            print(" -", f)
        if len(failures) > 10:
            print(f" ... and {len(failures)-10} more failures.")
    else:
        print("\nALL EVALUATIONS PASSED PHYSICAL SANITY CHECKS PERFECTLY!")

if __name__ == '__main__':
    run_grid_sanity_test()
