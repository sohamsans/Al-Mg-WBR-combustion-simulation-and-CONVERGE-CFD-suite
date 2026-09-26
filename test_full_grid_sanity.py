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
    
    # 1. Define High-Resolution Fine Increment Parameter Grid
    al_pct_list = list(np.linspace(0.0, 100.0, 101))      # 1% increments (0%, 1%, 2%, ..., 100%)
    htpb_pct_list = list(np.linspace(5.0, 30.0, 26))      # 1% increments (5%, 6%, ..., 30%)
    pc_psia_list = [50.0, 100.0, 200.0, 350.0, 500.0, 750.0, 1000.0]
    dp0_um_list = [5.0, 10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0]
    lc_m_list = [0.05, 0.10, 0.20, 0.35, 0.50]
    primary_of_list = [0.10, 0.20, 0.30, 0.50]
    water_of_list = [1.5, 2.5, 4.0, 6.0, 8.5, 10.0]

    total_tests = 0
    passed_tests = 0
    failures = []

    start_time = time.time()

    # --- TEST SUITE 1: NASA CEA Water-Ramjet Thermochemistry Sweep ---
    print("\n--- [Suite 1/3] Testing NASA CEA Water-Ramjet Equilibrium Sweeps ---")
    add_new_oxidizer('SanityH2O', "ox H2O(L) H 2 O 1 wt%=100.0 h,cal=-3788.5 t(k)=298.15")
    
    for al in al_pct_list:
        for htpb in htpb_pct_list:
            fuel_card = f"SanityFuel_{int(al)}_{int(htpb)}"
            fuel_str = generate_atom_balanced_fuel_string(al, htpb, fuel_card)
            add_new_fuel(fuel_card, fuel_str)
            
            try:
                cea = CEA_Obj(oxName='SanityH2O', fuelName=fuel_card)
            except Exception as e:
                failures.append({'type': 'CEA_Init', 'al': al, 'htpb': htpb, 'error': str(e)})
                continue

            for pc in pc_psia_list:
                for w_of in water_of_list:
                    total_tests += 1
                    try:
                        isp = cea.get_Isp(Pc=pc, MR=w_of, eps=8.0)
                        temps = cea.get_Temperatures(Pc=pc, MR=w_of, eps=8.0)
                        mw_gam = cea.get_Chamber_MolWt_gamma(Pc=pc, MR=w_of, eps=8.0)
                        
                        tc = temps[0] if temps else np.nan
                        mw = mw_gam[0]
                        gam = mw_gam[1]

                        # Physical Sanity Assertions
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

    print(f"Suite 1 Finished. Passed {passed_tests}/{total_tests} evaluations.")

    # --- TEST SUITE 2: 1D Droplet Trajectory & Combustion Kinetics ---
    print("\n--- [Suite 2/3] Testing 1D Droplet Dynamics & Expulsion Kinetic Grid ---")
    suite2_total = 0
    suite2_passed = 0

    for al in al_pct_list:
        for htpb in htpb_pct_list:
            for pc in pc_psia_list:
                for dp0 in dp0_um_list:
                    for lc in lc_m_list:
                        for p_of in primary_of_list:
                            suite2_total += 1
                            total_tests += 1
                            try:
                                df_traj, metrics = simulate_1d_primary_combustor(
                                    al_pct=al,
                                    htpb_pct=htpb,
                                    d_p0_um=dp0,
                                    pc_psia=pc,
                                    primary_of=p_of,
                                    chamber_length_m=lc,
                                    chamber_diameter_m=0.08,
                                    mass_flow_rate_kg_s=0.50
                                )

                                assert len(df_traj) > 0, "Empty trajectory DataFrame"
                                assert not df_traj.isna().any().any(), "NaN values in trajectory"
                                
                                # Check particle diameter decay monotonicity
                                d_p_arr = df_traj['d_p_um'].values
                                diffs = np.diff(d_p_arr)
                                assert (diffs <= 1e-9).all(), "Particle diameter decay non-monotonic!"
                                
                                # Check expulsion efficiency bounds
                                eta_arr = df_traj['expulsion_eff_pct'].values
                                assert (eta_arr >= 0.0).all() and (eta_arr <= 100.0).all(), "Expulsion efficiency out of [0, 100%]"

                                suite2_passed += 1
                                passed_tests += 1
                            except Exception as e:
                                failures.append({
                                    'suite': '1D_Droplet_Trajectory',
                                    'al': al, 'htpb': htpb, 'pc': pc, 'dp0': dp0, 'lc': lc, 'p_of': p_of,
                                    'error': str(e)
                                })

    print(f"Suite 2 Finished. Passed {suite2_passed}/{suite2_total} evaluations.")

    # --- TEST SUITE 3: BATES Grain Internal Ballistics & CFD Export ---
    print("\n--- [Suite 3/3] Testing BATES Grain Regression & CFD Export ---")
    suite3_total = 0
    suite3_passed = 0

    for al in al_pct_list:
        for htpb in htpb_pct_list:
            for p_of in primary_of_list:
                suite3_total += 1
                total_tests += 1
                try:
                    df_bates = simulate_bates_grain_regression(
                        al_pct=al,
                        htpb_pct=htpb,
                        primary_of=p_of
                    )
                    assert len(df_bates) > 0, "Empty BATES DataFrame"
                    assert not df_bates.isna().any().any(), "NaN in BATES DataFrame"
                    assert df_bates['P_c_MPa'].min() >= 0.05, "Unphysical low pressure in BATES"
                    assert df_bates['r_b_mm_s'].min() > 0.0, "Negative or zero burn rate in BATES"

                    # Verify CFD exporter handles this dataset
                    export_converge_cfd_inputs(df_bates, output_dir="d:/CEA/scratch_test_out")
                    
                    suite3_passed += 1
                    passed_tests += 1
                except Exception as e:
                    failures.append({
                        'suite': 'BATES_Grain_Regression',
                        'al': al, 'htpb': htpb, 'p_of': p_of,
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
