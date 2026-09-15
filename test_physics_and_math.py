"""
Automated Physics & Mathematical Edge-Case Test Suite
=====================================================
Tests and validates:
  1. Clamping & boundary limits for alloy properties (0% to 100% Al, out-of-bound inputs).
  2. Thermochemical equilibrium gas state under low/high pressures and extreme O/F ratios.
  3. 1D Droplet dynamics & expulsion efficiency:
     - Monotonic mass & diameter decrease (conservation of mass).
     - Expulsion efficiency strictly bounded in [0%, 100%].
     - Slag accumulation within stoichiometric bounds.
     - Complete absence of NaNs, Infs, and zero-division errors.
  4. BATES grain regression internal ballistics:
     - Realistic chamber pressure and burn rate curves.
     - Web burnout termination without division-by-zero.
     - Invalid geometry error detection (D_o <= d_i0, d_t <= 0).
  5. CONVERGE CFD export file generation & formatting verification.
"""

import sys
import os
import unittest
import numpy as np
import pandas as pd

from primary_combustion_model import (
    get_alloy_properties,
    calculate_primary_gas_state,
    simulate_1d_primary_combustor
)
from bates_converge_cfd_exporter import (
    simulate_bates_grain_regression,
    export_converge_cfd_inputs
)

class TestPhysicsAndMathHardening(unittest.TestCase):

    def test_alloy_properties_bounds(self):
        """Verify alloy density, cp, and ignition temp are bounded and handle invalid inputs."""
        # 1. Pure Al (1.0)
        rho, cp, t_ign = get_alloy_properties(1.0)
        self.assertAlmostEqual(rho, 2700.0, places=1)
        self.assertAlmostEqual(cp, 900.0, places=1)
        self.assertAlmostEqual(t_ign, 2030.0, places=1)

        # 2. Pure Mg (0.0)
        rho, cp, t_ign = get_alloy_properties(0.0)
        self.assertAlmostEqual(rho, 1738.0, places=1)
        self.assertAlmostEqual(cp, 1020.0, places=1)
        self.assertAlmostEqual(t_ign, 1100.0, places=1)

        # 3. 50/50 Alloy
        rho, cp, t_ign = get_alloy_properties(0.5)
        self.assertTrue(1738.0 < rho < 2700.0)
        self.assertTrue(900.0 < cp < 1020.0)
        self.assertTrue(1100.0 < t_ign < 2030.0)

        # 4. Out-of-bounds input clamping (< 0 and > 1)
        rho_neg, cp_neg, t_ign_neg = get_alloy_properties(-0.5)
        self.assertEqual(rho_neg, 1738.0)
        rho_pos, cp_pos, t_ign_pos = get_alloy_properties(1.5)
        self.assertEqual(rho_pos, 2700.0)

    def test_primary_gas_state_robustness(self):
        """Verify gas properties never divide by zero, produce NaN/Inf, or negative values."""
        pressures = [10.0, 50.0, 200.0, 1000.0]
        of_ratios = [0.05, 0.25, 1.0]

        for p in pressures:
            for of in of_ratios:
                gas = calculate_primary_gas_state(al_pct=75.0, htpb_pct=15.0, pc_psia=p, primary_of=of)
                self.assertTrue(np.isfinite(gas['Tc']), f"Tc not finite at P={p}, OF={of}")
                self.assertTrue(gas['Tc'] > 300.0, f"Tc too low at P={p}, OF={of}")
                self.assertTrue(np.isfinite(gas['rho_g']), f"rho_g not finite at P={p}")
                self.assertTrue(gas['rho_g'] > 0.0, f"rho_g non-positive at P={p}")
                self.assertTrue(gas['Gamma'] > 1.0, f"Gamma <= 1.0 at P={p}")
                self.assertTrue(gas['cp_g'] > 0.0, f"cp_g non-positive at P={p}")
                self.assertTrue(np.isfinite(gas['cp_g']), f"cp_g not finite at P={p}")

    def test_1d_droplet_combustion_bounds(self):
        """Verify droplet simulation outputs stay strictly in physical [0, 100%] bounds."""
        cases = [
            {'al_pct': 0.0, 'd_p0_um': 15.0, 'L_c': 0.10},   # Pure Mg, small droplet
            {'al_pct': 75.0, 'd_p0_um': 50.0, 'L_c': 0.20},  # Standard alloy, medium
            {'al_pct': 100.0, 'd_p0_um': 100.0, 'L_c': 0.40},# Pure Al, large droplet
        ]

        for c in cases:
            df, m = simulate_1d_primary_combustor(
                al_pct=c['al_pct'],
                d_p0_um=c['d_p0_um'],
                chamber_length_m=c['L_c']
            )

            self.assertFalse(df.empty, "Trajectory DataFrame should not be empty")
            self.assertFalse(df.isnull().values.any(), f"NaN found in trajectory for {c}")

            # Verify expulsion efficiency bounds: 0 <= eta <= 100%
            self.assertTrue((df['expulsion_eff_pct'] >= 0.0).all(), "Negative expulsion efficiency found")
            self.assertTrue((df['expulsion_eff_pct'] <= 100.0).all(), "Expulsion efficiency exceeded 100%")
            self.assertTrue(0.0 <= m['exit_expulsion_eff_pct'] <= 100.0, "Exit expulsion efficiency out of bounds")

            # Verify mass monotonicity: droplet mass must never increase
            mass_diff = df['m_p_kg'].diff().dropna()
            self.assertTrue((mass_diff <= 1e-25).all(), "Particle mass increased during combustion!")

            # Verify diameter monotonicity
            diam_diff = df['d_p_um'].diff().dropna()
            self.assertTrue((diam_diff <= 1e-12).all(), "Droplet diameter increased during combustion!")

            # Verify slag accumulation bounds
            self.assertTrue(0.0 <= m['total_slag_mass_pct'] <= 189.0, "Slag exceeded stoichiometric limit")

    def test_bates_grain_regression_bounds(self):
        """Verify BATES internal ballistics, web burnout, and geometry error handling."""
        # 1. Standard nominal run
        df = simulate_bates_grain_regression(
            D_o_m=0.10,
            d_i0_m=0.04,
            L_g0_m=0.15,
            N_seg=2,
            d_t_m=0.018
        )

        self.assertFalse(df.empty, "BATES regression DataFrame should not be empty")
        self.assertFalse(df.isnull().values.any(), "NaN found in BATES regression output")

        # Verify positive and realistic physical parameters
        self.assertTrue((df['P_c_MPa'] > 0.05).all(), "Chamber pressure non-positive")
        self.assertTrue((df['P_c_MPa'] < 50.0).all(), "Chamber pressure exceeds 50 MPa physical ceiling")
        self.assertTrue((df['r_b_mm_s'] > 0.1).all(), "Burn rate non-positive")
        self.assertTrue((df['mass_flux_wall_kg_m2s'] > 0.0).all(), "Wall mass flux non-positive")
        self.assertTrue((df['v_inflow_m_s'] > 0.0).all(), "Inflow velocity non-positive")

        # Verify geometric integrity: inner diameter never exceeds outer diameter
        self.assertTrue((df['d_i_m'] <= 0.10).all(), "Inner port diameter exceeded outer grain diameter")

        # 2. Geometry Error Handling: D_o <= d_i0
        with self.assertRaises(ValueError):
            simulate_bates_grain_regression(D_o_m=0.04, d_i0_m=0.04)

        # 3. Geometry Error Handling: zero throat diameter
        with self.assertRaises(ValueError):
            simulate_bates_grain_regression(d_t_m=0.0)

    def test_converge_cfd_export(self):
        """Verify CONVERGE CFD files are created cleanly and contain no NaNs or empty blocks."""
        df = simulate_bates_grain_regression()
        test_dir = "d:/CEA/test_converge_cfd_out"
        export_converge_cfd_inputs(df, output_dir=test_dir)

        expected_files = [
            "bates_internal_ballistics.csv",
            "converge_inflow_mass_flux.dat",
            "converge_bates_boundary.in",
            "converge_thermo.dat"
        ]

        for fname in expected_files:
            fpath = os.path.join(test_dir, fname)
            self.assertTrue(os.path.exists(fpath), f"Missing exported CONVERGE file: {fname}")
            self.assertTrue(os.path.getsize(fpath) > 0, f"Empty exported CONVERGE file: {fname}")

            # Check for NaN and Inf tokens using word boundaries (avoid matching 'inflow')
            import re
            with open(fpath, 'r') as f:
                content = f.read()
                self.assertFalse(re.search(r'\b(nan|inf|-inf)\b', content, re.IGNORECASE), f"Found NaN/Inf token in {fname}")

if __name__ == '__main__':
    unittest.main(verbosity=2)
