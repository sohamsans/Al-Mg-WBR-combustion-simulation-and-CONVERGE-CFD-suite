"""
CEA Sweep for Al-Mg-HTPB / Water-Ramjet Propellant
====================================================
Uses rocketcea built-in propellants wherever possible:
  - Oxidizer : 'H2O'   (built-in)
  - Fuel     : custom blend of AL (built-in) + MG (custom) + HTPB (built-in)
Only Mg needs add_new_fuel() since it has no native rocketcea entry.

Sweep: MR = 3 → 14, Pc = 435 psia, area ratio = 8.0
"""

import sys
import traceback
import numpy as np
import pandas as pd
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

# ── constants ──────────────────────────────────────────────────────────────
G0_SI       = 9.80665          # m/s²
PSI_TO_PA   = 6894.76          # psi → Pa
PC_PSIA     = 435.11           # chamber pressure [psia]
EPS         = 8.0              # nozzle area ratio
MR_MIN, MR_MAX, MR_STEP = 3.0, 14.0, 0.5

# ── propellant composition (weight fractions of fuel grain) ────────────────
FRAC_AL   = 0.65
FRAC_MG   = 0.25
FRAC_HTPB = 0.10
assert abs(FRAC_AL + FRAC_MG + FRAC_HTPB - 1.0) < 1e-9, "Fuel fractions must sum to 1"

# ── register Mg as a custom fuel (only species not in built-in library) ────
MG_CARD = """\
fuel Magnesium  MG 1.0   wt%=100.00
h,cal=0.0     t(k)=298.15   rho=1.738
"""

def register_Mg():
    try:
        add_new_fuel("Mg_metal", MG_CARD)
        print("[OK]  Mg custom fuel card registered.")
    except Exception as e:
        print(f"[WARN] Mg registration: {e}  (may already be registered)")

# ── build CEA object ────────────────────────────────────────────────────────
def build_cea_obj(al_frac: float, mg_frac: float, htpb_frac: float) -> CEA_Obj:
    """
    Build a CEA_Obj using built-in 'H2O' oxidizer and a blended fuel
    comprising built-in 'AL', custom 'Mg_metal', and built-in 'HTPB'.
    Mixing is done via the fuelName blend notation.
    """
    print(f"\nBuilding CEA object (Al={al_frac*100:.1f}%, "
          f"Mg={mg_frac*100:.1f}%, Binder={htpb_frac*100:.1f}%)…")

    # Blend string for the three-component fuel
    # rocketcea blended fuel: list of (name, wt%) tuples
    from rocketcea.blends import newFuelBlend
    blend_name = newFuelBlend(
        fuelL    = ["AL",       "Mg_metal",   "HTPB"],
        fuelPcentL = [al_frac*100, mg_frac*100, htpb_frac*100]
    )
    print(f"  Fuel blend name : {blend_name}")

    cea = CEA_Obj(
        oxName    = "H2O",
        fuelName  = blend_name,
        isp_units = "sec",
        cstar_units = "m/s",
        pressure_units = "psia",
        temperature_units = "K",
        sonic_velocity_units = "m/s",
        enthalpy_units = "cal/g",
        density_units = "g/cc",
        specific_heat_units = "cal/g-K",
        viscosity_units = "poise",
        thermal_cond_units = "cal/cm-K-s",
    )
    print(f"  CEA obj: {cea}")
    return cea

# ── probe a single MR point ─────────────────────────────────────────────────
def probe_mr(cea: CEA_Obj, mr: float) -> dict | None:
    """Return dict of performance properties at given MR, or None on failure."""
    try:
        isp_vac, cstar, tc = cea.get_IvacCstrTc(
            Pc=PC_PSIA, MR=mr, eps=EPS
        )
    except Exception as e:
        print(f"  [SKIP] MR={mr:.1f}: get_IvacCstrTc failed → {e}")
        return None

    if isp_vac == 0.0 and cstar == 0.0:
        # Dump raw output for diagnosis
        try:
            raw = cea.get_full_cea_output(Pc=PC_PSIA, MR=mr, eps=EPS,
                                          short_output=1, show_transport=0)
            print(f"\n── RAW CEA (MR={mr:.1f}) ──────────────────")
            print(raw[:3000])
            print("────────────────────────────────────────────\n")
        except Exception:
            pass
        print(f"  [ZERO] MR={mr:.1f}: all properties zero — check RAW output above")
        return None

    # Derived quantities
    cf   = (isp_vac * G0_SI) / cstar if cstar > 0 else 0.0

    # Optional extras – wrapped individually
    gamma = None
    try:
        gamma = cea.get_exit_MolWt_gamma(Pc=PC_PSIA, MR=mr, eps=EPS)[1]
    except Exception:
        pass

    mw_exit = None
    try:
        mw_exit = cea.get_exit_MolWt_gamma(Pc=PC_PSIA, MR=mr, eps=EPS)[0]
    except Exception:
        pass

    isp_sl = None
    try:
        isp_sl = cea.estimate_Ambient_Isp(Pc=PC_PSIA, MR=mr, eps=EPS,
                                           Pamb=14.696)[0]
    except Exception:
        pass

    return {
        "MR"       : round(mr, 2),
        "Isp_vac_s": round(isp_vac, 2),
        "Cstar_m/s": round(cstar, 1),
        "Tc_K"     : round(tc, 1),
        "Cf"       : round(cf, 4),
        "gamma"    : round(gamma, 4)   if gamma   else None,
        "MWe_g/mol": round(mw_exit, 2) if mw_exit else None,
        "Isp_SL_s" : round(isp_sl, 2)  if isp_sl  else None,
    }

# ── main sweep ──────────────────────────────────────────────────────────────
def run_sweep():
    register_Mg()
    cea = build_cea_obj(FRAC_AL, FRAC_MG, FRAC_HTPB)

    # Quick probe at MR=6 to verify object before full sweep
    print("\n── Quick probe at MR=6.0 ──")
    qp = probe_mr(cea, 6.0)
    if qp:
        print(f"  Isp={qp['Isp_vac_s']}s  Cstar={qp['Cstar_m/s']}m/s  Tc={qp['Tc_K']}K")
    else:
        print("  [ERROR] Quick probe still zero — dumping raw CEA output…")
        try:
            raw = cea.get_full_cea_output(Pc=PC_PSIA, MR=6.0, eps=EPS,
                                          short_output=1, show_transport=0)
            print(raw[:4000])
        except Exception as e:
            print(f"  get_full_cea_output also failed: {e}")
        print("\n[FATAL] Cannot continue sweep.  See raw output above.")
        return

    print(f"\n── Starting CEA sweep  MR={MR_MIN}→{MR_MAX} step={MR_STEP} ──")
    mr_vals = np.arange(MR_MIN, MR_MAX + MR_STEP/2, MR_STEP)
    rows = []

    for mr in mr_vals:
        result = probe_mr(cea, mr)
        if result:
            rows.append(result)
            print(f"  MR={mr:5.1f}  Isp={result['Isp_vac_s']:7.2f}s  "
                  f"Cstar={result['Cstar_m/s']:7.1f}m/s  Tc={result['Tc_K']:7.1f}K  "
                  f"Cf={result['Cf']:.4f}")

    if not rows:
        print("\n[ERROR] Sweep returned no data.")
        return

    df = pd.DataFrame(rows)
    print(f"\n{'='*65}")
    print(f"Sweep complete — {len(df)} converged points out of {len(mr_vals)}")
    print(f"{'='*65}")
    print(df.to_string(index=False))

    csv_path = "cea_almg_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n[SAVED] Results → {csv_path}")

    # ── Summary statistics ──
    best_isp = df.loc[df["Isp_vac_s"].idxmax()]
    print(f"\n★ Peak Isp_vac = {best_isp['Isp_vac_s']} s  at  MR = {best_isp['MR']}")
    best_tc  = df.loc[df["Tc_K"].idxmax()]
    print(f"★ Peak Tc      = {best_tc['Tc_K']} K     at  MR = {best_tc['MR']}")

if __name__ == "__main__":
    try:
        run_sweep()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
