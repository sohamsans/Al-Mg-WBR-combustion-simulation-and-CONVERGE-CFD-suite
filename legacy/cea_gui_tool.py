"""
cea_gui_tool.py
===============
Al-Mg Composite Propellant CEA Analysis Tool
- Tkinter GUI with table view and controls
- NASA CEA sweep via rocketcea
- Safe species parsing (percent vs fraction detection)
- Slag fraction computation
- Derived thermodynamic properties
- Plotting (MR, Isp, Tc, MW, gamma, slag)
- Export: result.csv + Propellant_Data.dat (CFD-ready)
"""

# ─────────────────────────────────────────────────────────────
# 1.  IMPORTS
# ─────────────────────────────────────────────────────────────
import os
import sys
import csv
import math
import traceback
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from rocketcea.cea_obj_w_units import CEA_Obj
    HAS_CEA = True
except ImportError:
    try:
        from rocketcea.cea_obj import CEA_Obj
        HAS_CEA = True
    except ImportError:
        HAS_CEA = False

# ─────────────────────────────────────────────────────────────
# 2.  CONSTANTS & DEFAULT PARAMETERS
# ─────────────────────────────────────────────────────────────
DEFAULT_PARAMS = {
    "Al_pct":      65.0,   # wt% Aluminum in fuel
    "Mg_pct":      25.0,   # wt% Magnesium in fuel
    "binder_pct":  10.0,   # wt% Binder (HTPB)
    "Pc_bar":      30.0,   # Chamber pressure [bar]
    "MR_min":       3.0,   # O/F mass ratio sweep start (Al/Mg/Water stoich ~6-8)
    "MR_max":      12.0,   # O/F mass ratio sweep end
    "MR_step":      0.5,   # O/F step size
    "eps":          8.0,   # Nozzle area ratio
}

R_UNIVERSAL = 8314.46   # J/(kmol·K)
SLAG_SPECIES = {"AL2O3", "MGO", "AL2O3(L)", "MGO(L)",
                "AL2O3(S)", "MGO(S)", "ALCL3", "MGCL2"}

# ─────────────────────────────────────────────────────────────
# 3.  CEA HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────

def build_cea_obj(al_pct: float, mg_pct: float, binder_pct: float):
    """
    Build a rocketcea CEA_Obj for Al-Mg/Water propellant.

    IMPORTANT: Always uses plain CEA_Obj (cea_obj, NOT cea_obj_w_units).
    Plain CEA_Obj works in psia/imperial. Unit conversion is handled
    in run_cea_sweep() after the fact.

    Card format follows the confirmed working example from:
      https://rocketcea.readthedocs.io/en/latest/new_propellants.html
    Rules:
      - wt%=  (with percent sign)
      - h,cal= and t(k)= on the line BELOW the species line
      - rho= (no .g/cc suffix for Al line in the docs example)
    """
    if not HAS_CEA:
        return None

    total = al_pct + mg_pct + binder_pct
    if abs(total - 100.0) > 0.5:
        raise ValueError(f"Fuel percentages must sum to 100 (got {total:.1f})")

    from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

    # Normalise to exactly 100 wt%
    scale = 100.0 / total
    al   = al_pct   * scale
    mg   = mg_pct   * scale
    bind = binder_pct * scale

    # ── Water oxidizer ────────────────────────────────────────
    # CRITICAL: first token after 'oxid' must EXACTLY match ox_name
    # Hf(H2O liq, 298K) = -68317.4 cal/mol
    ox_name = "MyWater"
    ox_card = (
        f"oxid {ox_name}  H 2.0 O 1.0   wt%=100.00\n"
        f"h,cal=-68317.4     t(k)=298.15   rho.g/cc=1.0\n"
    )
    add_new_oxidizer(ox_name, ox_card)

    # ── Al/Mg/HTPB fuel ───────────────────────────────────────
    # CRITICAL: first token after each 'fuel' does NOT need to match
    # fuel_name — it is a species label. fuel_name is just the
    # handle passed to CEA_Obj(fuelName=...).
    # Al(cr) Hf=0, Mg(cr) Hf=0 (reference-state elements)
    # HTPB (C7.337 H10.982 O0.103): Hf=-1255 cal/mol (JANNAF)
    fuel_name = f"AlMgFuel_{int(al):02d}_{int(mg):02d}"
    fuel_card = (
        f"fuel AL    AL 1.0                     wt%={al:.4f}\n"
        f"h,cal=0.0      t(k)=298.15   rho.g/cc=2.7\n"
        f"fuel MG    Mg 1.0                     wt%={mg:.4f}\n"
        f"h,cal=0.0      t(k)=298.15   rho.g/cc=1.74\n"
        f"fuel HTPB  C 7.337 H 10.982 O 0.103  wt%={bind:.4f}\n"
        f"h,cal=-1255.0  t(k)=298.15   rho.g/cc=0.93\n"
    )
    add_new_fuel(fuel_name, fuel_card)

    # Always use plain CEA_Obj (psia) — cea_obj_w_units silently returns
    # zeros for custom-registered propellants in many rocketcea builds.
    cea = CEA_Obj(oxName=ox_name, fuelName=fuel_name)
    return cea, ox_card, fuel_card   # return cards for diagnostic logging


def safe_extract_species_dict(cea_obj, Pc, MR, eps):
    """
    Return a {species: massFraction} dict, normalised to sum=1.
    Handles rocketcea versions that return list, dict, or tuple.
    """
    result = {}
    try:
        raw = cea_obj.get_SpeciesMassFractions(Pc=Pc, MR=MR, eps=eps)
        if isinstance(raw, dict):
            result = {k: float(v) for k, v in raw.items()}
        elif isinstance(raw, (list, tuple)):
            # Alternating [name, value, ...] or [(name, val), ...]
            flat = []
            for item in raw:
                if isinstance(item, (list, tuple)):
                    flat.extend(item)
                else:
                    flat.append(item)
            it = iter(flat)
            for name in it:
                try:
                    val = next(it)
                    result[str(name)] = float(val)
                except StopIteration:
                    break
    except Exception:
        pass
    return result


def detect_and_normalize_massfractions(species_dict: dict) -> dict:
    """
    Detect if values are in percent (0-100) or fraction (0-1).
    Returns dict with values guaranteed in 0..1 range.
    """
    if not species_dict:
        return {}
    total = sum(species_dict.values())
    if total > 1.5:           # likely percent
        return {k: v / 100.0 for k, v in species_dict.items()}
    return dict(species_dict)


def compute_slag_fraction(species_dict: dict) -> float:
    """Sum mass fractions of condensed oxide species → slag fraction."""
    normalised = detect_and_normalize_massfractions(species_dict)
    return sum(v for k, v in normalised.items()
               if k.upper().split("(")[0].strip() in SLAG_SPECIES)

# ─────────────────────────────────────────────────────────────
# 4.  CEA SWEEP ENGINE
# ─────────────────────────────────────────────────────────────

def run_cea_sweep(params: dict, log_fn=print) -> list[dict]:
    """
    Sweep O/F ratio and extract performance + species data.
    Returns list of row dicts suitable for the table and CSV.

    Notes on rocketcea API:
      - ALL calls use POSITIONAL arguments (not keyword).
      - cea_obj_w_units accepts Pc in the unit declared at construction (Bar).
      - Plain cea_obj always expects Pc in psia — convert Bar→psia if needed.
      - estimate_Ambient_Isp does not exist on all builds; wrapped safely.
    """
    if not HAS_CEA:
        log_fn("[ERROR] rocketcea not installed. Cannot run CEA sweep.")
        return []

    al   = params["Al_pct"]
    mg   = params["Mg_pct"]
    bind = params["binder_pct"]
    Pc_bar = params["Pc_bar"]
    eps    = params["eps"]

    # Plain CEA_Obj always uses psia. Convert once here.
    BAR_TO_PSIA  = 14.5038
    RANKINE_TO_K = 0.5556      # Tc from CEA is in Rankine → Kelvin
    FTPS_TO_MPS  = 0.3048      # Cstar from CEA is ft/s  → m/s
    GCCC_TO_KGCM = 1000.0      # g/cc → kg/m³
    Pc = Pc_bar * BAR_TO_PSIA  # psia

    mr_vals = []
    mr = params["MR_min"]
    while mr <= params["MR_max"] + 1e-9:
        mr_vals.append(round(mr, 6))
        mr += params["MR_step"]

    log_fn(f"Building CEA object (Al={al}%, Mg={mg}%, Binder={bind}%)…")
    try:
        cea, ox_card, fuel_card = build_cea_obj(al, mg, bind)
    except Exception as e:
        log_fn(f"[ERROR] CEA object build failed: {e}")
        log_fn(traceback.format_exc())
        return []

    log_fn(f"CEA obj: plain (psia)  Pc = {Pc:.2f} psia = {Pc_bar:.2f} bar")
    log_fn(f"--- OXIDIZER CARD ---\n{ox_card}")
    log_fn(f"--- FUEL CARD ---\n{fuel_card}")

    # Also dump raw CEA output at probe MR for full visibility
    try:
        raw = cea.get_full_cea_output(Pc=Pc, MR=mr_vals[0], eps=eps, short_output=1)
        log_fn(f"--- RAW CEA OUTPUT (MR={mr_vals[0]}) ---\n{raw}\n---")
    except Exception as e:
        log_fn(f"[WARN] get_full_cea_output failed: {e}")

    # ── List all available get_* methods for this CEA_Obj version ──
    available = [m for m in dir(cea) if m.startswith("get_")]
    log_fn(f"Available get_* methods: {', '.join(available)}")

    # ── Diagnostic probe at first MR ───────────────────────────
    probe_mr = mr_vals[0]
    try:
        _r    = cea.get_IvacCstrTc(Pc, probe_mr, eps)
        _isp  = float(_r[0])
        _cstr = float(_r[1]) * FTPS_TO_MPS
        _tc   = float(_r[2]) * RANKINE_TO_K
        log_fn(f"Probe MR={probe_mr}: Isp={_isp:.2f}s  "
               f"Cstar={_cstr:.1f}m/s  Tc={_tc:.0f}K")
        if _isp == 0.0 and _tc == 0.0:
            log_fn("[ERROR] Probe still zero — check RAW CEA OUTPUT above.")
            return []
    except Exception as e:
        log_fn(f"[ERROR] Probe failed at MR={probe_mr}: {e}")
        log_fn(traceback.format_exc())
        return []

    rows = []
    for MR in mr_vals:
        row  = {"MR": MR}
        errs = []
        try:
            # ── Primary call: get_IvacCstrTc ──────────────────
            # Most robust getter for difficult propellants — returns
            # (Isp_vac[s], Cstar[ft/s], Tc[R]) in a single CEA call.
            IvacCstrTc = cea.get_IvacCstrTc(Pc, MR, eps)
            row["Isp_vac_s"] = float(IvacCstrTc[0])
            row["cstar_ms"]  = float(IvacCstrTc[1]) * FTPS_TO_MPS
            row["Tc_K"]      = float(IvacCstrTc[2]) * RANKINE_TO_K

            # Cf = Isp_vac * g0 / Cstar  (derived — no get_Cf in this build)
            G0 = 9.80665  # m/s²
            cstar = row["cstar_ms"]
            row["Cf"] = (row["Isp_vac_s"] * G0 / cstar) if cstar > 0 else float("nan")

            # ── Delivered Isp (ambient = 1 atm = 14.696 psia) ─
            try:
                amb = cea.estimate_Ambient_Isp(Pc, MR, eps, 14.696)
                row["Isp_del_s"] = float(amb[0]) if isinstance(amb, (list, tuple)) else float(amb)
            except Exception:
                row["Isp_del_s"] = row["Isp_vac_s"]

            # ── Chamber MW and gamma ───────────────────────────
            try:
                mw_gam           = cea.get_Chamber_MolWt_gamma(Pc, MR, eps)
                row["MW_kgkmol"] = float(mw_gam[0])
                row["gamma"]     = float(mw_gam[1])
            except Exception:
                row["MW_kgkmol"] = float("nan")
                row["gamma"]     = float("nan")

            # ── Exit temperature (Rankine → K) ─────────────────
            try:
                temps       = cea.get_Temperatures(Pc, MR, eps)
                row["Te_K"] = float(temps[2]) * RANKINE_TO_K
            except Exception:
                row["Te_K"] = float("nan")

            # ── Chamber density (g/cc → kg/m³) ────────────────
            try:
                row["rho_kgm3"] = float(cea.get_Densities(Pc, MR, eps)[0]) * GCCC_TO_KGCM
            except Exception:
                row["rho_kgm3"] = float("nan")

            # ── Derived thermodynamic properties ──────────────
            MW  = row["MW_kgkmol"]
            gam = row["gamma"]
            # MW from plain CEA_Obj comes in g/mol; w_units gives kg/kmol (same numeric)
            R_s = R_UNIVERSAL / MW if MW > 0 else 0.0      # J/(kg·K)
            cp  = (gam / (gam - 1)) * R_s / 1000.0 if gam > 1 else 0.0   # kJ/(kg·K)
            cv  = cp / gam if gam > 0 else 0.0
            row["R_specific"] = R_s
            row["cp_kJkgK"]   = cp
            row["cv_kJkgK"]   = cv

            # ── Species & slag ─────────────────────────────────
            sp = safe_extract_species_dict(cea, Pc, MR, eps)
            row["SlagFraction"] = compute_slag_fraction(sp)

        except Exception as e:
            errs.append(traceback.format_exc())
            for key in ["Isp_vac_s","Isp_del_s","Tc_K","Te_K","MW_kgkmol",
                        "gamma","cstar_ms","Cf","rho_kgm3",
                        "R_specific","cp_kJkgK","cv_kJkgK","SlagFraction"]:
                row.setdefault(key, float("nan"))

        rows.append(row)
        isp  = row.get("Isp_vac_s", float("nan"))
        tc   = row.get("Tc_K",      float("nan"))
        slag = row.get("SlagFraction", float("nan"))
        if errs:
            log_fn(f"  [WARN] MR={MR:.2f} failed:\n{errs[0].strip()}")
        else:
            log_fn(f"  MR={MR:.2f}  Isp={isp:.1f}s  Tc={tc:.0f}K  Slag={slag:.4f}")

    return rows


# ─────────────────────────────────────────────────────────────
# 5.  CSV  I/O  &  CFD EXPORT
# ─────────────────────────────────────────────────────────────

CSV_COLUMNS = [
    "MR", "Isp_vac_s", "Isp_del_s", "Tc_K", "Te_K",
    "MW_kgkmol", "gamma", "cstar_ms", "Cf", "rho_kgm3",
    "R_specific", "cp_kJkgK", "cv_kJkgK", "SlagFraction",
]

CSV_UNITS = {
    "MR": "–", "Isp_vac_s": "s", "Isp_del_s": "s",
    "Tc_K": "K", "Te_K": "K", "MW_kgkmol": "kg/kmol",
    "gamma": "–", "cstar_ms": "m/s", "Cf": "–",
    "rho_kgm3": "kg/m³", "R_specific": "J/(kg·K)",
    "cp_kJkgK": "kJ/(kg·K)", "cv_kJkgK": "kJ/(kg·K)",
    "SlagFraction": "–",
}


def save_csv(rows: list[dict], filepath: str, log_fn=print):
    """Write rows to a labelled CSV with units header row."""
    if not rows:
        log_fn("[WARN] No data to save.")
        return
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        # Units comment row
        writer.writerow({c: f"[{CSV_UNITS.get(c,'?')}]" for c in CSV_COLUMNS})
        for row in rows:
            clean = {}
            for c in CSV_COLUMNS:
                val = row.get(c, "")
                if isinstance(val, float) and math.isnan(val):
                    clean[c] = "NaN"
                else:
                    clean[c] = f"{val:.6g}" if isinstance(val, float) else val
            writer.writerow(clean)
    log_fn(f"CSV saved → {filepath}")


def load_csv(filepath: str, log_fn=print) -> list[dict]:
    """Load a previously saved result CSV. Skips the units row."""
    rows = []
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i == 0 and "[" in str(list(row.values())):
                continue  # skip units header row
            parsed = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = v
            rows.append(parsed)
    log_fn(f"Loaded {len(rows)} rows from {filepath}")
    return rows


def save_propellant_dat(rows: list[dict], filepath: str,
                        params: dict, log_fn=print):
    """
    Write a CFD-ready Propellant_Data.dat file with:
      - File header with formulation metadata
      - Tabulated MR, Tc, gamma, R_specific, cp, rho, SlagFraction
    Format is compatible with SU2 custom inlet/combustion patches.
    """
    if not rows:
        log_fn("[WARN] No data — .dat not saved.")
        return
    with open(filepath, "w") as f:
        f.write("# Propellant_Data.dat — Al-Mg Water-Ramjet CEA Output\n")
        f.write(f"# Al={params.get('Al_pct','?')}%  "
                f"Mg={params.get('Mg_pct','?')}%  "
                f"Binder={params.get('binder_pct','?')}%\n")
        f.write(f"# Pc={params.get('Pc_bar','?')} bar  "
                f"eps={params.get('eps','?')}\n")
        f.write("# Columns: MR  Tc[K]  gamma  R[J/kgK]  "
                "cp[kJ/kgK]  rho[kg/m3]  SlagFrac\n")
        f.write(f"{'MR':>8s} {'Tc_K':>10s} {'gamma':>8s} "
                f"{'R_sp':>10s} {'cp':>10s} {'rho':>10s} {'Slag':>10s}\n")
        for row in rows:
            def g(k): return row.get(k, float("nan"))
            f.write(
                f"{g('MR'):>8.4f} {g('Tc_K'):>10.2f} {g('gamma'):>8.5f} "
                f"{g('R_specific'):>10.3f} {g('cp_kJkgK'):>10.5f} "
                f"{g('rho_kgm3'):>10.4f} {g('SlagFraction'):>10.6f}\n"
            )
    log_fn(f"Propellant_Data.dat saved → {filepath}")


# ─────────────────────────────────────────────────────────────
# 6.  PLOTTING
# ─────────────────────────────────────────────────────────────

PLOT_CONFIGS = [
    ("MR", "Isp_vac_s",    "O/F Ratio",  "Isp Vacuum [s]",          "Isp vs MR"),
    ("MR", "Isp_del_s",    "O/F Ratio",  "Isp Delivered [s]",       "Delivered Isp vs MR"),
    ("MR", "Tc_K",         "O/F Ratio",  "Chamber Temp [K]",        "Tc vs MR"),
    ("MR", "MW_kgkmol",    "O/F Ratio",  "Mol. Weight [kg/kmol]",   "MW vs MR"),
    ("MR", "gamma",        "O/F Ratio",  "Gamma [–]",               "γ vs MR"),
    ("MR", "SlagFraction", "O/F Ratio",  "Slag Fraction [–]",       "Slag Fraction vs MR"),
    ("MR", "cstar_ms",     "O/F Ratio",  "C* [m/s]",                "C* vs MR"),
    ("MR", "cp_kJkgK",     "O/F Ratio",  "cp [kJ/kg·K]",            "cp vs MR"),
]


def make_plots_window(rows: list[dict], parent: tk.Tk):
    """Open a Toplevel window with a grid of MR cross-correlation plots."""
    if not HAS_MPL:
        messagebox.showerror("Missing Library",
                             "matplotlib is not installed. Cannot plot.")
        return
    if not rows:
        messagebox.showwarning("No Data", "Run a CEA sweep or load a CSV first.")
        return

    # Extract numeric columns safely
    def col(key):
        vals = []
        for r in rows:
            v = r.get(key)
            try:
                f = float(v)
                vals.append(f if not math.isnan(f) else None)
            except (TypeError, ValueError):
                vals.append(None)
        return vals

    win = tk.Toplevel(parent)
    win.title("CEA Sweep — Cross-Correlation Plots")
    win.geometry("1200x800")

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle("Al–Mg/Water Ramjet — CEA Performance Sweep",
                 fontsize=13, fontweight="bold")
    axes_flat = axes.flatten()
    x_data = col("MR")

    for ax, (xk, yk, xl, yl, title) in zip(axes_flat, PLOT_CONFIGS):
        y_data = col(yk)
        pairs = [(x, y) for x, y in zip(x_data, y_data)
                 if x is not None and y is not None]
        if pairs:
            xs, ys = zip(*pairs)
            ax.plot(xs, ys, "o-", linewidth=1.8, markersize=4, color="#1f77b4")
            ax.fill_between(xs, ys, alpha=0.08, color="#1f77b4")
        ax.set_title(title, fontsize=9, fontweight="bold")
        ax.set_xlabel(xl, fontsize=8)
        ax.set_ylabel(yl, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    canvas = FigureCanvasTkAgg(fig, master=win)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Save-plots button
    def _save_all():
        path = filedialog.asksaveasfilename(
            parent=win, defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf")],
            title="Save plot as…",
        )
        if path:
            fig.savefig(path, dpi=180, bbox_inches="tight")
            messagebox.showinfo("Saved", f"Plot saved:\n{path}", parent=win)

    tk.Button(win, text="💾 Save Plot", command=_save_all,
              bg="#2c7bb6", fg="white", padx=8).pack(pady=4)


# ─────────────────────────────────────────────────────────────
# 7.  TKINTER GUI  –  AppGUI class
# ─────────────────────────────────────────────────────────────

class AppGUI:
    """Main application window."""

    BTN_STYLE = {"padx": 6, "pady": 3, "relief": tk.RAISED, "cursor": "hand2"}

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Al–Mg CEA Propellant Analysis Tool")
        self._rows: list[dict] = []
        self._params = dict(DEFAULT_PARAMS)
        self._sort_col: str | None = None
        self._sort_asc: bool = True
        self._sweep_thread: threading.Thread | None = None
        self._build_ui()
        self._build_menus()

    # ── Layout ────────────────────────────────────────────────
    def _build_ui(self):
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)

        self._build_toolbar(row=0)
        self._build_table(row=1)
        self._build_log(row=2)
        self._build_statusbar(row=3)

    def _build_toolbar(self, row):
        bar = tk.Frame(self.root, bg="#2c3e50", pady=4)
        bar.grid(row=row, column=0, sticky="ew")

        # ── Parameter entry widgets ──
        fields = [
            ("Al %",     "Al_pct",     5),
            ("Mg %",     "Mg_pct",     5),
            ("Binder %", "binder_pct", 5),
            ("Pc [bar]", "Pc_bar",     5),
            ("MR min",   "MR_min",     5),
            ("MR max",   "MR_max",     5),
            ("MR step",  "MR_step",    5),
            ("ε",        "eps",        4),
        ]
        self._vars: dict[str, tk.StringVar] = {}
        for label, key, width in fields:
            tk.Label(bar, text=label, bg="#2c3e50", fg="white",
                     font=("Arial", 8)).pack(side=tk.LEFT, padx=(6, 1))
            var = tk.StringVar(value=str(self._params[key]))
            self._vars[key] = var
            tk.Entry(bar, textvariable=var, width=width,
                     font=("Consolas", 9)).pack(side=tk.LEFT, padx=(0, 4))

        # ── Action buttons ──
        tk.Button(bar, text="▶ Run CEA", command=self._run_sweep,
                  bg="#27ae60", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=4)
        tk.Button(bar, text="📂 Load CSV", command=self._load_csv,
                  bg="#2980b9", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=2)
        tk.Button(bar, text="💾 Save CSV", command=self._save_csv,
                  bg="#8e44ad", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=2)
        tk.Button(bar, text="📄 Save .dat", command=self._save_dat,
                  bg="#d35400", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=2)
        tk.Button(bar, text="📊 Plot", command=self._plot,
                  bg="#16a085", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=2)
        tk.Button(bar, text="🗑 Clear", command=self._clear,
                  bg="#7f8c8d", fg="white", **self.BTN_STYLE).pack(side=tk.LEFT, padx=2)

    def _build_table(self, row):
        frame = tk.Frame(self.root)
        frame.grid(row=row, column=0, sticky="nsew", padx=4, pady=2)
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        cols = CSV_COLUMNS
        self._tree = ttk.Treeview(frame, columns=cols, show="headings",
                                  selectmode="browse")
        col_widths = {
            "MR": 55, "Isp_vac_s": 85, "Isp_del_s": 85, "Tc_K": 75,
            "Te_K": 70, "MW_kgkmol": 85, "gamma": 65, "cstar_ms": 80,
            "Cf": 55, "rho_kgm3": 80, "R_specific": 85,
            "cp_kJkgK": 80, "cv_kJkgK": 80, "SlagFraction": 95,
        }
        for c in cols:
            self._tree.heading(
                c,
                text=f"{c}\n[{CSV_UNITS.get(c,'?')}]",
                anchor=tk.CENTER,
                command=lambda _c=c: self._sort_by(_c),
            )
            self._tree.column(c, width=col_widths.get(c, 80),
                              anchor=tk.E, minwidth=50)

        vsb = ttk.Scrollbar(frame, orient=tk.VERTICAL,
                            command=self._tree.yview)
        hsb = ttk.Scrollbar(frame, orient=tk.HORIZONTAL,
                            command=self._tree.xview)
        self._tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self._tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        # Alternating row colours
        self._tree.tag_configure("odd",  background="#f0f4f8")
        self._tree.tag_configure("even", background="#ffffff")
        self._tree.tag_configure("warn", background="#fff3cd")  # NaN rows

        # Right-click context menu
        self._ctx_menu = tk.Menu(self.root, tearoff=0)
        self._ctx_menu.add_command(label="🔍 View Row Detail",
                                   command=self._show_row_detail)
        self._ctx_menu.add_separator()
        self._ctx_menu.add_command(label="📋 Copy Row as CSV",
                                   command=self._copy_row_csv)
        self._ctx_menu.add_separator()
        self._ctx_menu.add_command(label="📊 Plot from here",
                                   command=self._plot)
        self._tree.bind("<Button-3>", self._show_context_menu)
        self._tree.bind("<Double-1>", lambda e: self._show_row_detail())

    def _build_log(self, row):
        frame = tk.Frame(self.root)
        frame.grid(row=row, column=0, sticky="ew", padx=4)
        frame.columnconfigure(0, weight=1)

        hdr = tk.Frame(frame)
        hdr.grid(row=0, column=0, columnspan=2, sticky="ew")
        tk.Label(hdr, text="Log:", anchor="w",
                 font=("Arial", 8, "bold")).pack(side=tk.LEFT)
        tk.Button(hdr, text="✕ Clear Log", font=("Arial", 7),
                  command=self._clear_log, relief=tk.FLAT,
                  fg="#e74c3c", cursor="hand2").pack(side=tk.RIGHT, padx=4)
        tk.Button(hdr, text="💾 Save Log", font=("Arial", 7),
                  command=self._save_log, relief=tk.FLAT,
                  fg="#2980b9", cursor="hand2").pack(side=tk.RIGHT, padx=2)

        self._log_text = tk.Text(frame, height=5, state=tk.DISABLED,
                                 font=("Consolas", 8), bg="#1e1e2e", fg="#cdd6f4",
                                 insertbackground="white", wrap=tk.WORD)
        self._log_text.grid(row=1, column=0, sticky="ew")
        sb = ttk.Scrollbar(frame, command=self._log_text.yview)
        self._log_text.configure(yscrollcommand=sb.set)
        sb.grid(row=1, column=1, sticky="ns")

    def _build_statusbar(self, row):
        self._status_var = tk.StringVar(value="Ready.")
        bar = tk.Label(self.root, textvariable=self._status_var,
                       relief=tk.SUNKEN, anchor=tk.W,
                       font=("Arial", 8), bg="#ecf0f1")
        bar.grid(row=row, column=0, sticky="ew")

    # ── Helpers ───────────────────────────────────────────────
    def _log(self, msg: str):
        self._log_text.configure(state=tk.NORMAL)
        self._log_text.insert(tk.END, msg + "\n")
        self._log_text.see(tk.END)
        self._log_text.configure(state=tk.DISABLED)
        self._status_var.set(msg[:110])
        self.root.update_idletasks()

    def _read_params(self) -> dict:
        p = {}
        for key, var in self._vars.items():
            try:
                p[key] = float(var.get())
            except ValueError:
                p[key] = DEFAULT_PARAMS[key]
                self._log(f"[WARN] Invalid value for '{key}', using default.")
        return p

    def _populate_table(self, rows: list[dict]):
        self._tree.delete(*self._tree.get_children())
        for i, row in enumerate(rows):
            vals = []
            has_nan = False
            for c in CSV_COLUMNS:
                v = row.get(c, "")
                if isinstance(v, float):
                    if math.isnan(v):
                        has_nan = True
                        vals.append("NaN")
                    else:
                        vals.append(f"{v:.5g}")
                else:
                    vals.append(str(v))
            tag = "warn" if has_nan else ("odd" if i % 2 else "even")
            self._tree.insert("", tk.END, values=vals, tags=(tag,))

    # ── Button Callbacks ──────────────────────────────────────
    def _run_sweep(self):
        if not HAS_CEA:
            messagebox.showerror(
                "rocketcea Missing",
                "rocketcea is not installed.\n\n"
                "Install it with:\n  pip install rocketcea"
            )
            return
        if self._sweep_thread and self._sweep_thread.is_alive():
            messagebox.showwarning("Busy", "A sweep is already running. Please wait.")
            return

        self._log("─" * 60)
        self._log("Starting CEA sweep (background thread)…")
        params = self._read_params()
        self._params = params
        self._status_var.set("⏳ Running CEA sweep…")

        def _worker():
            try:
                rows = run_cea_sweep(params, log_fn=self._log)
                # UI updates must happen on main thread
                self.root.after(0, lambda: self._on_sweep_done(rows))
            except Exception as e:
                msg = traceback.format_exc()
                self.root.after(0, lambda: self._log(f"[ERROR] {e}\n{msg}"))
                self.root.after(0, lambda: self._status_var.set("Sweep failed."))

        self._sweep_thread = threading.Thread(target=_worker, daemon=True)
        self._sweep_thread.start()

    def _on_sweep_done(self, rows: list[dict]):
        if rows:
            self._rows = rows
            self._populate_table(rows)
            self._log(f"✔ Sweep complete — {len(rows)} MR points computed.")
            self._status_var.set(f"Sweep complete: {len(rows)} rows.")
        else:
            self._log("[ERROR] Sweep returned no data.")
            self._status_var.set("Sweep returned no data.")

    def _load_csv(self):
        path = filedialog.askopenfilename(
            title="Load CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            rows = load_csv(path, log_fn=self._log)
            if rows:
                self._rows = rows
                self._populate_table(rows)
                self._log(f"Loaded {len(rows)} rows from {os.path.basename(path)}")
        except Exception as e:
            self._log(f"[ERROR] Load failed: {e}")

    def _save_csv(self):
        if not self._rows:
            messagebox.showwarning("No Data", "Nothing to save — run a sweep first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save CSV",
            defaultextension=".csv",
            initialfile="result.csv",
            filetypes=[("CSV files", "*.csv")],
        )
        if path:
            try:
                save_csv(self._rows, path, log_fn=self._log)
            except Exception as e:
                self._log(f"[ERROR] Save CSV failed: {e}")

    def _save_dat(self):
        if not self._rows:
            messagebox.showwarning("No Data", "Nothing to export — run a sweep first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save Propellant_Data.dat",
            defaultextension=".dat",
            initialfile="Propellant_Data.dat",
            filetypes=[("DAT files", "*.dat"), ("All files", "*.*")],
        )
        if path:
            try:
                save_propellant_dat(self._rows, path,
                                    self._params, log_fn=self._log)
            except Exception as e:
                self._log(f"[ERROR] Save .dat failed: {e}")

    def _plot(self):
        make_plots_window(self._rows, self.root)

    def _clear(self):
        if self._rows and not messagebox.askyesno(
                "Clear", f"Clear all {len(self._rows)} rows from the table?"):
            return
        self._rows = []
        self._tree.delete(*self._tree.get_children())
        self._sort_col = None
        self._log("Table cleared.")
        self._status_var.set("Ready.")

    # ── Column Sort ───────────────────────────────────────────
    def _sort_by(self, col: str):
        """Toggle ascending/descending sort on a column header click."""
        if not self._rows:
            return
        if self._sort_col == col:
            self._sort_asc = not self._sort_asc
        else:
            self._sort_col = col
            self._sort_asc = True

        def _key(r):
            v = r.get(col, 0)
            try:
                return (0, float(v)) if not (isinstance(v, float) and math.isnan(v)) else (1, 0)
            except (TypeError, ValueError):
                return (2, str(v))

        self._rows.sort(key=_key, reverse=not self._sort_asc)
        self._populate_table(self._rows)
        arrow = "▲" if self._sort_asc else "▼"
        self._status_var.set(f"Sorted by {col} {arrow}")

    # ── Context Menu ──────────────────────────────────────────
    def _show_context_menu(self, event):
        """Post right-click menu at cursor position."""
        iid = self._tree.identify_row(event.y)
        if iid:
            self._tree.selection_set(iid)
            self._ctx_menu.post(event.x_root, event.y_root)

    def _get_selected_row(self) -> dict | None:
        """Return the data dict for the currently selected treeview row."""
        sel = self._tree.selection()
        if not sel:
            return None
        vals = self._tree.item(sel[0], "values")
        return dict(zip(CSV_COLUMNS, vals))

    # ── Row Detail Popup ──────────────────────────────────────
    def _show_row_detail(self):
        """Open a Toplevel showing all fields for the selected row."""
        row = self._get_selected_row()
        if not row:
            messagebox.showinfo("No Selection", "Select a row first.")
            return

        win = tk.Toplevel(self.root)
        win.title(f"Row Detail — MR = {row.get('MR', '?')}")
        win.geometry("420x460")
        win.resizable(False, False)

        tk.Label(win, text=f"MR = {row.get('MR','?')}",
                 font=("Arial", 11, "bold"), pady=6).pack()

        frm = tk.Frame(win)
        frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=4)

        for i, col in enumerate(CSV_COLUMNS):
            val  = row.get(col, "—")
            unit = CSV_UNITS.get(col, "")
            bg   = "#f0f4f8" if i % 2 == 0 else "#ffffff"
            r_frm = tk.Frame(frm, bg=bg)
            r_frm.pack(fill=tk.X)
            tk.Label(r_frm, text=col, width=16, anchor="w", bg=bg,
                     font=("Consolas", 9, "bold")).pack(side=tk.LEFT, padx=4)
            tk.Label(r_frm, text=str(val), width=18, anchor="e", bg=bg,
                     font=("Consolas", 9)).pack(side=tk.LEFT)
            tk.Label(r_frm, text=f"  {unit}", anchor="w", bg=bg,
                     fg="#7f8c8d", font=("Arial", 8)).pack(side=tk.LEFT)

        tk.Button(win, text="Close", command=win.destroy,
                  width=12).pack(pady=8)

    # ── Copy Row as CSV ───────────────────────────────────────
    def _copy_row_csv(self):
        """Copy the selected row to clipboard as a comma-separated string."""
        row = self._get_selected_row()
        if not row:
            return
        line = ",".join(str(row.get(c, "")) for c in CSV_COLUMNS)
        self.root.clipboard_clear()
        self.root.clipboard_append(line)
        self._status_var.set("Row copied to clipboard.")

    # ── Log Utilities ─────────────────────────────────────────
    def _clear_log(self):
        self._log_text.configure(state=tk.NORMAL)
        self._log_text.delete("1.0", tk.END)
        self._log_text.configure(state=tk.DISABLED)
        self._status_var.set("Log cleared.")

    def _save_log(self):
        path = filedialog.asksaveasfilename(
            title="Save Log",
            defaultextension=".txt",
            initialfile="cea_sweep_log.txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if path:
            content = self._log_text.get("1.0", tk.END)
            with open(path, "w") as f:
                f.write(content)
            self._status_var.set(f"Log saved → {os.path.basename(path)}")

    # ── Top-level Menu Bar ────────────────────────────────────
    def _build_menus(self):
        menubar = tk.Menu(self.root)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="▶  Run CEA Sweep",
                              accelerator="Ctrl+R", command=self._run_sweep)
        file_menu.add_separator()
        file_menu.add_command(label="📂  Load CSV…",
                              accelerator="Ctrl+O", command=self._load_csv)
        file_menu.add_command(label="💾  Save CSV…",
                              accelerator="Ctrl+S", command=self._save_csv)
        file_menu.add_command(label="📄  Save Propellant_Data.dat…",
                              command=self._save_dat)
        file_menu.add_separator()
        file_menu.add_command(label="🗑  Clear Table",  command=self._clear)
        file_menu.add_separator()
        file_menu.add_command(label="✕  Exit",
                              accelerator="Alt+F4",
                              command=self.root.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="📊  Open Plots",
                              accelerator="Ctrl+P", command=self._plot)
        view_menu.add_command(label="🔍  Row Detail",
                              command=self._show_row_detail)
        menubar.add_cascade(label="View", menu=view_menu)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="ℹ  About", command=self._show_about)
        help_menu.add_command(label="🔧  Dependency Check",
                              command=self._show_deps)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

        # Keyboard shortcuts
        self.root.bind("<Control-r>", lambda e: self._run_sweep())
        self.root.bind("<Control-o>", lambda e: self._load_csv())
        self.root.bind("<Control-s>", lambda e: self._save_csv())
        self.root.bind("<Control-p>", lambda e: self._plot())

    # ── About / Deps Dialogs ──────────────────────────────────
    def _show_about(self):
        messagebox.showinfo(
            "About",
            "Al–Mg CEA Propellant Analysis Tool\n"
            "─────────────────────────────────\n"
            "Project: Optimized Metal-Composite Propellant CFD Analysis\n\n"
            "Performs NASA CEA performance sweeps for Al–Mg/Water\n"
            "hydro-reactive ramjet propellants.\n\n"
            "Outputs: result.csv  |  Propellant_Data.dat (SU2-ready)\n"
            "Requires: rocketcea, matplotlib, tkinter"
        )

    def _show_deps(self):
        rcea = "✔ installed" if HAS_CEA else "✘ MISSING  →  pip install rocketcea"
        mpl  = "✔ installed" if HAS_MPL else "✘ MISSING  →  pip install matplotlib"
        messagebox.showinfo(
            "Dependency Check",
            f"rocketcea  : {rcea}\n"
            f"matplotlib : {mpl}\n"
            f"tkinter    : ✔ (built-in)\n\n"
            f"Python     : {sys.version.split()[0]}"
        )


# ─────────────────────────────────────────────────────────────
# 8.  DEPENDENCY DIAGNOSTICS  (printed on import)
# ─────────────────────────────────────────────────────────────

def _print_diagnostics():
    print("=" * 55)
    print("  Al–Mg CEA GUI Tool — Dependency Check")
    print("=" * 55)
    print(f"  rocketcea  : {'✔ available' if HAS_CEA else '✘ MISSING  →  pip install rocketcea'}")
    print(f"  matplotlib : {'✔ available' if HAS_MPL else '✘ MISSING  →  pip install matplotlib'}")
    print("=" * 55)


# ─────────────────────────────────────────────────────────────
# 9.  ENTRY POINT
# ─────────────────────────────────────────────────────────────

def main():
    _print_diagnostics()
    root = tk.Tk()
    root.geometry("1100x620")
    app = AppGUI(root)
    # Auto-load result.csv if it exists in the current directory
    default_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.csv")
    if os.path.isfile(default_csv):
        app._log(f"Auto-loading {default_csv}…")
        try:
            rows = load_csv(default_csv, log_fn=app._log)
            if rows:
                app._rows = rows
                app._populate_table(rows)
                app._log(f"Auto-load complete: {len(rows)} rows.")
        except Exception as e:
            app._log(f"[WARN] Auto-load failed: {e}")
    root.mainloop()


if __name__ == "__main__":
    main()