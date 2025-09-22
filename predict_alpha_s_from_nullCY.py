#!/usr/bin/env python3
"""
predict_alpha_s_from_nullCY.py

Workflow:
  1) read EW inputs from inputs/ew_targets.json
  2) PRE-DERIVE SUSY with a provisional μ0_guess to estimate MSUSY
  3) solve α1=α2 piecewise (SM→MSSM) to get μ0, αU
  4) DERIVE both spectra at this μ0 (calls derive_spectrum_from_nullCY.py --mu0 ... )
  5) run guard on μ0
  6) compute thresholds (Δ_GUT from heavy, Δ_MSUSY from SUSY), 2-loop RG with Yukawas
     • NEW: optionally override Δ_GUT via --thresholds-json or env NCY_GUT_THRESHOLDS_JSON
  7) print αs(MZ)

No αs input, no Δ-matching, no mass scaling scripts. Artifacts MUST carry provenance.
"""

import json, math, subprocess, shlex, os, importlib.util, sys
from pathlib import Path
import argparse  # NEW

# --- constants and helpers
pi = math.pi
MZ  = 91.1876
def g_to_alpha(g): return (g*g)/(4.0*pi)
def alpha_to_g(a): return (4.0*pi*a)**0.5
def ew_to_alpha1_alpha2(inv_alpha_EM_MZ: float, sin2thetaW_MZ: float):
    a_em = 1.0 / inv_alpha_EM_MZ; s2 = float(sin2thetaW_MZ); c2 = 1.0 - s2
    return ( (5.0/3.0) * a_em / c2,  a_em / s2 )

# 1-loop betas:
b1_SM, b2_SM, b3_SM = 41.0/10.0, -19.0/6.0, -7.0
b1_MS, b2_MS, b3_MS = 33.0/5.0, 1.0, -3.0
def inv_alpha_segment(inv_alpha_muA, b, muA, muB): return inv_alpha_muA - (b/(2.0*pi))*math.log(muB/muA)

# import betas (as you used before)
def _import():
    p = os.path.abspath("predict_alpha_s_from_geometry.py")
    spec = importlib.util.spec_from_file_location("mod_fp_noast", p)
    mod = importlib.util.module_from_spec(spec); sys.modules["mod_fp_noast"]=mod
    spec.loader.exec_module(mod)  # type: ignore
    return mod
v = _import()
beta_gauge_SM_withY = v.beta_gauge_SM_withY
beta_gauge_MS_withY = v.beta_gauge_MS_withY
beta_yt_SM = v.beta_yt_SM; beta_yb_SM = v.beta_yb_SM; beta_ytau_SM = v.beta_ytau_SM
beta_yt_MS = v.beta_yt_MS; beta_yb_MS = v.beta_yb_MS; beta_ytau_MS = v.beta_ytau_MS
yt_yb_ytau_UV_from_geometry = v.yt_yb_ytau_UV_from_geometry

# thresholds (your patched exact b-vectors)
import thresholds_kl as TH

def inv_alpha_DR_to_MS(inv_a_DR, i):
    if i==0: return inv_a_DR
    if i==1: return inv_a_DR + 2.0/(12.0*pi)  # SU(2)
    return inv_a_DR + 3.0/(12.0*pi)           # SU(3)

def rk4_step_full(g, y, dln, f_g, f_y):
    def F(_g,_y):
        dg = f_g(_g[0],_g[1],_g[2], _y[0],_y[1],_y[2])
        dy = f_y(_g,_y)
        return dg, dy
    k1g,k1y = F(g,y)
    g2  = [g[i] + 0.5*dln*k1g[i] for i in range(3)]
    y2  = [y[i] + 0.5*dln*k1y[i] for i in range(3)]
    k2g,k2y = F(g2,y2)
    g3  = [g[i] + 0.5*dln*k2g[i] for i in range(3)]
    y3  = [y[i] + 0.5*dln*k2y[i] for i in range(3)]
    k3g,k3y = F(g3,y3)
    g4  = [g[i] + dln*k3g[i] for i in range(3)]
    y4  = [y[i] + dln*k3y[i] for i in range(3)]
    k4g,k4y = F(g4,y4)
    g_next = [g[i] + (dln/6.0)*(k1g[i] + 2*k2g[i] + 2*k3g[i] + k4g[i]) for i in range(3)]
    y_next = [y[i] + (dln/6.0)*(k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i]) for i in range(3)]
    return g_next, y_next

def evolve_block(mu_hi, g_hi, y_hi, mu_lo, regime):
    if regime=="MSSM":
        def f_g(g1,g2,g3, yt,yb,ytau): return beta_gauge_MS_withY(g1,g2,g3, yt,yb,ytau)
        def f_y(G,Y): g1,g2,g3=G; yt,yb,ytau=Y; return [
            beta_yt_MS(g1,g2,g3, yt,yb,ytau), beta_yb_MS(g1,g2,g3, yt,yb,ytau), beta_ytau_MS(g1,g2,g3, yt,yb,ytau)]
    else:
        def f_g(g1,g2,g3, yt,yb,ytau): return beta_gauge_SM_withY(g1,g2,g3, yt,yb,ytau)
        def f_y(G,Y): g1,g2,g3=G; yt,yb,ytau=Y; return [
            beta_yt_SM(g1,g2,g3, yt,yb,ytau), beta_yb_SM(g1,g2,g3, yt,yb,ytau), beta_ytau_SM(g1,g2,g3, yt,yb,ytau)]
    t_hi,t_lo = math.log(mu_hi), math.log(mu_lo); dln=(t_lo-t_hi)/4600.0
    g=list(g_hi); y=list(y_hi); mu=mu_hi
    for _ in range(4600):
        g,y = rk4_step_full(g,y,dln,f_g,f_y); mu *= math.exp(dln)
    return g,y

def solve_mu0_and_alphaU_piecewise_1loop(a1_MZ, a2_MZ, MSUSY):
    inv_a1_MZ, inv_a2_MZ = 1.0/a1_MZ, 1.0/a2_MZ
    inv_a1_MS = inv_alpha_segment(inv_a1_MZ, b1_SM, MZ, MSUSY)
    inv_a2_MS = inv_alpha_segment(inv_a2_MZ, b2_SM, MZ, MSUSY)
    t = (inv_a1_MS - inv_a2_MS)/((b1_MS - b2_MS)/(2.0*pi))
    mu0 = MSUSY * math.exp(t)
    inv_a1_mu0 = inv_a1_MS - (b1_MS/(2.0*pi))*t
    return mu0, 1.0/inv_a1_mu0

def _call(cmd: str):
    try:
        return subprocess.check_output(shlex.split(cmd), stderr=subprocess.STDOUT).decode()
    except subprocess.CalledProcessError as e:
        print("\n--- subprocess failed ----------------------------------")
        print("CMD:", cmd)
        print("EXIT:", e.returncode)
        print("STDOUT/ERR:")
        try:
            print(e.output.decode())
        except Exception:
            print(e.output)
        print("--------------------------------------------------------\n")
        raise

def _load_artifact_entries(path: Path):
    data = json.loads(path.read_text())
    return data["entries"], data["provenance"]

# NEW: optional thresholds override
def _maybe_load_thresholds_override(path: Path):
    if path.exists():
        try:
            obj = json.loads(path.read_text())
            D = obj.get("Delta_GUT_inv_alpha") or obj  # allow raw dict or the wrapper
            d1 = float(D["d1"]); d2 = float(D["d2"]); d3 = float(D["d3"])
            return [d1, d2, d3]
        except Exception as e:
            print(f"[warn] Failed to parse thresholds override at {path}: {e}")
    return None

def main():
    # CLI: allow --thresholds-json; also honor env NCY_GUT_THRESHOLDS_JSON
    ap = argparse.ArgumentParser()
    ap.add_argument("--thresholds-json", type=str, default=os.environ.get("NCY_GUT_THRESHOLDS_JSON", ""))
    args = ap.parse_args()
    override_path = Path(args.thresholds_json) if args.thresholds_json else None

    # 1) EW inputs
    ew = json.loads(Path("inputs/ew_targets.json").read_text())
    inv_alpha_EM = float(ew["inv_alpha_EM_MZ"]); sin2w = float(ew["sin2thetaW_MZ"])
    a1_MZ, a2_MZ = ew_to_alpha1_alpha2(inv_alpha_EM, sin2w)

    # --- pre-derive SUSY with a provisional μ0 to estimate MSUSY ---
    mu0_guess = 1.5e16
    _call(f"python3 derive_spectrum_from_nullCY.py --mu0 {mu0_guess}")
    susy_entries, _ = _load_artifact_entries(Path("inputs/susy_spectrum.json"))
    masses_sp = [float(e["mass_GeV"]) for e in susy_entries if float(e["mass_GeV"])>0]
    MSUSY = math.exp(sum(math.log(m) for m in masses_sp)/len(masses_sp)) if masses_sp else 1000.0

    # 2) Solve μ0, αU
    mu0, alphaU = solve_mu0_and_alphaU_piecewise_1loop(a1_MZ, a2_MZ, MSUSY)

    # 3) Derive BOTH spectra at this μ0; then guard
    _call(f"python3 derive_spectrum_from_nullCY.py --mu0 {mu0}")
    _call(f"python3 tools/guard.py {mu0}")

    # 4) Load spectra (entries only)
    heavy_entries, _ = _load_artifact_entries(Path("inputs/heavy_spectrum.json"))
    susy_entries,  _ = _load_artifact_entries(Path("inputs/susy_spectrum.json"))

    # 5) Build Δ thresholds (GUT: from spectrum unless overridden; SUSY: from spectrum)
    Delta_GUT = None
    if override_path:
        Delta_GUT = _maybe_load_thresholds_override(override_path)
        if Delta_GUT:
            print(f"[info] Using Δ_GUT override from {override_path}")
    if Delta_GUT is None:
        # default: compute from heavy spectrum via thresholds_kl
        Delta_GUT = TH.delta_inverse_alpha_from_spectrum(mu0, heavy_entries)

    Delta_MS  = TH.delta_MSUSY_inverse_alpha(MSUSY, susy_entries)

    # 6) UV boundary at μ0 including Δ_GUT (inverse-alpha)
    inv_a_mu0 = 1.0/alphaU
    inv_a_mu0_plus = [inv_a_mu0 + Delta_GUT[i] for i in range(3)]
    g_mu0 = [alpha_to_g(1.0/x) for x in inv_a_mu0_plus]

    # UV Yukawas (same deterministic mapping as before)
    yt0, yb0, ytau0, tb = yt_yb_ytau_UV_from_geometry(json.loads(Path("inputs/geometry.json").read_text()))

    # MSSM down to MSUSY
    g_hi, y_hi = evolve_block(mu0, g_mu0, [yt0,yb0,ytau0], MSUSY, regime="MSSM")

    # DR->MS + SUSY thresholds
    inv_a_MS_plus  = [inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[0]),0),
                      inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[1]),1),
                      inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[2]),2)]
    inv_a_MS_minus = [inv_a_MS_plus[i] + Delta_MS[i] for i in range(3)]
    g_MS_minus = [alpha_to_g(1.0/x) for x in inv_a_MS_minus]

    # SM down to MZ
    g_lo, y_lo = evolve_block(MSUSY, g_MS_minus, y_hi, MZ, regime="SM")

    out = {
      "mode": "first_principles_thresholds (no alpha_s input)",
      "inputs": ew,
      "mu0_GeV": mu0,
      "alphaU": alphaU,
      "MSUSY_GeV": MSUSY,
      "Delta_GUT_inv_alpha": {"d1": Delta_GUT[0], "d2": Delta_GUT[1], "d3": Delta_GUT[2]},
      "Delta_MSUSY_inv_alpha": {"d1": Delta_MS[0], "d2": Delta_MS[1], "d3": Delta_MS[2]},
      "prediction": {"alpha_s_MZ": g_to_alpha(g_lo[2])}
    }
    Path("results").mkdir(exist_ok=True)
    Path("results/alpha_s_from_nullCY.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

if __name__=="__main__":
    main()

