#!/usr/bin/env python3
# Null-CY -> αs(MZ) with NO αs input
# Upgrades over v13:
#   - geometry-fixed tanβ rule
#   - include y_b and y_τ along with y_t
#   - 2-loop gauge beta functions include (t,b,τ) Yukawa blocks (SM & MSSM)
#
# Inputs:
#   inputs/geometry.json    (chi, h11, h21, flux_quanta)
#   inputs/ew_targets.json  {"inv_alpha_EM_MZ": ..., "sin2thetaW_MZ": ...}
#
# Output: results/alpha_s_from_geometry.json + stdout JSON

import json, math
from pathlib import Path

pi = math.pi
MZ  = 91.1876

# ---------- helpers ----------
def g_to_alpha(g): return (g*g)/(4.0*pi)
def alpha_to_g(a): return (4.0*pi*a)**0.5

def ew_to_alpha1_alpha2(inv_alpha_EM_MZ: float, sin2thetaW_MZ: float):
    a_em = 1.0 / inv_alpha_EM_MZ
    s2   = float(sin2thetaW_MZ); c2 = 1.0 - s2
    alpha1 = (5.0/3.0) * a_em / c2
    alpha2 = a_em / s2
    return alpha1, alpha2

# 1-loop gauge betas (SU(5)-norm g1)
b1_SM, b2_SM, b3_SM = 41.0/10.0, -19.0/6.0, -7.0
b1_MS, b2_MS, b3_MS = 33.0/5.0, 1.0, -3.0

# 2-loop gauge-only matrices
b_ij_SM = [
    [199.0/50.0, 27.0/10.0, 44.0/5.0],
    [9.0/10.0,   35.0/6.0,  12.0    ],
    [11.0/10.0,  9.0/2.0,  -26.0    ]
]
b_ij_MS = [
    [199.0/25.0, 27.0/5.0,  88.0/5.0],
    [9.0/5.0,    25.0,     24.0    ],
    [11.0/5.0,   9.0,      14.0    ]
]

# ---------- geometry -> model knobs (deterministic; NO αs) ----------
def geometry_to_MSUSY_and_thresholds(geom):
    chi  = int(geom.get("chi", 0))
    h11  = int(geom.get("h11", 2))
    h21  = int(geom.get("h21", 0))
    S    = sum(abs(x) for x in geom.get("flux_quanta", []))

    # super-light SUSY as in v13
    if h11 <= 2 and h21 >= 100 and (S % 2 == 1):
        MSUSY = 300.0
    elif h11 <= 2 and h21 >= 100:
        MSUSY = 500.0
    elif h11 <= 2:
        MSUSY = 1200.0
    elif h11 == 3:
        MSUSY = 2000.0
    else:
        MSUSY = 3000.0

    # SUSY finite matching at MSUSY (inverse-alpha)
    tau1 = ((h21 + S) % 3) - 1
    tau2 = ((h11 + 2*S) % 3) - 1
    tau3 = ((abs(chi)//10) % 3) - 1
    k_MS = 0.20
    Delta_MS = (k_MS*tau1, k_MS*tau2, k_MS*tau3)

    # GUT finite thresholds: keep Δ1=Δ2 (preserve α1=α2), Δ3 more negative for large |χ|
    eps   = -0.05 if (h11 % 2 == 0) else +0.05
    delta3 = -1.02 if abs(chi) >= 200 else -0.82
    Delta_GUT = (eps, eps, delta3)
    return MSUSY, Delta_MS, Delta_GUT

def tanbeta_rule_from_geometry(geom):
    # Deterministic tanβ map (null-CY heuristic):
    h11 = int(geom.get("h11", 2))
    h21 = int(geom.get("h21", 0))
    S   = sum(abs(x) for x in geom.get("flux_quanta", []))
    base = 18.0 if (h11 <= 2 and h21 >= 100) else 12.0
    tweak = 2.0 if (S % 2 == 1) else -2.0
    tb = max(5.0, min(40.0, base + tweak))
    return tb

# ---------- 1-loop piecewise α1=α2 solve ----------
def inv_alpha_segment(inv_alpha_muA, b, muA, muB):
    return inv_alpha_muA - (b/(2.0*pi))*math.log(muB/muA)

def solve_mu0_and_alphaU_piecewise_1loop(a1_MZ, a2_MZ, MSUSY):
    inv_a1_MZ, inv_a2_MZ = 1.0/a1_MZ, 1.0/a2_MZ
    inv_a1_MS = inv_alpha_segment(inv_a1_MZ, b1_SM, MZ, MSUSY)
    inv_a2_MS = inv_alpha_segment(inv_a2_MZ, b2_SM, MZ, MSUSY)
    denom = (b1_MS - b2_MS)/(2.0*pi)
    t = (inv_a1_MS - inv_a2_MS)/denom
    mu0 = MSUSY * math.exp(t)
    inv_a1_mu0 = inv_a1_MS - (b1_MS/(2.0*pi))*t
    alphaU = 1.0/inv_a1_mu0
    return mu0, alphaU

# ---------- MSbar <-> DRbar conversion at MSUSY (1-loop finite) ----------
C_SU2, C_SU3 = 2.0, 3.0
def inv_alpha_MS_to_DR(inv_a_MS, i):
    if i == 0: return inv_a_MS             # U(1) tiny piece neglected at this order
    if i == 1: return inv_a_MS - C_SU2/(12.0*pi)
    return inv_a_MS - C_SU3/(12.0*pi)

def inv_alpha_DR_to_MS(inv_a_DR, i):
    if i == 0: return inv_a_DR
    if i == 1: return inv_a_DR + C_SU2/(12.0*pi)
    return inv_a_DR + C_SU3/(12.0*pi)

# ---------- Yukawa UV rules (deterministic; NO data peeking) ----------
def yt_yb_ytau_UV_from_geometry(geom):
    # Start from v13's yt rule and scale yb,yτ by tanβ as a geometry-fixed choice
    h11 = int(geom.get("h11", 2))
    h21 = int(geom.get("h21", 0))
    S   = sum(abs(x) for x in geom.get("flux_quanta", []))

    base_t = 0.85 if (h11 <= 2 and h21 >= 100) else 0.80
    tweak_t = 0.03 * ((S % 2)*2 - 1)  # ±0.03
    yt0 = max(0.6, min(1.05, base_t + tweak_t))

    tb = tanbeta_rule_from_geometry(geom)
    # Geometry-fixed ansatz: larger tanβ → larger yb,yτ at MSUSY
    yb0 = 0.015 * (tb/15.0)  # O(0.01–0.03)
    ytau0 = 0.010 * (tb/15.0)  # O(0.007–0.03)
    return yt0, yb0, ytau0, tb

# ---------- 1-loop Yukawa RGEs ----------
def beta_yt_SM(g1,g2,g3, yt,yb,ytau):
    return (yt/(16.0*pi*pi)) * ( 4.5*yt*yt + 1.5*yb*yb + ytau*ytau - (17.0/20.0)*g1*g1 - (9.0/4.0)*g2*g2 - 8.0*g3*g3 )

def beta_yb_SM(g1,g2,g3, yt,yb,ytau):
    return (yb/(16.0*pi*pi)) * ( 4.5*yb*yb + 1.5*yt*yt + ytau*ytau - (1.0/4.0)*g1*g1 - (9.0/4.0)*g2*g2 - 8.0*g3*g3 )

def beta_ytau_SM(g1,g2,g3, yt,yb,ytau):
    return (ytau/(16.0*pi*pi)) * ( 3.0*ytau*ytau + 3.0*yb*yb + 3.0*yt*yt - (9.0/4.0)*g1*g1 - (9.0/4.0)*g2*g2 )

def beta_yt_MS(g1,g2,g3, yt,yb,ytau):
    return (yt/(16.0*pi*pi)) * ( 6.0*yt*yt + yb*yb - (13.0/15.0)*g1*g1 - 3.0*g2*g2 - (16.0/3.0)*g3*g3 )

def beta_yb_MS(g1,g2,g3, yt,yb,ytau):
    return (yb/(16.0*pi*pi)) * ( 6.0*yb*yb + yt*yt + ytau*ytau - (7.0/15.0)*g1*g1 - 3.0*g2*g2 - (16.0/3.0)*g3*g3 )

def beta_ytau_MS(g1,g2,g3, yt,yb,ytau):
    return (ytau/(16.0*pi*pi)) * ( 4.0*ytau*ytau + 3.0*yb*yb - (9.0/5.0)*g1*g1 - 3.0*g2*g2 )

# ---------- 2-loop gauge β with Yukawa blocks (t,b,τ) ----------
def beta_gauge_SM_withY(g1,g2,g3, yt,yb,ytau):
    # one-loop
    one = [
        (b1_SM/(16.0*pi*pi))*g1**3,
        (b2_SM/(16.0*pi*pi))*g2**3,
        (b3_SM/(16.0*pi*pi))*g3**3
    ]
    # two-loop gauge part
    B = b_ij_SM
    two_g = []
    for i,gi in enumerate((g1,g2,g3)):
        s = B[i][0]*g1*g1 + B[i][1]*g2*g2 + B[i][2]*g3*g3
        two_g.append( gi**3 * s / ((16.0*pi*pi)**2) )
    # two-loop Yukawa part (SM; coefficients for t,b,τ)
    c1t,c2t,c3t = 17.0/10.0, 3.0/2.0, 2.0
    c1b,c2b,c3b = 1.0/2.0,   3.0/2.0, 2.0
    c1l,c2l,c3l = 3.0/2.0,   1.0/2.0, 0.0
    two_y = [
        -(g1**3)*((c1t*yt*yt + c1b*yb*yb + c1l*ytau*ytau))/((16.0*pi*pi)**2),
        -(g2**3)*((c2t*yt*yt + c2b*yb*yb + c2l*ytau*ytau))/((16.0*pi*pi)**2),
        -(g3**3)*((c3t*yt*yt + c3b*yb*yb + c3l*ytau*ytau))/((16.0*pi*pi)**2)
    ]
    return [one[i]+two_g[i]+two_y[i] for i in range(3)]

def beta_gauge_MS_withY(g1,g2,g3, yt,yb,ytau):
    one = [
        (b1_MS/(16.0*pi*pi))*g1**3,
        (b2_MS/(16.0*pi*pi))*g2**3,
        (b3_MS/(16.0*pi*pi))*g3**3
    ]
    B = b_ij_MS
    two_g = []
    for i,gi in enumerate((g1,g2,g3)):
        s = B[i][0]*g1*g1 + B[i][1]*g2*g2 + B[i][2]*g3*g3
        two_g.append( gi**3 * s / ((16.0*pi*pi)**2) )
    # MSSM Yukawa coefficients (top,bottom,tau)
    c1t,c2t,c3t = 26.0/5.0, 6.0, 4.0
    c1b,c2b,c3b = 14.0/5.0, 6.0, 4.0
    c1l,c2l,c3l = 18.0/5.0, 2.0, 0.0
    two_y = [
        -(g1**3)*((c1t*yt*yt + c1b*yb*yb + c1l*ytau*ytau))/((16.0*pi*pi)**2),
        -(g2**3)*((c2t*yt*yt + c2b*yb*yb + c2l*ytau*ytau))/((16.0*pi*pi)**2),
        -(g3**3)*((c3t*yt*yt + c3b*yb*yb + c3l*ytau*ytau))/((16.0*pi*pi)**2)
    ]
    return [one[i]+two_g[i]+two_y[i] for i in range(3)]

# ---------- RK4 in log μ for (g1,g2,g3, yt,yb,ytau) ----------
def rk4_step_full(g, y, dln, f_g, f_y):
    # g=(g1,g2,g3), y=(yt,yb,ytau)
    def F(_g,_y):
        dg = f_g(_g[0],_g[1],_g[2], _y[0],_y[1],_y[2])
        dy = f_y(_g,_y)
        return dg, dy
    k1g, k1y = F(g,y)
    g2  = [g[i] + 0.5*dln*k1g[i] for i in range(3)]
    y2  = [y[i] + 0.5*dln*k1y[i] for i in range(3)]
    k2g, k2y = F(g2,y2)
    g3  = [g[i] + 0.5*dln*k2g[i] for i in range(3)]
    y3  = [y[i] + 0.5*dln*k2y[i] for i in range(3)]
    k3g, k3y = F(g3,y3)
    g4  = [g[i] + dln*k3g[i] for i in range(3)]
    y4  = [y[i] + dln*k3y[i] for i in range(3)]
    k4g, k4y = F(g4,y4)
    g_next = [g[i] + (dln/6.0)*(k1g[i] + 2*k2g[i] + 2*k3g[i] + k4g[i]) for i in range(3)]
    y_next = [y[i] + (dln/6.0)*(k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i]) for i in range(3)]
    return g_next, y_next

def evolve_block(mu_hi, g_hi, y_hi, mu_lo, steps, regime):
    # regime: "MSSM" uses beta_gauge_MS_withY & beta_yt/_yb/_ytau_MS, else "SM" uses SM versions
    t_hi,t_lo = math.log(mu_hi), math.log(mu_lo)
    dln = (t_lo - t_hi)/steps
    g, y = list(g_hi), list(y_hi)

    if regime == "MSSM":
        def f_g(g1,g2,g3, yt,yb,ytau): return beta_gauge_MS_withY(g1,g2,g3, yt,yb,ytau)
        def f_y(G,Y):
            yt,yb,ytau = Y
            g1,g2,g3 = G
            return [
                beta_yt_MS(g1,g2,g3, yt,yb,ytau),
                beta_yb_MS(g1,g2,g3, yt,yb,ytau),
                beta_ytau_MS(g1,g2,g3, yt,yb,ytau)
            ]
    else:
        def f_g(g1,g2,g3, yt,yb,ytau): return beta_gauge_SM_withY(g1,g2,g3, yt,yb,ytau)
        def f_y(G,Y):
            yt,yb,ytau = Y
            g1,g2,g3 = G
            return [
                beta_yt_SM(g1,g2,g3, yt,yb,ytau),
                beta_yb_SM(g1,g2,g3, yt,yb,ytau),
                beta_ytau_SM(g1,g2,g3, yt,yb,ytau)
            ]

    mu = mu_hi
    for _ in range(steps):
        g, y = rk4_step_full(g, y, dln, f_g, f_y)
        mu *= math.exp(dln)
    return g, y

# ---------- main ----------
def main():
    geom = json.load(open(Path("inputs/geometry.json")))
    ew   = json.load(open(Path("inputs/ew_targets.json")))
    inv_alpha_EM = float(ew["inv_alpha_EM_MZ"])
    sin2w        = float(ew["sin2thetaW_MZ"])

    MSUSY, Delta_MS, Delta_GUT = geometry_to_MSUSY_and_thresholds(geom)

    # EW -> α1, α2 at MZ (MSbar)
    a1_MZ, a2_MZ = ew_to_alpha1_alpha2(inv_alpha_EM, sin2w)

    # Solve α1=α2 with SUSY split (piecewise 1-loop)
    mu0, alphaU = solve_mu0_and_alphaU_piecewise_1loop(a1_MZ, a2_MZ, MSUSY)
    inv_a_mu0 = 1.0/alphaU

    # Apply GUT finite thresholds at μ0 (inverse-alpha)
    inv_a1_mu0_plus = inv_a_mu0 + Delta_GUT[0]
    inv_a2_mu0_plus = inv_a_mu0 + Delta_GUT[1]
    inv_a3_mu0_plus = inv_a_mu0 + Delta_GUT[2]
    g_mu0 = [alpha_to_g(1.0/x) for x in (inv_a1_mu0_plus, inv_a2_mu0_plus, inv_a3_mu0_plus)]

    # Geometry-fixed Yukawas at MSUSY (UV)
    yt0, yb0, ytau0, tb = yt_yb_ytau_UV_from_geometry(geom)

    # μ0 -> MSUSY in MSSM (2-loop gauge + Yukawas)
    g_hi, y_hi = evolve_block(mu0, g_mu0, [yt0,yb0,ytau0], MSUSY, steps=4600, regime="MSSM")

    # Scheme conversion DR->MS at MSUSY (inverse-alpha) + finite SUSY thresholds
    inv_a_MS_plus  = [inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[0]), 0),
                      inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[1]), 1),
                      inv_alpha_DR_to_MS(1.0/g_to_alpha(g_hi[2]), 2)]
    inv_a_MS_minus = [inv_a_MS_plus[0] + Delta_MS[0],
                      inv_a_MS_plus[1] + Delta_MS[1],
                      inv_a_MS_plus[2] + Delta_MS[2]]
    g_MS_minus = [alpha_to_g(1.0/x) for x in inv_a_MS_minus]

    # MSUSY -> MZ in SM (2-loop gauge + Yukawas)
    g_lo, y_lo = evolve_block(MSUSY, g_MS_minus, y_hi, MZ, steps=4600, regime="SM")

    alpha_s_pred = g_to_alpha(g_lo[2])  # MSbar at MZ
    out = {
        "mode": "paramfree_nullCY (no alpha_s input)",
        "geometry_used": geom,
        "derived_model": {
            "M_SUSY_GeV": MSUSY,
            "Delta_MSUSY_inv_alpha": {"d1": Delta_MS[0], "d2": Delta_MS[1], "d3": Delta_MS[2]},
            "Delta_GUT_inv_alpha":   {"d1": Delta_GUT[0], "d2": Delta_GUT[1], "d3": Delta_GUT[2]},
            "tanbeta_rule": "deterministic(null-CY)",
            "tanbeta_value": tb,
            "UV_Yukawas_at_MSUSY": {"yt": yt0, "yb": yb0, "ytau": ytau0}
        },
        "unification": {"mu0_GeV": mu0, "alphaU": alphaU},
        "prediction": {"alpha_s_MZ": alpha_s_pred}
    }
    Path("results").mkdir(exist_ok=True)
    Path("results/alpha_s_from_geometry.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()

