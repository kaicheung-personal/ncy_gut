# --- Group theory helpers (exact presets for the families we use) ---
# All b_i below are the one-loop beta-coefficient contributions per supermultiplet,
# in SU(5) normalization for U(1)_Y. Threshold shifts use c_i = b_i / (2π) * ln(μ/M).

import math
pi = math.pi

def b_vector_for_entry(entry):
    t = entry["type"]
    rep = entry["sm_rep"]

    # 1) Heavy XY gauge bosons (vector multiplet): canonical SU(5) result
    if t == "vector" and rep.get("adjoint_of") == "GUT_XY":
        # contributes negatively to all three:
        return (-5.0/3.0, -3.0, -2.0)

    # 2) Adjoint chiral fragments we use in the fitter:
    #    SU(2) triplet (color singlet), SU(3) octet (weak singlet)
    name = (rep.get("name") or "").lower()
    if t == "chiral":
        # Explicit tag by structure (preferred)
        if (rep.get("dim_SU2") == 3 and rep.get("dim_SU3") == 1 and float(rep.get("Y", 0.0)) == 0.0):
            # SU(2) triplet chiral
            return (0.0, +2.0, 0.0)
        if (rep.get("dim_SU2") == 1 and rep.get("dim_SU3") == 8 and float(rep.get("Y", 0.0)) == 0.0):
            # SU(3) octet chiral
            return (0.0, 0.0, +3.0)

    # 3) Fallback: very generic estimator (use sparingly)
    # This is only for reps not covered above.
    if t == "chiral":
        Y  = float(rep.get("Y", 0.0))
        d2 = int(rep.get("dim_SU2", 1))
        d3 = int(rep.get("dim_SU3", 1))
        # U(1): SU(5) normalization
        b1 = (3.0/5.0) * (Y*Y) * d2 * d3
        # SU(2): T(R) for fundamental doublet is 1/2; for triplet adjoint is 2
        if d2 == 1: T2 = 0.0
        elif d2 == 2: T2 = 0.5 * d3
        elif d2 == 3: T2 = 2.0 * d3
        else: T2 = 0.0
        # SU(3): T(R) for fund is 1/2; for octet adjoint is 3
        if d3 == 1: T3 = 0.0
        elif d3 == 3: T3 = 0.5 * d2
        elif d3 == 8: T3 = 3.0 * d2
        else: T3 = 0.0
        return (b1, T2, T3)

    # 4) Any other vector types default to zero (extend as needed)
    return (0.0, 0.0, 0.0)

def delta_inverse_alpha_from_spectrum(mu_ref, entries):
    d1=d2=d3=0.0
    for e in entries:
        M = float(e["mass_GeV"])
        if M <= 0: continue
        b1,b2,b3 = b_vector_for_entry(e)
        ln = math.log(mu_ref/M)
        d1 += (b1/(2.0*pi))*ln
        d2 += (b2/(2.0*pi))*ln
        d3 += (b3/(2.0*pi))*ln
        if "finite_correction" in e and e["finite_correction"]:
            f1,f2,f3 = e["finite_correction"]
            d1 += f1; d2 += f2; d3 += f3
    return (d1,d2,d3)

def delta_MSUSY_inverse_alpha(MSUSY, sparticle_entries):
    return delta_inverse_alpha_from_spectrum(MSUSY, sparticle_entries)

