# Null Calabi–Yau → Parameter-Free Prediction of α<sub>s</sub>(M<sub>Z</sub>)

This repository reproduces a **parameter-free** prediction of the strong coupling at the Z mass,
derived from a fixed Null Calabi–Yau geometry.  
There are **no free knobs, no tuning, and no α<sub>s</sub> input**: the result follows uniquely from the geometry.

---

## Background

- The **Null Calabi–Yau (NCY)** is a special 6-dimensional geometry consistent with string theory.
- From this geometry we derive the **GUT threshold corrections** and **SUSY spectrum thresholds**.
- The only experimental inputs are two electroweak observables at the Z pole:
  - Inverse electromagnetic coupling: 1/α<sub>EM</sub>(M<sub>Z</sub>)
  - Weak mixing angle: sin²θ<sub>W</sub>(M<sub>Z</sub>)
- Everything else is determined by:
  1. Solving μ₀ from the physical condition α₁(μ₀) = α₂(μ₀),
  2. Deterministically deriving heavy and SUSY spectra from the fixed geometry,
  3. Computing thresholds and evolving couplings with the renormalization group.

The pipeline then predicts α<sub>s</sub>(M<sub>Z</sub>) with no freedom to retune.

---

## Requirements

- Python ≥ 3.6  
- Install dependencies:
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

---

## Run (single command)

```bash
python3 predict_alpha_s_from_nullCY.py
```

This will:
1. Solve μ₀ from α₁ = α₂ (piecewise SM→MSSM),
2. Derive heavy & SUSY spectra from `inputs/geometry.json`,
3. Check provenance with `tools/guard.py`,
4. Compute thresholds and run RGEs,
5. Print α<sub>s</sub>(M<sub>Z</sub>) and write results to `results/alpha_s_from_nullCY.json`.

---

## Inputs

- `inputs/geometry.json` – fixed Null Calabi–Yau data  
- `inputs/ew_targets.json` – electroweak inputs at M<sub>Z</sub> (1/α<sub>EM</sub>, sin²θ<sub>W</sub>)

---

## Generated artifacts

These are written automatically and should not be version-controlled:

- `inputs/heavy_spectrum.json`
- `inputs/susy_spectrum.json`

They are deterministically derived from the geometry and μ₀, and their provenance is enforced by a guard.

---

## Verification checklist

1. **Predict α<sub>s</sub>:**
```bash
python3 predict_alpha_s_from_nullCY.py
```
You should see JSON output containing roughly:

```json
{
  "mu0_GeV": 1.44e16,
  "alphaU": 0.0386,
  "Delta_GUT_inv_alpha": { "d1": -0.05, "d2": -0.05, "d3": -1.016 },
  "prediction": { "alpha_s_MZ": 0.1197 }
}
```

2. **Guard provenance:**
```bash
python3 tools/guard.py <mu0_GeV_from_output>
```
Expected:  
```json
{"ok": true, "mu0_checked_GeV": ...}
```

If the artifacts were hand-edited or not derived from the current geometry+code+μ₀, the guard will fail.

---

## Design principles (why this is knob-free)

- **No α<sub>s</sub> input** ever.
- **μ₀ fixed** by the physical condition α₁ = α₂.
- **Heavy thresholds** from a deterministic basis match in SU(5) b-space (XY, SU(2) triplet, SU(3) octet).
- **SUSY thresholds** from deterministic geometry rules.
- **Provenance guard** bans hand-edits or retuning of artifacts.

---

## Result

This pipeline produces a **zero-knob, reproducible prediction** of α<sub>s</sub>(M<sub>Z</sub>) from the Null Calabi–Yau theory alone.  
The central value is within experimental range, with quantified electroweak and spectrum uncertainties.

---

## License

Choose a license (MIT / BSD-3-Clause / Apache-2.0) and add a `LICENSE` file here.

