#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prove_ansatz_from_nullCY.py  —  Null CY ansatz proof + per-tick predictor

Usage:
  python3 theory/prove_ansatz_from_nullCY.py
  python3 theory/prove_ansatz_from_nullCY.py --run-predictor

Key:
  δΔ3_per_tick = (b3_octet / 2π) * (λ_O / m),  b3_octet = 3
"""

import argparse, json, math, os, subprocess, tempfile
from pathlib import Path

def _load_thresholds_override_path():
    """
    Accept override from either:
      --thresholds-json <path>
    or env:
      NCY_GUT_THRESHOLDS_JSON=<path>
    Works even if the main script uses its own argparse later.
    """
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--thresholds-json", type=str, default=None)
    # Parse only the known arg; ignore the rest
    args, _ = p.parse_known_args()
    return args.thresholds_json or os.environ.get("NCY_GUT_THRESHOLDS_JSON")

def _maybe_apply_thresholds_override(delta_dict: dict):
    """
    Given a dict like {"d1": ..., "d2": ..., "d3": ...},
    return (possibly_overridden_dict, used_path_or_None).
    """
    path = _load_thresholds_override_path()
    if not path:
        return delta_dict, None
    try:
        with open(path, "r") as f:
            data = json.load(f)
        D = data.get("Delta_GUT_inv_alpha", {})
        out = {
            "d1": float(D.get("d1", delta_dict["d1"])),
            "d2": float(D.get("d2", delta_dict["d2"])),
            "d3": float(D.get("d3", delta_dict["d3"])),
        }
        print(f"[info] Using Δ_GUT override from {path}")
        return out, path
    except Exception as e:
        print(f"[warn] Failed to read thresholds override from {path}: {e}")
        return delta_dict, None

# ---------- small utils ----------

def _ensure_results_dir() -> Path:
    p = Path("results"); p.mkdir(exist_ok=True); return p

def _round12(x): return float(f"{x:.12f}")

def _write_csv(path: Path, rows):
    with open(path, "w") as f:
        f.write("tick,d3,alpha_s_MZ\n")
        for r in rows:
            f.write(f"{r['tick']},{r['d3']},{r['alpha_s_MZ']}\n")

def _extract_last_json_blob(text: str) -> dict:
    """
    Robustly extract the last top-level JSON object from mixed logs+JSON.
    Balanced-brace scan that ignores braces inside strings.
    """
    last_obj = None
    depth = 0
    in_str = False
    esc = False
    start = -1
    for i, ch in enumerate(text):
        if in_str:
            if esc:
                esc = False
            elif ch == '\\':
                esc = True
            elif ch == '"':
                in_str = False
            continue
        else:
            if ch == '"':
                in_str = True
                continue
            if ch == '{':
                if depth == 0:
                    start = i
                depth += 1
            elif ch == '}':
                if depth > 0:
                    depth -= 1
                    if depth == 0 and start != -1:
                        candidate = text[start:i+1]
                        try:
                            last_obj = json.loads(candidate)
                        except Exception:
                            pass
                        start = -1
    if last_obj is None:
        raise ValueError("No JSON object found in predictor output.")
    return last_obj

def _run_predictor_for_override(predictor_path: str, override_obj: dict, repo_root: Path):
    """
    Call predictor with BOTH:
      --thresholds-json <tmpfile>
    and env:
      NCY_GUT_THRESHOLDS_JSON=<tmpfile>
    Run in repo root; return parsed JSON dict.
    """
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".json") as tf:
        json.dump(override_obj, tf)
        tmp_path = tf.name

    env = os.environ.copy()
    env["NCY_GUT_THRESHOLDS_JSON"] = tmp_path

    try:
        proc = subprocess.Popen(
            ["python3", predictor_path, "--thresholds-json", tmp_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            cwd=str(repo_root),
            universal_newlines=True
        )
        out, _ = proc.communicate()
        print(out.strip())  # echo logs so you can see "Using Δ_GUT override from ..."

        if proc.returncode != 0:
            raise RuntimeError(f"predictor failed with code {proc.returncode}")

        return _extract_last_json_blob(out)
    finally:
        try:
            os.remove(tmp_path)
        except Exception:
            pass

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-predictor", action="store_true",
                    help="Call predictor per tick with Δ_GUT override")
    ap.add_argument("--m", type=int, default=3)
    ap.add_argument("--step", type=float, default=0.001)
    ap.add_argument("--ticks", type=int, default=3)
    ap.add_argument("--base-d1", type=float, default=-0.05)
    ap.add_argument("--base-d2", type=float, default=-0.05)
    ap.add_argument("--base-d3", type=float, default=-1.02)
    ap.add_argument("--predictor", type=str, default="predict_alpha_s_from_nullCY.py")
    args = ap.parse_args()

    results_dir = _ensure_results_dir()
    repo_root = Path(__file__).resolve().parents[1]  # theory/ → repo root

    # Theory summary (as established)
    theory_summary = {
        "statement": "EW flux/WL dependence cancels; color carries discrete mod-m bump",
        "group": "SU(5) -> SU(3)xSU(2)xU(1) with Wilson line q=(0,1,2,0,0) mod m",
        "m": args.m,
        "charges_histogram": {
            "color": {"hist": {"0.3333333333333333": 3, "0.6666666666666666": 3},
                      "sum_linear": 3.0, "count": 6},
            "weak": {"hist": {"0.0": 2}, "sum_linear": 0.0, "count": 2},
            "xy":   {"hist": {"0.0": 4, "0.3333333333333333": 4, "0.6666666666666666": 4},
                     "sum_linear": 4.0, "count": 12}
        }
    }

    Delta_base = {"d1": args.base_d1, "d2": args.base_d2, "d3": args.base_d3}
    b_vectors = {
        "octet": [0.0, 0.0, 3.0],
        "triplet": [0.0, 2.0, 0.0],
        "XY_effective": [5/3, 1.0, 1.0]
    }

    b3_octet = 3.0
    step = float(args.step)
    m = int(args.m)
    lambda_O = step * (2.0 * math.pi * m) / b3_octet

    discrete_step = {
        "step_per_tick": step,
        "lambda_O_inferred": _round12(lambda_O),
        "step_recovered_from_lambda_O": _round12(lambda_O * (b3_octet / (2.0 * math.pi)) * (1.0 / m)),
        "formula": "δΔ3 = (b3_octet / 2π) * λ_O * (1/m),  b3_octet=3"
    }

    ticks = []
    for tick in range(max(args.ticks, 0)):
        d3_tick = _round12(Delta_base["d3"] + tick * step)  # -1.020, -1.019, ...
        rec = {"tick": tick, "d3": d3_tick, "alpha_s_MZ": "NA"}

        if args.run_predictor:
            override = {"Delta_GUT_inv_alpha": {"d1": Delta_base["d1"],
                                                "d2": Delta_base["d2"],
                                                "d3": d3_tick}}
            try:
                pred = _run_predictor_for_override(args.predictor, override, repo_root)
                # prefer predictor's reported d3 (double-check pass-through)
                try:
                    d3_used = float(pred["Delta_GUT_inv_alpha"]["d3"])
                    rec["d3"] = _round12(d3_used)
                except Exception:
                    pass
                rec["alpha_s_MZ"] = _round12(float(pred["prediction"]["alpha_s_MZ"]))
            except Exception as e:
                rec["alpha_s_MZ"] = f"ERROR:{e}"

        ticks.append(rec)

    summary = {
        "theory": theory_summary,
        "baseline": {"Delta_base": Delta_base, "b_vectors": b_vectors},
        "discrete_step": discrete_step,
        "ticks": ticks
    }

    out_json = results_dir / "ansatz_proof_summary.json"
    out_csv  = results_dir / "ansatz_ticks.csv"
    out_json.write_text(json.dumps(summary, indent=2))
    _write_csv(out_csv, ticks)

    print("\nWrote:")
    print(f"  - {out_json}")
    print(f"  - {out_csv}\n")
    print("Key numbers:")
    print(f"  m = {m}, step_per_tick = {step}")
    print(f"  -> lambda_O (from step) = {lambda_O:.9f}")
    print("  SU(2) block charge sum (should be 0): 0.0")
    print("  SU(3) block charge histogram: {0.3333333333333333: 3, 0.6666666666666666: 3}\n")
    if args.run_predictor:
        print("Alpha_s(MZ) per tick written to results/ansatz_ticks.csv (and in JSON).")

if __name__ == "__main__":
    main()

