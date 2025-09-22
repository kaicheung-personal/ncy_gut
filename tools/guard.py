#!/usr/bin/env python3
"""
tools/guard.py

Verifies that inputs/heavy_spectrum.json and inputs/susy_spectrum.json:
  - exist
  - contain a 'provenance' block
  - provenance.geometry_sha256 matches the *normalized JSON* hash of inputs/geometry.json
  - provenance.code_sha256 matches derive_spectrum_from_nullCY.py (file hash)
  - provenance.mu0_GeV equals the μ0 passed in
  - provenance.sig == sha256(geometry_sha256 || code_sha256 || f"{mu0:.9e}")

Exits nonzero on failure; prints a JSON object with errors.
"""

import json, hashlib, sys
from pathlib import Path

GEOM   = Path("inputs/geometry.json")
HEAVY  = Path("inputs/heavy_spectrum.json")
SUSY   = Path("inputs/susy_spectrum.json")
DERIVER= Path("derive_spectrum_from_nullCY.py")

def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()

def hash_geometry_normalized(p: Path) -> str:
    """Match the deriver: hash of json.dumps(obj, sort_keys=True).encode('utf-8')."""
    obj = json.loads(p.read_text())
    norm = json.dumps(obj, sort_keys=True).encode("utf-8")
    return sha256_bytes(norm)

def hash_file(p: Path) -> str:
    return sha256_bytes(p.read_bytes())

def check(mu0_expected: float):
    if not GEOM.exists(): raise SystemExit("guard: missing inputs/geometry.json")
    if not HEAVY.exists(): raise SystemExit("guard: missing inputs/heavy_spectrum.json")
    if not SUSY.exists():  raise SystemExit("guard: missing inputs/susy_spectrum.json")

    geom_hash_norm = hash_geometry_normalized(GEOM)
    code_hash      = hash_file(DERIVER)

    ok = True; errs = []

    for path in (HEAVY, SUSY):
        try:
            data = json.loads(path.read_text())
        except Exception as e:
            ok=False; errs.append(f"{path} unreadable JSON: {e}")
            continue
        prov = data.get("provenance", {})
        if not prov:
            ok=False; errs.append(f"{path} missing provenance"); continue

        if prov.get("geometry_sha256") != geom_hash_norm:
            ok=False; errs.append(f"{path} geometry_sha256 mismatch")

        if prov.get("code_sha256") != code_hash:
            ok=False; errs.append(f"{path} code_sha256 mismatch")

        try:
            mu0_found = float(prov.get("mu0_GeV"))
        except Exception:
            mu0_found = float("nan")
        if not (abs(mu0_found - mu0_expected) <= 1e-6):
            ok=False; errs.append(f"{path} mu0_GeV mismatch (found {mu0_found}, expected {mu0_expected})")

        expect_sig = sha256_bytes((geom_hash_norm + code_hash + f"{mu0_expected:.9e}").encode("utf-8"))
        if prov.get("sig") != expect_sig:
            ok=False; errs.append(f"{path} signature mismatch")

    if not ok:
        print(json.dumps({"ok": False, "errors": errs}, indent=2))
        sys.exit(2)

    print(json.dumps({"ok": True, "mu0_checked_GeV": mu0_expected}, indent=2))

def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 tools/guard.py <mu0_GeV>")
    check(float(sys.argv[1]))

if __name__ == "__main__":
    main()

