#!/usr/bin/env python3
"""
compare_mco.py — Compare two MCML output (.mco) files.

Usage:
    python3 compare_mco.py <cpu_file.mco> <gpu_file.mco>

Reports:
  - RAT values (Rsp, Rd, A, Tt) side by side with relative error
  - Energy conservation check (Rsp + Rd + A + Tt ≈ 1)
  - Pearson correlation of Rd_r and Tt_r radial profiles
  - Pass/Fail verdict
"""

import sys
import re
import numpy as np


def parse_mco(filepath):
    """Parse a minimal set of quantities from a .mco ASCII file."""
    result = {}
    with open(filepath, "r") as fh:
        lines = fh.readlines()

    i = 0
    n = len(lines)

    def next_data_line(start):
        """Skip blank / comment lines and return next non-trivial line index."""
        j = start
        while j < n:
            s = lines[j].strip()
            if s and not s.startswith("#"):
                return j
            j += 1
        return j

    while i < n:
        line = lines[i].strip()

        # --- RAT section ---
        if line.startswith("RAT"):
            i += 1
            vals = []
            while len(vals) < 4 and i < n:
                m = re.match(r"\s*([-+\d.eE]+)", lines[i])
                if m:
                    vals.append(float(m.group(1)))
                i += 1
            if len(vals) == 4:
                result["Rsp"] = vals[0]
                result["Rd"]  = vals[1]
                result["A"]   = vals[2]
                result["Tt"]  = vals[3]
            continue

        # --- Rd_r radial reflectance ---
        if line.startswith("Rd_r"):
            arr = []
            i += 1
            while i < n:
                s = lines[i].strip()
                if not s or (s and not re.match(r"^\s*[-+\d.eE]", s)):
                    break
                m = re.match(r"\s*([-+\d.eE]+)", s)
                if m:
                    arr.append(float(m.group(1)))
                i += 1
            result["Rd_r"] = np.array(arr)
            continue

        # --- Tt_r radial transmittance ---
        if line.startswith("Tt_r"):
            arr = []
            i += 1
            while i < n:
                s = lines[i].strip()
                if not s or (s and not re.match(r"^\s*[-+\d.eE]", s)):
                    break
                m = re.match(r"\s*([-+\d.eE]+)", s)
                if m:
                    arr.append(float(m.group(1)))
                i += 1
            result["Tt_r"] = np.array(arr)
            continue

        i += 1

    return result


def rel_err(a, b):
    """Relative error of b w.r.t. a (in %)."""
    if a == 0:
        return float("inf") if b != 0 else 0.0
    return abs(b - a) / abs(a) * 100.0


def pearson(a, b):
    """Pearson correlation coefficient."""
    if len(a) != len(b) or len(a) == 0:
        return float("nan")
    c = np.corrcoef(a, b)
    return c[0, 1]


def verdict(passed):
    return "✓ PASS" if passed else "✗ FAIL"


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_mco.py <cpu.mco> <gpu.mco>")
        sys.exit(1)

    cpu_path = sys.argv[1]
    gpu_path = sys.argv[2]

    print(f"\nParsing CPU result : {cpu_path}")
    cpu = parse_mco(cpu_path)
    print(f"Parsing GPU result : {gpu_path}")
    gpu = parse_mco(gpu_path)

    print()
    print("=" * 60)
    print("  MCML CPU vs GPU Comparison Report")
    print("=" * 60)

    # ---- RAT scalars ----
    rat_keys = ["Rsp", "Rd", "A", "Tt"]
    print(f"\n{'Quantity':<10} {'CPU':>14} {'GPU':>14} {'Rel.Err(%)':>12}  Status")
    print("-" * 60)

    all_pass = True
    thresholds = {"Rsp": 1e-6, "Rd": 2.0, "A": 2.0, "Tt": 2.0}  # % relative error
    for k in rat_keys:
        cv = cpu.get(k, float("nan"))
        gv = gpu.get(k, float("nan"))
        re_pct = rel_err(cv, gv)
        thr = thresholds.get(k, 2.0)
        ok = re_pct < thr
        if k == "Rsp":
            # Rsp is deterministic — should be identical
            ok = abs(cv - gv) < 1e-9
        all_pass = all_pass and ok
        print(f"  {k:<8} {cv:>14.6g} {gv:>14.6g} {re_pct:>11.4f}%  {verdict(ok)}")

    # ---- Energy conservation ----
    print()
    for label, data in [("CPU", cpu), ("GPU", gpu)]:
        rsp = data.get("Rsp", 0)
        rd  = data.get("Rd",  0)
        a   = data.get("A",   0)
        tt  = data.get("Tt",  0)
        total = rsp + rd + a + tt
        ok = abs(total - 1.0) < 1e-3
        all_pass = all_pass and ok
        print(f"  Energy ({label}): Rsp+Rd+A+Tt = {total:.6f}  {verdict(ok)}")

    # ---- Radial profiles ----
    print()
    for arr_key in ["Rd_r", "Tt_r"]:
        ca = cpu.get(arr_key)
        ga = gpu.get(arr_key)
        if ca is not None and ga is not None and len(ca) > 0 and len(ga) > 0:
            corr = pearson(ca, ga)
            ok = corr > 0.995
            all_pass = all_pass and ok
            print(f"  {arr_key} Pearson r = {corr:.6f}  {verdict(ok)}")
        else:
            print(f"  {arr_key}: not found in one or both files")

    # ---- Overall verdict ----
    print()
    print("=" * 60)
    if all_pass:
        print("  OVERALL: ✓ PASS — GPU results are statistically consistent with CPU")
    else:
        print("  OVERALL: ✗ FAIL — Significant discrepancy detected")
    print("=" * 60)
    print()

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
