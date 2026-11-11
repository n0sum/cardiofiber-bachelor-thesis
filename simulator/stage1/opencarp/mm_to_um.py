#!/usr/bin/env python3
"""
AugmentA mesh (mm) -> openCARP mesh (um).

Only the coordinates of the .pts are multiplied by 1000; .elem and .lon are
copied unchanged.
"""

import shutil
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------- paths
ROOT    = Path(__file__).resolve().parents[3]                    # repository root
SRC     = ROOT / "data" / "augmenta"  / "LA_endo_with_fiber"     # without extension
DST     = ROOT / "data" / "reference" / "LA_um"
# --------------------------------------------------------------------------


def main():
    pts = np.loadtxt(f"{SRC}.pts", skiprows=1, dtype=np.float64) * 1000.0
    DST.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(f"{DST}.pts", pts, fmt="%.4f", header=str(len(pts)), comments="")
    shutil.copyfile(f"{SRC}.elem", f"{DST}.elem")
    shutil.copyfile(f"{SRC}.lon", f"{DST}.lon")


if __name__ == "__main__":
    main()
