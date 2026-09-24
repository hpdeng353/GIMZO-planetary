#!/usr/bin/env python3
"""Convert REMIX test initial conditions from SPH-EXA file-init format
(square test, iron/rock KHI) to planetg (GIZMO HDF5) initial conditions.

Conventions handled here:
  - SPH-EXA boxes are centred on the origin (attr "box"); GIZMO expects
    coordinates inside [0, BoxSize). Both REMIX tests span 1 code length in x
    (and y for KHI), so we shift by -box_min per axis and set BoxSize = 1.
    The thin z slab of the KHI keeps its own extent (<< 1); with NOGRAVITY and
    h << gap the z wrap never interacts, so a cubic periodic box is exact.
  - SPH-EXA materialId (62=rock, 63=iron) becomes planetg's imat index into
    EosTableMatIds (sorted unique ids -> 0,1,...); the raw id is preserved in
    EosMaterialId for provenance, mirroring makeplanet_planetg.py.
  - SmoothingLength is informational in planetg (recomputed at startup unless
    READ_HSML is enabled); we carry over the SPH-EXA h as the initial guess.

Usage:
    python remix_ic_to_gizmo.py INPUT.h5 OUTPUT.h5
"""
from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    with h5py.File(args.input, "r") as f:
        g = f["Step#0"]
        x, y, z = (np.asarray(g[k][:], dtype=np.float64) for k in ("x", "y", "z"))
        vx, vy, vz = (np.asarray(g[k][:], dtype=np.float64) for k in ("vx", "vy", "vz"))
        m = np.asarray(g["m"][:], dtype=np.float64)
        u = np.asarray(g["u"][:], dtype=np.float64)
        ids = np.asarray(g["id"][:], dtype=np.uint64)
        h = np.asarray(g["h"][:], dtype=np.float64) if "h" in g else None
        mat = np.asarray(g["materialId"][:], dtype=np.uint32) if "materialId" in g else None
        box = np.asarray(g.attrs["box"], dtype=np.float64) if "box" in g.attrs else None

    n = len(x)
    if box is None:
        box = np.array([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5])
    coords = np.column_stack((x - box[0], y - box[2], z - box[4]))
    lx, ly, lz = box[1] - box[0], box[3] - box[2], box[5] - box[4]
    if abs(lx - ly) > 1e-12 or abs(lx - 1.0) > 1e-9:
        print(f"warning: non-unit/non-square box spans Lx={lx}, Ly={ly}; BoxSize is set to {lx}", file=sys.stderr)
    boxsize = float(lx)

    # sanity
    assert n == len(m) == len(u) == len(vx), "ragged IC arrays"
    if not (np.all(np.isfinite(coords)) and np.all(coords >= -1e-12) and np.all(coords <= boxsize + 1e-12)):
        raise SystemExit("coordinates outside [0, BoxSize] after shift")
    if not (np.all(m > 0) and np.all(np.isfinite(m))):
        raise SystemExit("non-positive or non-finite masses")
    if not (np.all(u > 0) and np.all(np.isfinite(u))):
        raise SystemExit("non-positive or non-finite internal energies")

    imat = None
    if mat is not None:
        unique = np.unique(mat)
        imat = np.searchsorted(unique, mat).astype(np.uint32)
        print("material map (imat <- materialId):", {int(i): int(v) for i, v in enumerate(unique)})

    with h5py.File(args.output, "w") as out:
        header = out.create_group("Header")
        header.attrs["NumPart_ThisFile"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.int32)
        header.attrs["NumPart_Total"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.uint32)
        header.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
        header.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
        header.attrs["Time"] = 0.0
        header.attrs["Redshift"] = 0.0
        header.attrs["BoxSize"] = boxsize
        header.attrs["NumFilesPerSnapshot"] = 1
        header.attrs["Omega0"] = 0.0
        header.attrs["OmegaLambda"] = 0.0
        header.attrs["HubbleParam"] = 1.0
        for flag in ("Flag_Sfr", "Flag_Cooling", "Flag_StellarAge", "Flag_Metals",
                     "Flag_Feedback", "Flag_DoublePrecision", "Flag_IC_Info"):
            header.attrs[flag] = 0
        header.attrs["Lz_slab_code"] = lz  # informational (thin z extent for KHI)

        p = out.create_group("PartType0")
        p.create_dataset("Coordinates", data=coords, dtype=np.float64)
        p.create_dataset("Velocities", data=np.column_stack((vx, vy, vz)), dtype=np.float64)
        p.create_dataset("Masses", data=m, dtype=np.float64)
        p.create_dataset("InternalEnergy", data=u, dtype=np.float64)
        p.create_dataset("ParticleIDs", data=ids, dtype=np.uint64)
        if mat is not None:
            p.create_dataset("Materials", data=imat, dtype=np.uint32)
            p.create_dataset("EosMaterialId", data=mat, dtype=np.uint32)
        if h is not None:
            p.create_dataset("SmoothingLength", data=h, dtype=np.float64)

    print(f"wrote {args.output}: {n} particles, BoxSize={boxsize}, slab z=[0,{lz:.6g}]"
          + (", with Materials" if mat is not None else ", no Materials (ideal gas)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
