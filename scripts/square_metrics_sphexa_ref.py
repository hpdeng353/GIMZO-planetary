#!/usr/bin/env python3
"""Square-test metrics for sphexa dumps (final_tests/square_em, square_es),
same definitions as planetg/scripts/square_metrics_planetg.py."""
import sys, glob, os
import numpy as np
import h5py

run_dir = sys.argv[1]

def load(path):
    with h5py.File(path, "r") as h:
        g = h["Step#0"] if "Step#" in list(h.keys())[0] or "Step#0" in h else h[list(h.keys())[0]]
        d = {k: g[k][:] for k in ("x", "y", "z", "vx", "vy", "vz", "rho")}
        try:
            t = float(np.atleast_1d(g.attrs["time"])[0])
        except Exception:
            t = float("nan")
        d["t"] = t
    return d

def squareness(x, y):
    th = np.arctan2(y, x)
    return float(np.mean(np.cos(4.0 * th)))

snaps = sorted(glob.glob(os.path.join(run_dir, "dump_square.*")))
snaps = [s for s in snaps if not s.endswith(".h5") or True]
print(f"{run_dir}: {len(snaps)} snapshots")
s0 = None
for sp in snaps:
    d = load(sp)
    cx = 0.5 * (d["x"].min() + d["x"].max())
    cy = 0.5 * (d["y"].min() + d["y"].max())
    v = np.sqrt(d["vx"] ** 2 + d["vy"] ** 2 + d["vz"] ** 2)
    inner = d["rho"] > 2.5
    s4 = squareness(d["x"][inner] - cx, d["y"][inner] - cy)
    if s0 is None:
        s0 = s4
    print(f"  {os.path.basename(sp)}  t={d['t']:5.2f}  mean|v|={v.mean():.4f}  "
          f"max|v|={v.max():.4f}  S4={s4:+.4f}  S4/S4(0)={s4/s0:.3f}  N={len(d['x'])}")
