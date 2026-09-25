#!/usr/bin/env python3
"""Square-test metrics + density maps for planetg (GIZMO) REMIX square runs.

Metrics (REMIX-style, mirrors remix_tests/square_metrics.py):
- particle speeds at t=3 (perfect equilibrium -> ~0)
- squareness S4 = <cos(4 theta)> of dense particles (rho > 2.5), relative to
  the square center (box center); square -> positive, circular blob -> ~0.

Usage: python square_metrics_planetg.py RUN_DIR OUT_PNG
"""
import sys, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

run_dir = sys.argv[1]
out_png = sys.argv[2]

def load(path):
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        d = dict(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2],
                 v=g["Velocities"][:], rho=g["Density"][:],
                 t=float(h["Header"].attrs.get("Time", np.nan)),
                 box=float(h["Header"].attrs.get("BoxSize", 1.0)))
    return d

def squareness(x, y):
    th = np.arctan2(y, x)
    return float(np.mean(np.cos(4.0 * th)))

snaps = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.hdf5")))
if not snaps:
    raise SystemExit(f"no snapshots in {run_dir}")
print(f"{run_dir}: {len(snaps)} snapshots")

rows = []
frames = []
s0 = None
for sp in snaps:
    d = load(sp)
    cx = cy = 0.5 * d["box"]
    x, y = d["x"] - cx, d["y"] - cy
    v = np.sqrt((d["v"] ** 2).sum(axis=1))
    inner = d["rho"] > 2.5
    s4 = squareness(x[inner], y[inner])
    if s0 is None:
        s0 = s4
    rows.append((d["t"], v.mean(), v.max(), s4, s4 / s0))
    frames.append(d)
    print(f"  t={d['t']:5.2f}  mean|v|={v.mean():.4f}  max|v|={v.max():.4f}  "
          f"S4={s4:+.4f}  S4/S4(0)={s4 / s0:.3f}  N={len(x)}")

# density maps, sphexa REMIX style: |z|<0.05 slice, RdYlBu_r, square markers,
# centered coords, no ticks, one column per snapshot (t = 0, 1.5, 3.0)
# (see D:\sphexa-geos\PLOT_STYLE.md)
ncol = len(frames)
fig, axes = plt.subplots(1, ncol, figsize=(2.9 * ncol, 3.2),
                         sharex=True, sharey=True, squeeze=False)
axes = axes[0]
for ax, d in zip(axes, frames):
    cx = cy = 0.5 * d["box"]
    m = np.abs(d["z"] - 0.5 * d["box"]) < 0.05
    ax.scatter(d["x"][m] - cx, d["y"][m] - cy, c=d["rho"][m], s=3,
               cmap="RdYlBu_r", vmin=0.9, vmax=4.2, marker="s", linewidths=0)
    ax.set_title(f"t = {d['t']:.2f}", fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-0.52, 0.52); ax.set_ylim(-0.52, 0.52)
    ax.set_xticks([]); ax.set_yticks([])
axes[0].set_ylabel(run_dir.rstrip("/").split("/")[-1], fontsize=10)
name = run_dir.rstrip("/").split("/")[-1]
t_end, vmean, vmax, s4, s4rel = rows[-1]
fig.suptitle(f"planetg square test: {name}   "
             f"t={t_end:g}: mean|v|={vmean:.4f}, max|v|={vmax:.4f}, "
             f"S4/S4(0)={s4rel:.3f}", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(out_png, dpi=150)
print("saved", out_png)
