#!/usr/bin/env python3
"""RTI (REMIX Fig. 12) comparison: planetg (MFM) vs sphexa reference.

Metric (mirrors remix_tests/rti_metrics.py): spike = deepest iron penetration
below the interface (1st pct of y_rel / 0.05), bubble = highest rock above it
(99th pct / 0.05); also rms/max |vy|.

Outputs:
  rti_compare.png - density slice panels, rows = runs, cols = t/tau_RT
  rti_growth.png  - spike/bubble growth curves

Run on the cluster (both data sets live there).
"""
import glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU_RT = 0.014203
LY_HALF = 0.05
RHO_MIN, RHO_MAX = 1000.0, 8000.0

ROWS = [("sphexa em", "/home/hpdeng/sphexa-geos/final_tests/rti_em", "sphexa"),
        ("planetg em", "/home/hpdeng/planetg/runs/remix/rti_em", "planetg"),
        ("sphexa es", "/home/hpdeng/sphexa-geos/final_tests/rti_es", "sphexa"),
        ("planetg es", "/home/hpdeng/planetg/runs/remix/rti_es", "planetg")]


def load_pg(path):
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        return dict(x=pos[:, 0] - 0.025, y=pos[:, 1] - 0.05, z=pos[:, 2],
                    vy=g["Velocities"][:, 1], rho=g["Density"][:] * 368.38,
                    imat=g["Materials"][:],
                    t=float(h["Header"].attrs["Time"]))


def load_sx(path):
    with h5py.File(path, "r") as h:
        g = h["Step#0"]
        return dict(x=g["x"][:], y=g["y"][:], z=g["z"][:], vy=g["vy"][:],
                    rho=g["rho"][:] * 368.38, imat=g["materialId"][:],
                    t=float(np.atleast_1d(g.attrs["time"])[0]))


def snaps_pg(d):
    return sorted(glob.glob(os.path.join(d, "snapshot_*.hdf5")))


def snaps_sx(d):
    return sorted(glob.glob(os.path.join(d, "dump_rti.*")))


print(f"{'run':<11} {'t/tau':>6} {'spike':>8} {'bubble':>8} {'rms_vy':>9} {'max|vy|':>9}")
allrows = {}
for label, dpath, kind in ROWS:
    snaps = snaps_pg(dpath) if kind == "planetg" else snaps_sx(dpath)
    if not snaps:
        print(f"{label}: NO SNAPSHOTS in {dpath}")
        continue
    rows = []
    z0 = None
    for sp in snaps:
        d = load_pg(sp) if kind == "planetg" else load_sx(sp)
        if z0 is None:
            z0 = float(np.median(d["z"]))
        d["z"] = d["z"] - z0
        iron = d["imat"] == (1 if kind == "planetg" else 63)
        rock = ~iron
        spike = float(-np.percentile(d["y"][iron], 1.0) / LY_HALF)
        bubble = float(np.percentile(d["y"][rock], 99.0) / LY_HALF)
        rms = float(np.sqrt(np.mean(d["vy"] ** 2)))
        mx = float(np.abs(d["vy"]).max())
        rows.append((d["t"], spike, bubble, rms, mx, d))
        print(f"{label:<11} {d['t'] / TAU_RT:6.2f} {spike:8.4f} {bubble:8.4f} {rms:9.5f} {mx:9.5f}")
    allrows[label] = rows
    print()

# ---- panels: rows = runs, cols = t/tau = 0..5 ----
want = [0, 1, 2, 3, 4, 5]
fig, axes = plt.subplots(4, len(want), figsize=(2.0 * len(want), 2.9 * 4),
                         squeeze=False)
for r, (label, dpath, kind) in enumerate(ROWS):
    if label not in allrows:
        continue
    rows = allrows[label]
    for c, w in enumerate(want):
        ax = axes[r][c]
        t, spike, bubble, rms, mx, d = min(rows, key=lambda rr: abs(rr[0] / TAU_RT - w))
        m = np.abs(d["z"]) < 0.35 * (d["z"].max() - d["z"].min() + 1e-9)
        ax.scatter(d["x"][m], d["y"][m], c=d["rho"][m], s=1.5, cmap="RdYlBu_r",
                   vmin=RHO_MIN, vmax=RHO_MAX, marker="s", linewidths=0)
        ax.set_xlim(-0.027, 0.027); ax.set_ylim(-0.053, 0.053)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        if c == 0:
            ax.set_ylabel(label, fontsize=10)
        if r == 0:
            ax.set_title(f"{t / TAU_RT:.0f}" + r" $\tau_{RT}$", fontsize=10)
fig.suptitle("RTI iron/rock (REMIX Fig. 12): planetg MFM vs sphexa", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("/home/hpdeng/planetg/runs/remix/logs/rti_compare.png", dpi=150)
print("saved rti_compare.png")

# ---- growth curves ----
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
sty = {"sphexa em": ("o-", "tab:blue"), "planetg em": ("s--", "tab:cyan"),
       "sphexa es": ("o-", "tab:red"), "planetg es": ("s--", "tab:orange")}
for label, rows in allrows.items():
    tt = [r[0] / TAU_RT for r in rows]
    axes[0].plot(tt, [r[1] for r in rows], sty[label][0], color=sty[label][1],
                 label=label, ms=3.5, lw=1.3)
    axes[1].plot(tt, [r[2] for r in rows], sty[label][0], color=sty[label][1],
                 label=label, ms=3.5, lw=1.3)
axes[0].set_title("iron spike depth / (Ly/2)")
axes[1].set_title("rock bubble height / (Ly/2)")
for ax in axes:
    ax.set_xlabel(r"$t / \tau_{RT}$"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
fig.suptitle("RTI growth: planetg vs sphexa")
fig.tight_layout()
fig.savefig("/home/hpdeng/planetg/runs/remix/logs/rti_growth.png", dpi=150)
print("saved rti_growth.png")
