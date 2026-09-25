#!/usr/bin/env python3
"""Relaxation QC for planetg, sphexa qc_relax_profile.py style:
shell-averaged rho/T/S/p profiles at t = 0, 1, 2, 3 with 1-sigma bars.

GIZMO snapshot layout: PartType0/{Coordinates,Density,Temperature,Entropy,
Pressure}. Units: UnitLength = R_earth, UnitDensity = 0.36838 g/cm^3,
P_code * 0.36838 = GPa, entropy code * 1e10 = erg g^-1 K^-1.

Usage: python qc_profile_planetg.py RUN_DIR OUT_PNG
"""
import sys, glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO_UNIT_CGS = 0.36838
P_UNIT_GPA = 0.36838
S_UNIT = 1.0e10

run_dir = sys.argv[1]
out_png = sys.argv[2]

targets = [0.0, 1.0, 2.0, 3.0]
colors = {0.0: "tab:blue", 1.0: "tab:orange", 2.0: "tab:green", 3.0: "tab:red"}

# map target times to nearest snapshots
snaps = sorted(glob.glob(os.path.join(run_dir, "snapshot_*.hdf5")))
if not snaps:
    raise SystemExit(f"no snapshots in {run_dir}")
stab = {}
for p in snaps:
    with h5py.File(p, "r") as h:
        stab[p] = float(h["Header"].attrs.get("Time", np.nan))
picked = {}
for tt in targets:
    p = min(stab, key=lambda k: abs(stab[k] - tt))
    picked[tt] = p
    print(f"t={tt:g} -> {os.path.basename(p)} (Time={stab[p]:.4f})")

data = {}
for tt, path in picked.items():
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        data[tt] = dict(
            r=np.sqrt((pos ** 2).sum(axis=1)),
            rho=g["Density"][:] * RHO_UNIT_CGS,
            temp=g["Temperature"][:],
            S=g["Entropy"][:] * S_UNIT,
            p=g["Pressure"][:] * P_UNIT_GPA,
        )
print("particles:", len(data[0.0]["r"]))

r_max = max(float(np.percentile(d["r"], 99.9)) for d in data.values())
edges = np.linspace(0.0, r_max, 81)
centers = 0.5 * (edges[1:] + edges[:-1])

def binned(r, v):
    ib = np.clip(np.digitize(r, edges) - 1, 0, len(centers) - 1)
    mean = np.full(len(centers), np.nan)
    std = np.full(len(centers), np.nan)
    for i in range(len(centers)):
        sel = ib == i
        if np.count_nonzero(sel) > 20:
            mean[i] = np.mean(v[sel])
            std[i] = np.std(v[sel])
    return mean, std

panels = [("rho", r"Density [g cm$^{-3}$]", False),
          ("temp", "Temperature [K]", False),
          ("S", r"Entropy [erg g$^{-1}$ K$^{-1}$]", False),
          ("p", "Pressure [GPa]", True)]

fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
for ax, (key, label, logy) in zip(axes.ravel(), panels):
    for tt in targets:
        mean, std = binned(data[tt]["r"], data[tt][key])
        ax.errorbar(centers, mean, yerr=std, color=colors[tt], lw=1.4,
                    errorevery=4, capsize=2, elinewidth=0.7, label=f"t = {tt:g}")
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)
    if key == "rho":
        ax.set_ylim(bottom=0)
    if key == "S":
        ax.ticklabel_format(axis="y", scilimits=(0, 0))
    if logy:
        ax.set_yscale("log")
axes[0, 0].legend(loc="best", fontsize=10)
for ax in axes[1]:
    ax.set_xlabel(r"Radius [$R_\oplus$]")
name = run_dir.rstrip("/").split("/")[-1]
fig.suptitle(f"planetg relaxation QC: {name}  (shell mean +/- 1 sigma)")
fig.tight_layout()
fig.savefig(out_png, dpi=140)
print("saved", out_png)

# numeric summary: drift at t=1,2,3 vs t=0
for key in ("rho", "temp", "S", "p"):
    m0, _ = binned(data[0.0]["r"], data[0.0][key])
    for tt in (1.0, 2.0, 3.0):
        m1, _ = binned(data[tt]["r"], data[tt][key])
        good = np.isfinite(m0) & np.isfinite(m1) & (np.abs(m0) > 0)
        drift = np.abs(m1[good] - m0[good]) / np.abs(m0[good])
        print(f"{key}: median |drift| t={tt:g} vs t=0 = {np.median(drift)*100:.2f}%"
              f"  max = {np.max(drift)*100:.1f}%")
print("done")
