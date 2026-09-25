#!/usr/bin/env python3
"""KHI (REMIX Fig. 9) comparison: planetg (MFM) vs sphexa reference.

Metric (mirrors remix_tests/khi_metrics.py): bin iron particles (mat 63 ->
planetg imat=1) in x, envelope = 99th pct of |y - y0|, amplitude = RMS of
(envelope - 0.25). Also rms/max |vy|.

Outputs:
  khi_compare.png  - density slice panels, rows = runs, cols = t/tau = 0, 1, 2
  khi_growth.png   - billow amplitude growth curves

Run on the cluster (both data sets live there).
"""
import glob, os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU_KH = 0.826226
NBINS = 64
RHO_MIN, RHO_MAX = 1000.0, 8000.0   # kg/m^3-ish display range (code density * 0.36838 g/cm3 -> *368.38)

PG = {"planetg em": "/home/hpdeng/planetg/runs/remix/khi_em",
      "planetg es": "/home/hpdeng/planetg/runs/remix/khi_es"}
SX = {"sphexa em": "/home/hpdeng/sphexa-geos/final_tests/khi_em",
      "sphexa es": "/home/hpdeng/sphexa-geos/final_tests/khi_es"}
# display rows: sphexa/planetg pairs per IC type
ROWS = [("em", SX["sphexa em"], "sphexa em", "sphexa"),
        ("em", PG["planetg em"], "planetg em", "planetg"),
        ("es", SX["sphexa es"], "sphexa es", "sphexa"),
        ("es", PG["planetg es"], "planetg es", "planetg")]


def load_pg(path):
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        return dict(x=pos[:, 0] - 0.5, y=pos[:, 1] - 0.5, z=pos[:, 2],
                    vy=g["Velocities"][:, 1], rho=g["Density"][:] * 368.38,
                    imat=g["Materials"][:],
                    t=float(h["Header"].attrs["Time"]))


def load_sx(path):
    with h5py.File(path, "r") as h:
        g = h["Step#0"]
        return dict(x=g["x"][:], y=g["y"][:], z=g["z"][:], vy=g["vy"][:],
                    rho=g["rho"][:] * 368.38, imat=g["materialId"][:],
                    t=float(np.atleast_1d(g.attrs["time"])[0]))


def envelope_amplitude(x, y):
    edges = np.linspace(-0.5, 0.5, NBINS + 1)
    env = np.full(NBINS, np.nan)
    for b in range(NBINS):
        sel = (x >= edges[b]) & (x < edges[b + 1])
        if np.count_nonzero(sel) > 20:
            env[b] = np.percentile(np.abs(y[sel]), 99.0)
    env = env[np.isfinite(env)]
    return np.sqrt(np.mean((env - 0.25) ** 2))


def snaps_pg(d):
    return sorted(glob.glob(os.path.join(d, "snapshot_*.hdf5")))


def snaps_sx(d):
    return sorted(glob.glob(os.path.join(d, "dump_khi.*")))


# ---- metrics table ----
print(f"{'run':<11} {'t/tau':>6} {'amplitude':>10} {'rms_vy':>9} {'max|vy|':>9}")
allrows = {}
for key, dpath, label, kind in ROWS:
    snaps = snaps_pg(dpath) if kind == "planetg" else snaps_sx(dpath)
    rows = []
    z0 = None
    for sp in snaps:
        d = load_pg(sp) if kind == "planetg" else load_sx(sp)
        if z0 is None:
            z0 = float(np.median(d["z"]))  # slab centre from first snapshot
        d["z"] = d["z"] - z0
        iron = d["imat"] == (1 if kind == "planetg" else 63)
        amp = envelope_amplitude(d["x"][iron], d["y"][iron])
        rms = float(np.sqrt(np.mean(d["vy"] ** 2)))
        mx = float(np.abs(d["vy"]).max())
        rows.append((d["t"], amp, rms, mx, d))
        print(f"{label:<11} {d['t'] / TAU_KH:6.2f} {amp:10.5f} {rms:9.5f} {mx:9.5f}")
    allrows[label] = rows
    print()

# ---- panel figure: rows = runs, cols = t/tau = 0, 1, 2 ----
want = [0.0, 1.0, 2.0]
fig, axes = plt.subplots(4, 3, figsize=(9.5, 3.1 * 4), squeeze=False)
for r, (ic, dpath, label, kind) in enumerate(ROWS):
    rows = allrows[label]
    for c, w in enumerate(want):
        ax = axes[r][c]
        t, amp, rms, mx, d = min(rows, key=lambda rr: abs(rr[0] / TAU_KH - w))
        m = np.abs(d["z"]) < 0.05
        ax.scatter(d["x"][m], d["y"][m], c=d["rho"][m], s=2, cmap="RdYlBu_r",
                   vmin=RHO_MIN, vmax=RHO_MAX, marker="s", linewidths=0)
        ax.set_xlim(-0.52, 0.52); ax.set_ylim(-0.52, 0.52)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        if c == 0:
            ax.set_ylabel(label, fontsize=10)
        if r == 0:
            ax.set_title(f"t = {t / TAU_KH:.1f} " + r"$\tau_{KH}$", fontsize=11)
fig.suptitle("KHI iron/rock (REMIX Fig. 9): planetg MFM vs sphexa", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig("/home/hpdeng/planetg/runs/remix/logs/khi_compare.png", dpi=150)
print("saved khi_compare.png")

# ---- growth curves ----
fig, ax = plt.subplots(figsize=(7, 5))
sty = {"sphexa em": ("o-", "tab:blue"), "planetg em": ("s--", "tab:cyan"),
       "sphexa es": ("o-", "tab:red"), "planetg es": ("s--", "tab:orange")}
for label, rows in allrows.items():
    tt = [r[0] / TAU_KH for r in rows]
    aa = [r[1] for r in rows]
    ax.plot(tt, aa, sty[label][0], color=sty[label][1], label=label, ms=4, lw=1.3)
ax.set_xlabel(r"$t / \tau_{KH}$"); ax.set_ylabel("billow amplitude [code]")
ax.grid(alpha=0.3); ax.legend(fontsize=9)
ax.set_title("KHI growth: planetg vs sphexa")
fig.tight_layout()
fig.savefig("/home/hpdeng/planetg/runs/remix/logs/khi_growth.png", dpi=150)
print("saved khi_growth.png")
