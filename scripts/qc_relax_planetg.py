#!/usr/bin/env python3
"""Relaxation QC for planetg (GIZMO) runs — mirrors qc_relax_091_011.py style.

Radial scatter profiles (rho, T, s, p) for every available snapshot, plus a
|v| panel and shell-median density drift relative to snapshot_000.
GIZMO snapshot layout: PartType0/{Coordinates,Velocities,Density,Temperature,
Entropy,Pressure,Masses,Materials,...}.

Units (from the run log): UnitLength = 1 R_earth, UnitVelocity = 1 km/s,
UnitDensity = 0.36838 g/cm^3, so P_code * 0.36838 = GPa and entropy code * 1e10
= erg g^-1 K^-1.

Usage: python qc_relax_planetg.py --run-dir PATH --out PNG [--title STR]
"""
import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RHO_UNIT_CGS = 0.36838          # g cm^-3 per code density
P_UNIT_GPA = 0.36838            # GPa per code pressure
S_UNIT = 1.0e10                 # erg g^-1 K^-1 per code entropy
MAXPTS = 250_000

rng = np.random.default_rng(0)
colors = ["#222222", "#457b9d", "#2a9d8f", "#e76f51", "#b5838d", "#6d597a"]


def load(path):
    f = h5py.File(path, "r")
    g = f["PartType0"]
    d = {k: g[k][:] for k in
         ("Coordinates", "Velocities", "Density", "Temperature", "Entropy",
          "Pressure", "Masses", "Materials")}
    t = float(f["Header"].attrs.get("Time", np.nan))
    f.close()
    r = np.sqrt((d["Coordinates"] ** 2).sum(axis=1))
    vmag = np.sqrt((d["Velocities"] ** 2).sum(axis=1))
    return dict(r=r, rho=d["Density"] * RHO_UNIT_CGS, T=d["Temperature"],
                s=d["Entropy"] * S_UNIT, p=d["Pressure"] * P_UNIT_GPA,
                v=vmag, m=d["Masses"], imat=d["Materials"], t=t,
                pos=d["Coordinates"], vel=d["Velocities"])


def shell_medians(r, rho, nshell=40):
    r_max = np.percentile(r, 99.5)
    edges = np.linspace(0, r_max, nshell + 1)
    shell = np.digitize(r, edges) - 1
    med = np.array([np.median(rho[shell == i]) if np.any(shell == i) else np.nan
                    for i in range(nshell)])
    return edges, med


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    snaps = sorted(glob.glob(os.path.join(args.run_dir, "snapshot_*.hdf5")))
    if not snaps:
        raise SystemExit(f"no snapshots in {args.run_dir}")
    print(f"=== {args.title or args.run_dir}: {len(snaps)} snapshot(s)")

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axdrift = axes[1][2]
    shell_ref = None
    edges = None

    for isnap, snap in enumerate(snaps):
        d = load(snap)
        color = colors[isnap % len(colors)]
        n = len(d["r"])
        idx = rng.choice(n, min(MAXPTS, n), replace=False)
        label = f"t = {d['t']:.2f}"
        for ax, vals in ((axes[0][0], d["rho"]), (axes[0][1], d["T"]),
                         (axes[0][2], d["s"]), (axes[1][0], d["p"]),
                         (axes[1][1], d["v"])):
            ax.scatter(d["r"][idx], vals[idx], s=0.15, alpha=0.35, color=color,
                       linewidths=0, label=label, rasterized=True)

        edges, med = shell_medians(d["r"], d["rho"])
        if shell_ref is None:
            shell_ref = med
            print(f"  t={d['t']:.3f}: reference shell profile")
        else:
            interior = (edges[:-1] < 0.85 * edges[-1])
            rel = (med[interior] - shell_ref[interior]) / shell_ref[interior]
            axdrift.plot(0.5 * (edges[:-1] + edges[1:])[interior], 100 * rel,
                         color=color, label=label)
            print(f"  t={d['t']:.3f}: median |rho drift| interior = "
                  f"{np.nanmedian(np.abs(rel)) * 100:.3f}%  (max {np.nanmax(np.abs(rel)) * 100:.2f}%)")

        # global health stats
        mtot = d["m"].sum()
        com = (d["pos"] * d["m"][:, None]).sum(axis=0) / mtot
        bulkv = (d["vel"] * d["m"][:, None]).sum(axis=0) / mtot
        core = d["imat"] == 1
        print(f"  t={d['t']:.3f}: N={n} Mtot={mtot:.5f} |COM|={np.linalg.norm(com):.2e} "
              f"|bulk v|={np.linalg.norm(bulkv):.2e} km/s")
        print(f"      |v| med/95%/max = {np.median(d['v']):.2e}/"
              f"{np.percentile(d['v'], 95):.2e}/{d['v'].max():.2e} km/s")
        print(f"      rho min/max = {d['rho'].min():.3f}/{d['rho'].max():.1f} g/cm3; "
              f"p min = {d['p'].min():.3e} GPa")
        print(f"      mantle (imat=0): {np.sum(~core)} particles, median r={np.median(d['r'][~core]):.3f}; "
              f"core (imat=1): {np.sum(core)} particles, median r={np.median(d['r'][core]):.3f}")

    axes[0][0].set_ylabel(r"Density [g cm$^{-3}$]"); axes[0][0].set_title("Density")
    axes[0][0].set_ylim(bottom=0)
    axes[0][1].set_ylabel("Temperature [K]"); axes[0][1].set_title("Temperature")
    axes[0][2].set_ylabel(r"Entropy [erg g$^{-1}$ K$^{-1}$]"); axes[0][2].set_title("Specific entropy")
    axes[1][0].set_yscale("log"); axes[1][0].set_ylabel("Pressure [GPa]"); axes[1][0].set_title("Pressure")
    axes[1][1].set_ylabel(r"$|v|$ [km s$^{-1}$]"); axes[1][1].set_title("Speed (should damp out)")
    axes[1][1].set_yscale("log")
    axdrift.set_ylabel(r"$\Delta\rho/\rho_0$ [%]"); axdrift.set_title("Shell-median density drift vs t=0")
    axdrift.axhline(0, color="k", lw=0.6)
    for ax in axes.flat:
        ax.set_xlabel(r"Radius [$R_\oplus$]"); ax.grid(alpha=0.3)
    axes[0][0].legend(markerscale=20, loc="best", fontsize=9)
    fig.suptitle(f"planetg relaxation QC: {args.title}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=200)
    print("saved", args.out)


if __name__ == "__main__":
    main()
