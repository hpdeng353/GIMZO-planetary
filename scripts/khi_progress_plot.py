import glob
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU = 0.826226
RHO_MIN, RHO_MAX = 1000.0, 8000.0

def load_pg(path):
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        return dict(x=pos[:, 0] - 0.5, y=pos[:, 1] - 0.5, z=pos[:, 2],
                    rho=g["Density"][:] * 368.38,
                    t=float(h["Header"].attrs["Time"]))

def load_sx(path):
    with h5py.File(path, "r") as h:
        g = h["Step#0"]
        return dict(x=g["x"][:], y=g["y"][:], z=g["z"][:],
                    rho=g["rho"][:] * 368.38,
                    t=float(np.atleast_1d(g.attrs["time"])[0]))

rows = [("sphexa em", load_sx, sorted(glob.glob("/home/hpdeng/sphexa-geos/final_tests/khi_em/dump_khi.*"))),
        ("planetg em", load_pg, sorted(glob.glob("/home/hpdeng/planetg/runs/remix/khi_em/snapshot_*.hdf5")))]
want = [0.0, 0.5, 1.0]

fig, axes = plt.subplots(2, len(want), figsize=(3.0 * len(want), 3.2 * 2), squeeze=False)
for r, (label, loader, snaps) in enumerate(rows):
    frames = [loader(p) for p in snaps]
    z0 = float(np.median(frames[0]["z"]))
    for c, w in enumerate(want):
        ax = axes[r][c]
        d = min(frames, key=lambda f: abs(f["t"] / TAU - w))
        m = np.abs(d["z"] - z0) < 0.05
        ax.scatter(d["x"][m], d["y"][m], c=d["rho"][m], s=2, cmap="RdYlBu_r",
                   vmin=RHO_MIN, vmax=RHO_MAX, marker="s", linewidths=0)
        ax.set_xlim(-0.52, 0.52); ax.set_ylim(-0.52, 0.52)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if c == 0:
            ax.set_ylabel(label, fontsize=10)
        if r == 0:
            ax.set_title("t = %.1f " % (d["t"] / TAU) + r"$\tau_{KH}$", fontsize=11)
fig.suptitle("KHI iron/rock: planetg MFM vs sphexa (in progress)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("/home/hpdeng/planetg/runs/remix/logs/khi_progress.png", dpi=150)
print("saved khi_progress.png")
