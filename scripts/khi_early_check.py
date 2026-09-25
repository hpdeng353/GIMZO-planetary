import h5py, numpy as np
TAU = 0.826226
def amp(path):
    with h5py.File(path, "r") as h:
        g = h["PartType0"]
        pos = g["Coordinates"][:]
        x = pos[:, 0] - 0.5
        y = pos[:, 1] - 0.5
        vy = g["Velocities"][:, 1]
        imat = g["Materials"][:]
        t = float(h["Header"].attrs["Time"])
    iron = imat == 1
    edges = np.linspace(-0.5, 0.5, 65)
    env = np.full(64, np.nan)
    for b in range(64):
        sel = iron & (x >= edges[b]) & (x < edges[b + 1])
        if sel.sum() > 20:
            env[b] = np.percentile(np.abs(y[sel]), 99.0)
    env = env[np.isfinite(env)]
    a = np.sqrt(np.mean((env - 0.25) ** 2))
    print("%s t/tau=%.2f amp=%.5f rms_vy=%.5f max|vy|=%.5f"
          % (path.split("/")[-1], t / TAU, a, np.sqrt(np.mean(vy ** 2)), np.abs(vy).max()))

amp("/home/hpdeng/planetg/runs/remix/khi_em/snapshot_000.hdf5")
amp("/home/hpdeng/planetg/runs/remix/khi_em/snapshot_001.hdf5")
