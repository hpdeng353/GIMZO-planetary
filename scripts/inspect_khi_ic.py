import h5py, numpy as np

print("=== planetg khi_em snapshot_000 ===")
with h5py.File("/home/hpdeng/planetg/runs/remix/khi_em/snapshot_000.hdf5") as h:
    g = h["PartType0"]
    print("keys:", list(g.keys()))
    pos = g["Coordinates"][:]
    print("x range", pos[:, 0].min(), pos[:, 0].max())
    print("y range", pos[:, 1].min(), pos[:, 1].max())
    print("z range", pos[:, 2].min(), pos[:, 2].max())
    print("N", len(pos))
    rho = g["Density"][:]; P = g["Pressure"][:]
    u = g["InternalEnergy"][:]; T = g["Temperature"][:]
    imat = g["Materials"][:]
    for m in np.unique(imat):
        s = imat == m
        print("imat=%d: N=%d rho med=%.4f  P med=%.4f  u med=%.5f  T med=%.1f"
              % (m, s.sum(), np.median(rho[s]), np.median(P[s]),
                 np.median(u[s]), np.median(T[s])))

print("=== sphexa dump_khi.0000 ===")
with h5py.File("/home/hpdeng/sphexa-geos/final_tests/khi_em/dump_khi.0000") as h:
    g = h["Step#0"]
    print("keys:", list(g.keys()))
    print("x range", g["x"][:].min(), g["x"][:].max(),
          " z range", g["z"][:].min(), g["z"][:].max())
    rho = g["rho"][:]
    u = g["u"][:]
    p = g["p"][:] if "p" in g else None
    mat = g["materialId"][:]
    for m in np.unique(mat):
        s = mat == m
        line = "mat=%d: N=%d rho med=%.4f  u med=%.5f" % (m, s.sum(),
               np.median(rho[s]), np.median(u[s]))
        if p is not None:
            line += "  p med=%.4f" % np.median(p[s])
        print(line)
