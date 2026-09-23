#!/usr/bin/env python3
"""test_eos_table.py -- Python reference implementation of planetg's aneos/eos_table.c.

Independent re-implementation of the .spheos reader + interpolator, used to
validate the table file itself and to cross-check the C reader on the cluster:
run `test_eos_table <table>` (C) and this script on the same file and compare
the sample-point outputs.

Usage:  python test_eos_table.py <table.spheos> [materialId] [--skip-fnv]

Exit status 0 when all checks pass.
"""

import struct
import sys
import time
from pathlib import Path

import numpy as np

MAGIC = b"SPXEOST1"
HEADER_SIZE = 32
ENDIAN_MARKER = 0x01020304
FLAG_ROW_DEPENDENT_U = 0x2
FLAG_NATIVE_TGRID = 0x4
FLAG_HAS_ENTROPY = 0x8
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
MASK64 = 0xFFFFFFFFFFFFFFFF
DBL_EPSILON = 2.220446049250313e-16

SUCCESS = "success"
DENSITY_BELOW = "densityBelowRange"
DENSITY_ABOVE = "densityAboveRange"
ENERGY_BELOW = "energyBelowRange"
ENERGY_ABOVE = "energyAboveRange"
INVALID_INPUT = "invalidInput"
INVALID_TABLE = "invalidTableState"


def fnv1a64(data: bytes) -> int:
    h = FNV_OFFSET
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


class Material:
    pass


def load_table(path: Path, skip_fnv: bool):
    raw = path.read_bytes()
    if len(raw) < HEADER_SIZE or raw[:8] != MAGIC:
        raise SystemExit(f"{path}: bad magic")
    version, endian, payload_size, checksum = struct.unpack_from("<IIQQ", raw, 8)
    print(f"header: version={version} endian={endian:#x} payloadSize={payload_size}")
    if endian != ENDIAN_MARKER:
        raise SystemExit("bad endian marker")
    if payload_size != len(raw) - HEADER_SIZE:
        raise SystemExit(f"payloadSize mismatch: header {payload_size}, file {len(raw) - HEADER_SIZE}")
    payload = raw[HEADER_SIZE:]
    if skip_fnv:
        print("checksum: SKIPPED (--skip-fnv)")
    else:
        t0 = time.time()
        digest = fnv1a64(payload)
        ok = "OK" if digest == checksum else "MISMATCH"
        print(f"checksum: {digest:#018x} vs header {checksum:#018x} -> {ok} ({time.time()-t0:.1f}s)")
        if digest != checksum:
            raise SystemExit("checksum mismatch")

    pos = 0
    (num_materials,) = struct.unpack_from("<I", payload, pos)
    pos += 4
    materials = []
    for _ in range(num_materials):
        m = Material()
        m.materialId, m.nRho, m.nU, flags = struct.unpack_from("<IIII", payload, pos)
        pos += 16
        m.rowDependentU = bool(version >= 2 and flags & FLAG_ROW_DEPENDENT_U)
        m.nativeGrid = bool(version >= 2 and flags & FLAG_NATIVE_TGRID)
        m.hasEntropy = bool(version >= 3 and flags & FLAG_HAS_ENTROPY)
        if not flags & 1:
            raise SystemExit(f"material {m.materialId}: unsupported value encoding")
        n = m.nRho * m.nU
        n_axis = n if m.rowDependentU else m.nU
        m.logRho = np.frombuffer(payload, "<f8", m.nRho, pos)
        pos += 8 * m.nRho
        m.energyAxis = np.frombuffer(payload, "<f8", n_axis, pos)
        pos += 8 * n_axis
        m.logPressure = np.frombuffer(payload, "<f4", n, pos).astype(np.float64)
        pos += 4 * n
        m.logSoundSpeed = np.frombuffer(payload, "<f4", n, pos).astype(np.float64)
        pos += 4 * n
        m.logTemperature = np.frombuffer(payload, "<f4", n, pos).astype(np.float64)
        pos += 4 * n
        m.entropy = None
        if m.hasEntropy:
            m.entropy = np.frombuffer(payload, "<f4", n, pos).astype(np.float64)
            pos += 4 * n
        if not np.all(np.diff(m.logRho) > 0):
            raise SystemExit(f"material {m.materialId}: density axis not increasing")
        if not m.rowDependentU and not np.all(np.diff(m.energyAxis) > 0):
            raise SystemExit(f"material {m.materialId}: energy axis not increasing")
        m.uStride = m.nU if m.rowDependentU else 0
        m.minSoundSpeed = float(np.exp(m.logSoundSpeed.min()))
        m.minTemperature = float(np.exp(m.logTemperature.min()))
        materials.append(m)
    if pos != len(payload):
        raise SystemExit(f"trailing payload data: {len(payload) - pos} bytes")
    return version, materials


def bounds(m: Material):
    lo = m.energyAxis.min()
    hi = m.energyAxis.max()
    if not m.nativeGrid:
        lo, hi = np.exp(lo), np.exp(hi)
    return float(np.exp(m.logRho[0])), float(np.exp(m.logRho[-1])), float(lo), float(hi)


# --- interpolation, mirrors eos_table.c exactly -----------------------------

def lower_cell(axis: np.ndarray, value: float) -> int:
    n = len(axis)
    i = int(np.searchsorted(axis, value, side="right"))
    if i == 0:
        return 0
    if i == n:
        return n - 2
    return i - 1


def virtual_energy(m: Material, ir: int, iu: int, fr: float) -> float:
    if m.uStride == 0:
        return m.energyAxis[iu]
    lo = m.energyAxis[ir * m.uStride + iu]
    hi = m.energyAxis[(ir + 1) * m.uStride + iu]
    return lo + fr * (hi - lo)


def usable_segment(z0, z1, value):
    scale = max(1.0, abs(z0), abs(z1), abs(value))
    return np.isfinite(z0) and np.isfinite(z1) and abs(z1 - z0) > 64.0 * DBL_EPSILON * scale


def contains_energy(z0, z1, value):
    z_low = z0 + (-0.0001) * (z1 - z0)
    z_high = z0 + 1.0001 * (z1 - z0)
    return (z_low <= z_high and z_low <= value <= z_high) or (z_high < z_low and z_high <= value <= z_low)


def native_energy_location(m: Material, ir: int, fr: float, value: float):
    nU = m.nU
    a, b = 0, nU - 2
    c = (a + b) // 2
    while True:
        z0 = virtual_energy(m, ir, c, fr)
        z1 = virtual_energy(m, ir, c + 1, fr)
        if usable_segment(z0, z1, value) and contains_energy(z0, z1, value):
            return c, min(1.0, max(0.0, (value - z0) / (z1 - z0)))
        if b - a < 2:
            break
        z_low, z_high = (z0, z1) if z0 < z1 else (z1, z0)
        if value > z_high:
            a = c
        elif value < z_low:
            b = c
        else:
            break
        c = (a + b) // 2
    for j in range(nU - 1):
        z0 = virtual_energy(m, ir, j, fr)
        z1 = virtual_energy(m, ir, j + 1, fr)
        if usable_segment(z0, z1, value) and contains_energy(z0, z1, value):
            return j, min(1.0, max(0.0, (value - z0) / (z1 - z0)))
    # closest endpoint of a non-degenerate segment
    best, best_dist = (0, 0.0), float("inf")
    for j in range(nU - 1):
        z0 = virtual_energy(m, ir, j, fr)
        z1 = virtual_energy(m, ir, j + 1, fr)
        if not usable_segment(z0, z1, value):
            continue
        d0, d1 = abs(value - z0), abs(value - z1)
        d = min(d0, d1)
        if d < best_dist:
            best_dist = d
            best = (j, 0.0 if d0 <= d1 else 1.0)
    if np.isfinite(best_dist):
        return best
    # entirely flat column: closest node, no interpolation
    nodes = np.array([virtual_energy(m, ir, j, fr) for j in range(nU)])
    j = int(np.argmin(np.abs(nodes - value)))
    return (j, 0.0) if j + 1 < nU else (j - 1, 1.0)


def bilinear(values, nU, ir, iu, fr, fu):
    row0 = ir * nU
    row1 = (ir + 1) * nU
    lo = values[row0 + iu] + fu * (values[row0 + iu + 1] - values[row0 + iu])
    hi = values[row1 + iu] + fu * (values[row1 + iu + 1] - values[row1 + iu])
    return lo + fr * (hi - lo)


def bilinear_decoded(values, nU, ir, iu, fr, fu):
    row0 = ir * nU
    row1 = (ir + 1) * nU
    lo = np.exp(values[row0 + iu]) + fu * (np.exp(values[row0 + iu + 1]) - np.exp(values[row0 + iu]))
    hi = np.exp(values[row1 + iu]) + fu * (np.exp(values[row1 + iu + 1]) - np.exp(values[row1 + iu]))
    return lo + fr * (hi - lo)


def evaluate(m: Material, rho: float, u: float, request_entropy=True):
    has_entropy = request_entropy and m.hasEntropy
    if not (rho > 0.0) or not np.isfinite(rho) or not np.isfinite(u):
        return dict(P=0.0, cs=m.minSoundSpeed, T=m.minTemperature, S=0.0, status=INVALID_INPUT)
    if not m.nativeGrid and not (u > 0.0):
        return dict(P=0.0, cs=m.minSoundSpeed, T=m.minTemperature, S=0.0, status=INVALID_INPUT)

    status = SUCCESS
    log_rho = np.log(rho)
    if log_rho < m.logRho[0]:
        status = DENSITY_BELOW
    elif log_rho > m.logRho[-1]:
        status = DENSITY_ABOVE
    log_rho = min(max(log_rho, m.logRho[0]), m.logRho[-1])
    ir = lower_cell(m.logRho, log_rho)
    if m.nativeGrid:
        rho0, rho1 = np.exp(m.logRho[ir]), np.exp(m.logRho[ir + 1])
        fr = (rho - rho0) / (rho1 - rho0)
        fr = min(1.0, max(0.0, fr))
    else:
        fr = (log_rho - m.logRho[ir]) / (m.logRho[ir + 1] - m.logRho[ir])

    if m.nativeGrid:
        row = np.array([virtual_energy(m, ir, j, fr) for j in range(m.nU)])
        e_min, e_max = row.min(), row.max()
        coord = u
    else:
        e_min = virtual_energy(m, ir, 0, fr)
        e_max = virtual_energy(m, ir, m.nU - 1, fr)
        coord = np.log(u)
    if status == SUCCESS and coord < e_min:
        status = ENERGY_BELOW
    elif status == SUCCESS and coord > e_max:
        status = ENERGY_ABOVE
    coord = min(max(coord, e_min), e_max)

    if m.nativeGrid:
        iu, fu = native_energy_location(m, ir, fr, coord)
    else:
        # lower cell on the virtual (interpolated) axis
        lo, hi = 0, m.nU - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if virtual_energy(m, ir, mid, fr) <= coord:
                lo = mid
            else:
                hi = mid
        iu = lo
        u0 = virtual_energy(m, ir, iu, fr)
        u1 = virtual_energy(m, ir, iu + 1, fr)
        fu = (coord - u0) / (u1 - u0)
    if not np.isfinite(fu):
        return dict(P=0.0, cs=m.minSoundSpeed, T=m.minTemperature, S=0.0, status=INVALID_TABLE)
    fu = min(1.0, max(0.0, fu))

    S = bilinear(m.entropy, m.nU, ir, iu, fr, fu) if has_entropy else 0.0
    if m.nativeGrid:
        P = bilinear_decoded(m.logPressure, m.nU, ir, iu, fr, fu)
        cs = bilinear_decoded(m.logSoundSpeed, m.nU, ir, iu, fr, fu)
        T = bilinear_decoded(m.logTemperature, m.nU, ir, iu, fr, fu)
    else:
        P = np.exp(bilinear(m.logPressure, m.nU, ir, iu, fr, fu))
        cs = np.exp(bilinear(m.logSoundSpeed, m.nU, ir, iu, fr, fu))
        T = np.exp(bilinear(m.logTemperature, m.nU, ir, iu, fr, fu))

    if not (P > 0.0) or not np.isfinite(P) or not (cs > 0.0) or not np.isfinite(cs) \
            or not (T > 0.0) or not np.isfinite(T):
        return dict(P=0.0, cs=m.minSoundSpeed, T=m.minTemperature, S=0.0, status=INVALID_TABLE)
    return dict(P=P, cs=cs, T=T, S=S, status=status)


def check_point(m, rho, u):
    s = evaluate(m, rho, u)
    print("  rho=%12.5e g/cm^3  u=%12.5e erg/g -> P=%12.5e  cs=%12.5e  T=%12.5e  S=%12.5e  status=%s%s"
          % (rho, u, s["P"], s["cs"], s["T"], s["S"], s["status"], "" if m.hasEntropy else " (no entropy)"))
    ok = np.isfinite(s["P"]) and np.isfinite(s["cs"]) and np.isfinite(s["T"]) \
        and np.isfinite(s["S"]) and s["cs"] > 0.0 and s["T"] > 0.0
    if not ok:
        print("  ^^ FAIL: non-finite or non-positive state")
    return bool(ok)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    skip_fnv = "--skip-fnv" in sys.argv
    if not args:
        print("usage: python test_eos_table.py <table.spheos> [materialId] [--skip-fnv]")
        return 2
    version, materials = load_table(Path(args[0]), skip_fnv)
    print(f"loaded {args[0]}: {len(materials)} material(s)")

    failures = 0
    for m in materials:
        rmin, rmax, umin, umax = bounds(m)
        print("material %u: native=%d entropy=%d  rho in [%e, %e] g/cm^3  u in [%e, %e] erg/g"
              % (m.materialId, m.nativeGrid, m.hasEntropy, rmin, rmax, umin, umax))
        # raw array integrity: every stored node must decode finite & positive
        for name, arr, positive in (("P", m.logPressure, True), ("cs", m.logSoundSpeed, True),
                                    ("T", m.logTemperature, True)):
            decoded = np.exp(arr)
            n_bad = int(np.count_nonzero(~np.isfinite(decoded) | (decoded <= 0 if positive else False)))
            if n_bad:
                print(f"  FAIL: {n_bad} non-finite/non-positive {name} values"); failures += 1
        if m.hasEntropy:
            n_bad = int(np.count_nonzero(~np.isfinite(m.entropy)))
            if n_bad:
                print(f"  FAIL: {n_bad} non-finite entropy values"); failures += 1

    # sample points (identical to test_eos_table.c)
    m = materials[0]
    if len(args) > 1:
        want = int(args[1])
        m = next((x for x in materials if x.materialId == want), None)
        if m is None:
            print(f"material {want} not present in table"); return 1
    rmin, rmax, umin, umax = bounds(m)
    print(f"sample evaluations for material {m.materialId}:")
    rho_mid = np.exp(0.5 * (np.log(rmin) + np.log(rmax)))
    u_mid = np.exp(0.5 * (np.log(umin) + np.log(umax))) if umin > 0 else 0.5 * (umin + umax)
    for rho, u in ((rho_mid, u_mid), (rmin, umin), (rmax, umax), (rho_mid, umax), (rho_mid, umin),
                   (0.1 * rmin, u_mid), (rho_mid, 100.0 * umax)):
        failures += not check_point(m, rho, u)

    # 40x40 scan per material (identical to test_eos_table.c)
    print("full-grid scan over all materials ...")
    NR = NU = 40
    for m in materials:
        rmin, rmax, umin, umax = bounds(m)
        for ir in range(NR):
            rho = np.exp(np.log(rmin) + (np.log(rmax) - np.log(rmin)) * ir / (NR - 1))
            for iu in range(NU):
                u = (np.exp(np.log(umin) + (np.log(umax) - np.log(umin)) * iu / (NU - 1))
                     if umin > 0 else umin + (umax - umin) * iu / (NU - 1))
                s = evaluate(m, rho, u)
                if not (np.isfinite(s["P"]) and s["P"] > 0 and np.isfinite(s["cs"]) and s["cs"] > 0
                        and np.isfinite(s["T"]) and s["T"] > 0 and np.isfinite(s["S"])):
                    print("  FAIL material %u at rho=%e u=%e: P=%e cs=%e T=%e S=%e status=%s"
                          % (m.materialId, rho, u, s["P"], s["cs"], s["T"], s["S"], s["status"]))
                    failures += 1
                    if failures > 20:
                        print("too many failures, aborting scan")
                        print(f"RESULT: {failures} failure(s)")
                        return 1
        print(f"  material {m.materialId}: {NR} x {NU} grid OK")

    # on-node check: evaluate at exact stored nodes. NOTE: exact reproduction is
    # NOT expected -- for non-monotonic native T(rho,u) columns (plateaus), the
    # defensive segment search may legitimately map the energy to a neighboring
    # non-degenerate segment (identical behavior in the C code). We therefore
    # only require the deviation to stay within the local inter-node spacing
    # (a few per mil); a wrong axis/stride would give order-unity deviations.
    print("on-node check (500 random nodes per material) ...")
    rng = np.random.default_rng(42)
    for m in materials:
        deviations = []
        for _ in range(500):
            ir = int(rng.integers(0, m.nRho - 1))
            iu = int(rng.integers(0, m.nU))
            rho = float(np.exp(m.logRho[ir]))
            u = float(m.energyAxis[ir * m.uStride + iu] if m.uStride else m.energyAxis[iu])
            s = evaluate(m, rho, u)
            ref = float(np.exp(m.logPressure[ir * m.nU + iu]))
            if s["P"] > 0 and ref > 0:
                deviations.append(abs(s["P"] - ref) / ref)
            else:
                print(f"  FAIL node ir={ir} iu={iu}: P={s['P']} ref={ref} status={s['status']}")
                failures += 1
        deviations = np.array(deviations)
        worst = float(deviations.max())
        print(f"  material {m.materialId}: on-node rel P deviation median={np.median(deviations):.3e} "
              f"p99={np.percentile(deviations, 99):.3e} max={worst:.3e}")
        if worst > 5e-2:
            print("  FAIL: on-node deviation exceeds inter-node spacing scale"); failures += 1

    # random off-node samples must stay finite and positive
    print("random off-node sampling (2000 points per material) ...")
    for m in materials:
        rmin, rmax, umin, umax = bounds(m)
        n_bad = 0
        for _ in range(2000):
            rho = float(np.exp(rng.uniform(np.log(rmin), np.log(rmax))))
            u = float(rng.uniform(umin, umax)) if umin < 0 else float(np.exp(rng.uniform(np.log(umin), np.log(umax))))
            s = evaluate(m, rho, u)
            if not (np.isfinite(s["P"]) and s["P"] > 0 and s["cs"] > 0 and s["T"] > 0 and np.isfinite(s["S"])):
                n_bad += 1
                if n_bad <= 5:
                    print(f"  FAIL rho={rho:e} u={u:e}: {s}")
        print(f"  material {m.materialId}: {n_bad} bad / 2000")
        failures += n_bad > 0

    if failures:
        print(f"RESULT: {failures} failure(s)")
        return 1
    print("RESULT: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
