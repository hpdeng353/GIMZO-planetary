# coding: utf-8
# Kegerreis et al. 2022 (ApJL 937:L40) impact angular momentum.
# Paper base scenario: M_t = 0.887 M_E, M_i = 0.133 M_E, beta = 45 deg,
# v_c = 1 v_esc, L_EM = 3.5e34 kg m^2/s; suite L = 1.19-1.37 L_EM.
# User asked for 0.887 + 0.113; we compute both impactor masses.
# EOS: paper uses ANEOS Fe85Si15 core + forsterite mantle, 30/70, T_s = 2000 K.
import woma
import numpy as np

R_earth = 6.371e6
M_earth = 5.9724e24
G = 6.67430e-11
L_EM = 3.5e34  # kg m^2/s, value used in Kegerreis+2022

def build(mass_frac, rmin, rmax, name):
    m = mass_frac * M_earth
    planet = woma.Planet(
        name=name,
        A1_mat_layer=["ANEOS_Fe85Si15", "ANEOS_forsterite"],
        A1_T_rho_type=["adiabatic", "adiabatic"],
        A1_M_layer=[0.3 * m, 0.7 * m],
        P_s=1e5,
        T_s=2000,
    )
    planet.gen_prof_L2_find_R_R1_given_M1_M2(R_min=rmin * R_earth,
                                             R_max=rmax * R_earth)
    return planet

target = build(0.887, 0.90, 1.05, "Target_0.887")
imp113 = build(0.113, 0.49, 0.55, "Impactor_0.113")
imp133 = build(0.133, 0.52, 0.60, "Impactor_0.133")

for p in (target, imp113, imp133):
    print(f"{p.name:16s}: M = {p.M/M_earth:.4f} M_E, R = {p.R/R_earth:.4f} R_E "
          f"= {p.R:.4e} m, R_core = {p.A1_r[p.A1_idx_layer[0]]/R_earth:.4f} R_E")

def impact_L(M_t, R_t, M_i, R_i, beta_deg, vfac=1.0):
    M_tot = M_t + M_i
    mu = M_t * M_i / M_tot
    R_sum = R_t + R_i
    v_esc = np.sqrt(2 * G * M_tot / R_sum)
    v_c = vfac * v_esc
    L = mu * v_c * R_sum * np.sin(np.radians(beta_deg))
    return v_esc, L

print()
for imp, tag in [(imp113, "0.887 + 0.113 (user)"), (imp133, "0.887 + 0.133 (paper)")]:
    v_esc, L45 = impact_L(target.M, target.R, imp.M, imp.R, 45.0)
    _, L_graze = impact_L(target.M, target.R, imp.M, imp.R, 90.0)
    print(f"--- {tag} ---")
    print(f"  v_esc            = {v_esc/1e3:.3f} km/s   (paper: ~9 km/s)")
    print(f"  L(beta=45, vesc) = {L45:.3e} kg m^2/s = {L45/L_EM:.3f} L_EM")
    print(f"  L(grazing, vesc) = {L_graze:.3e} kg m^2/s = {L_graze/L_EM:.3f} L_EM")
    for beta in (43, 44, 46, 47, 48):
        _, L = impact_L(target.M, target.R, imp.M, imp.R, beta)
        print(f"  L(beta={beta}, vesc) = {L/L_EM:.3f} L_EM")
    print()
