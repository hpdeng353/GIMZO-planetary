#!/usr/bin/env python3
"""Build one or two two-resolution WoMa rocky planets as GIZMO-format HDF5 ICs for planetg.

Ported from sphexa's makeplanet_sphexa.py. The WoMa construction (adiabatic
two-layer profiles, CMB particle-mass matching, complete-realization search,
thinning) is unchanged; only the unit system and the output writer differ.

Output format: GADGET-3/GIZMO HDF5 initial conditions with a Header group and
PartType0 datasets Coordinates, Velocities, Masses, InternalEnergy, ParticleIDs,
Materials, Temperature and SmoothingLength. planetg (noon1.conf) reads the first
five plus Materials (READ_IMAT) and Temperature (EOS_CARRIES_TEMPERATURE);
SmoothingLength is informational (enable READ_HSML to consume it).

Materials convention: the Materials dataset carries small integer imat indices
(default: mantle = 0, core = 1), matching the EosTableMatIds parameter order in
the planetg parameter file. With the defaults, set
    EosTableMatIds    62,63
so imat 0 -> material 62 (ANEOS forsterite) and imat 1 -> material 63 (ANEOS iron).

Units: planetg code units are read from the parameter file values
(UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s); the defaults below
match noon1.params. All particle arrays are written in these code units. planetg
takes GRAVITY = 6.672e-8 cgs, so G = 1.00004 in the default units -- WoMa builds
the profile with its own SI G; the resulting ~0.03% gravity mismatch is far below
SPH discretization noise and is absorbed by relaxation.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import sys

import numpy as np


# WoMa works in SI; only the output unit system differs from the sphexa script.
M_EARTH_KG = 5.9724e24
EARTH_RADIUS_M = 6.37869e6
ROCKY_MASS_RADIUS_EXPONENT = 0.27

# Defaults match noon1.params.
DEFAULT_UNIT_LENGTH_CM = 6.37869e8
DEFAULT_UNIT_MASS_G = 9.56072e25
DEFAULT_UNIT_VELOCITY_CMS = 1.0e5


def planetg_units(length_cm: float, mass_g: float, velocity_cms: float) -> dict[str, float]:
    """Return the planetg code-unit system with SI conversion factors."""
    if not (math.isfinite(length_cm) and length_cm > 0.0):
        raise ValueError("length unit must be finite and positive")
    if not (math.isfinite(mass_g) and mass_g > 0.0):
        raise ValueError("mass unit must be finite and positive")
    if not (math.isfinite(velocity_cms) and velocity_cms > 0.0):
        raise ValueError("velocity unit must be finite and positive")

    length_m = length_cm * 1.0e-2
    mass_kg = mass_g * 1.0e-3
    velocity_ms = velocity_cms * 1.0e-2
    energy_jkg = velocity_ms**2
    return {
        "length_m": length_m,
        "velocity_ms": velocity_ms,
        "time_s": length_m / velocity_ms,
        "mass_kg": mass_kg,
        "density_kgm3": mass_kg / length_m**3,
        "specific_energy_jkg": energy_jkg,
        "pressure_pa": mass_kg / length_m**3 * energy_jkg,
        "earth_mass_code": M_EARTH_KG / mass_kg,
        "unit_length_cm": length_cm,
        "unit_mass_g": mass_g,
        "unit_velocity_cms": velocity_cms,
    }


def cmb_density_ratio(planet) -> tuple[float, int]:
    """Read the one-sided density contrast from a two-layer WoMa profile."""
    boundary = int(np.asarray(planet.A1_idx_layer).reshape(-1)[0])
    density = np.asarray(planet.A1_rho, dtype=np.float64)
    if boundary < 0 or boundary + 1 >= density.size:
        raise ValueError("WoMa returned an invalid core--mantle boundary index")
    ratio = float(density[boundary] / density[boundary + 1])
    if not math.isfinite(ratio) or ratio <= 0.0:
        raise ValueError("WoMa returned an invalid CMB density contrast")
    return ratio, boundary


def particle_set(woma, planet, count: int, neighbors: int):
    """Support both the older ParticlePlanet and current ParticleSet API."""
    if count < 1:
        raise ValueError("particle count must be positive")
    constructor = getattr(woma, "ParticlePlanet", None)
    if constructor is not None:
        try:
            return constructor(planet, int(count), N_ngb=neighbors, verbosity=0)
        except TypeError:
            return constructor(planet, int(count), verbosity=0)

    constructor = getattr(woma, "ParticleSet", None)
    if constructor is None:
        raise RuntimeError("WoMa provides neither ParticlePlanet nor ParticleSet")
    try:
        return constructor(planet, int(count), N_ngb=neighbors)
    except TypeError:
        return constructor(planet, int(count))


def layer_masses(mass_earth: float, core_mass_fraction: float) -> tuple[float, float, float]:
    """Return total, iron-core, and rock-mantle masses in kg."""
    total_mass_kg = mass_earth * M_EARTH_KG
    core_mass_kg = core_mass_fraction * total_mass_kg
    return total_mass_kg, core_mass_kg, total_mass_kg - core_mass_kg


def make_planet_profile(woma, name: str, mass_earth: float, args):
    """Construct one two-layer WoMa radial profile."""
    total_mass_kg, core_mass_kg, mantle_mass_kg = layer_masses(mass_earth, args.core_mass_fraction)
    radius_guess_earth = mass_earth**ROCKY_MASS_RADIUS_EXPONENT
    is_secondary = name.lower() == "secondary"
    explicit_min = args.secondary_radius_min_earth if is_secondary else args.radius_min_earth
    explicit_max = args.secondary_radius_max_earth if is_secondary else args.radius_max_earth
    if explicit_min is None:
        radius_min_earth = 0.95 * radius_guess_earth
        radius_max_earth = 1.05 * radius_guess_earth
        radius_source = "mass-scaled automatic"
    else:
        radius_min_earth = explicit_min
        radius_max_earth = explicit_max
        radius_source = "explicit"

    # WoMa requires usable trial planets at both bounds before bisection. Refine
    # only the rejected endpoint, retaining the rest of the original bracket.
    max_attempts = 1 if radius_source == "explicit" else 8
    last_error = None
    succeeded = False
    rejected_min = None
    rejected_max = None
    usable_min = None
    usable_max = None
    for attempt in range(max_attempts):
        print(
            f"{name} WoMa radius bracket ({radius_source}, attempt {attempt + 1}/{max_attempts}): "
            f"{radius_min_earth:.6g}--{radius_max_earth:.6g} Earth radii"
        )
        planet = woma.Planet(
            name=name,
            A1_mat_layer=["ANEOS_iron", "ANEOS_forsterite"],
            A1_T_rho_type=["adiabatic", "adiabatic"],
            A1_M_layer=[core_mass_kg, mantle_mass_kg],
            P_s=args.surface_pressure_pa,
            T_s=args.surface_temperature_k,
        )
        try:
            planet.gen_prof_L2_find_R_R1_given_M1_M2(
                R_min=radius_min_earth * EARTH_RADIUS_M,
                R_max=radius_max_earth * EARTH_RADIUS_M,
            )
            succeeded = True
            break
        except (RuntimeError, ValueError) as error:
            last_error = error
            message = str(error)
            if radius_source == "explicit" or attempt + 1 == max_attempts:
                break
            if "Could not build" in message and "R=R_max" in message:
                rejected_max = radius_max_earth
                lower = usable_max if usable_max is not None else radius_min_earth
                radius_max_earth = 0.5 * (lower + rejected_max)
                print(f"{name}: WoMa rejected R_max; moving only R_max halfway inward.")
            elif "tends to R_max" in message:
                usable_max = radius_max_earth
                if rejected_max is None:
                    radius_max_earth *= 1.05
                else:
                    radius_max_earth = 0.5 * (usable_max + rejected_max)
                print(f"{name}: solution tends to R_max; increasing only R_max.")
            elif "Could not build" in message and "R=R_min" in message:
                rejected_min = radius_min_earth
                upper = usable_min if usable_min is not None else radius_max_earth
                radius_min_earth = 0.5 * (rejected_min + upper)
                print(f"{name}: WoMa rejected R_min; moving only R_min halfway inward.")
            elif "tends to R_min" in message:
                usable_min = radius_min_earth
                if rejected_min is None:
                    radius_min_earth *= 0.95
                else:
                    radius_min_earth = 0.5 * (rejected_min + usable_min)
                print(f"{name}: solution tends to R_min; decreasing only R_min.")
            else:
                break

    if not succeeded:
        role = "secondary" if is_secondary else "primary"
        raise RuntimeError(
            f"WoMa could not construct the {role} ({mass_earth:.6g} Earth masses); "
            "last bracket was "
            f"{radius_min_earth:.6g}--{radius_max_earth:.6g} Earth radii. Override both "
            f"--{role}-radius-min-earth and --{role}-radius-max-earth. "
            f"Original WoMa error: {last_error}"
        ) from last_error
    ratio, _ = cmb_density_ratio(planet)
    return planet, total_mass_kg, core_mass_kg, mantle_mass_kg, ratio


def primary_resolution(
    total_mass_kg: float,
    core_mass_kg: float,
    mantle_mass_kg: float,
    cmb_ratio: float,
    approximate_particle_count: int,
) -> tuple[int, int, float, float]:
    """Choose the two full-realization sizes for an approximate merged count."""
    if approximate_particle_count < 2:
        raise ValueError("primary-particle-count must be at least two")

    mantle_particle_mass_kg = (core_mass_kg / cmb_ratio + mantle_mass_kg) / approximate_particle_count
    core_particle_mass_kg = cmb_ratio * mantle_particle_mass_kg
    high_count = max(1, round(total_mass_kg / mantle_particle_mass_kg))
    low_count = max(1, round(total_mass_kg / core_particle_mass_kg))
    return low_count, high_count, core_particle_mass_kg, mantle_particle_mass_kg


def secondary_layer_counts(
    core_mass_kg: float,
    mantle_mass_kg: float,
    core_particle_mass_kg: float,
    mantle_particle_mass_kg: float,
) -> tuple[int, int]:
    """Round requested secondary layer masses to primary particle masses."""
    return (
        max(1, round(core_mass_kg / core_particle_mass_kg)),
        max(1, round(mantle_mass_kg / mantle_particle_mass_kg)),
    )


def resolution_trial_counts(target_layer_count: int, layer_mass_fraction: float, trials: int) -> list[int]:
    """Return center-out full-planet resolutions likely to bracket the requested layer count."""
    if target_layer_count < 1 or not 0.0 < layer_mass_fraction < 1.0:
        raise ValueError("invalid target layer count or mass fraction")
    if trials < 1:
        raise ValueError("secondary-resolution-trials must be positive")

    center = max(1, round(target_layer_count / layer_mass_fraction))
    # Four binomial standard deviations cover normal realization-to-realization count fluctuations. Keep a useful
    # minimum range for deterministic particle generators and cap the range for very poorly resolved layers.
    fractional_span = min(
        0.25,
        max(0.02, 4.0 * math.sqrt((1.0 - layer_mass_fraction) / target_layer_count)),
    )
    levels = max(1, math.ceil((trials - 1) / 2))
    candidates = [center]
    for level in range(1, levels + 1):
        fractional_offset = fractional_span * level / levels
        candidates.extend(
            (
                max(1, round(center * (1.0 - fractional_offset))),
                max(1, round(center * (1.0 + fractional_offset))),
            )
        )

    unique = []
    for count in candidates:
        if count not in unique:
            unique.append(count)
        if len(unique) == trials:
            break
    next_count = center
    while len(unique) < trials:
        next_count += 1
        if next_count not in unique:
            unique.append(next_count)
    return unique


def best_complete_material_realization(
    woma,
    planet,
    target_layer_count: int,
    layer_mass_fraction: float,
    neighbors: int,
    woma_material_id: int,
    trials: int,
    material_name: str,
    thin_to_target: bool = False,
    overshoot_max: float = 0.0,
    rng: np.random.Generator | None = None,
) -> tuple[dict[str, np.ndarray], int]:
    """Search nearby resolutions and retain one complete material subset.

    With thin_to_target enabled, only oversized realizations are accepted
    (acceptance ladder 1%, 2%, ... of the target per trial, capped at
    overshoot_max), and the smallest-overshoot candidate is randomly thinned to
    exactly target_layer_count particles so the secondary mass is exact.
    """
    pending = resolution_trial_counts(target_layer_count, layer_mass_fraction, trials)
    attempted: set[int] = set()
    best_material = None
    best_source_count = 0
    best_score = None

    for trial in range(1, trials + 1):
        while pending and pending[0] in attempted:
            pending.pop(0)
        if not pending:
            candidate = max(attempted) + 1
            while candidate in attempted:
                candidate += 1
        else:
            candidate = pending.pop(0)

        attempted.add(candidate)
        realization = particle_set(woma, planet, candidate, neighbors)
        material = extract_material(realization, planet, woma_material_id)
        del realization
        actual_count = len(material["u_jkg"])
        count_error = actual_count - target_layer_count
        print(
            f"  {material_name} trial {trial}/{trials}: full realization {candidate:,}, "
            f"retained whole layer {actual_count:,}, target {target_layer_count:,}, "
            f"difference {count_error:+,}"
        )

        if thin_to_target:
            if count_error < 0:
                # Undersized realizations cannot be thinned to the target.
                del material
            else:
                overshoot = count_error / target_layer_count
                score = (overshoot, abs(candidate - round(target_layer_count / layer_mass_fraction)), candidate)
                if best_score is None or score < best_score:
                    if best_material is not None:
                        del best_material
                    best_material = material
                    best_source_count = candidate
                    best_score = score
                else:
                    del material
                if count_error == 0 or overshoot <= min(0.01 * trial, overshoot_max):
                    break
        else:
            score = (abs(count_error), abs(candidate - round(target_layer_count / layer_mass_fraction)), candidate)
            if best_score is None or score < best_score:
                best_material = material
                best_source_count = candidate
                best_score = score
            else:
                del material

            if count_error == 0:
                break

        # The observed layer fraction provides a better resolution estimate for the next attempt. The precomputed
        # center-out candidates remain as fallbacks if this correction repeats a resolution already tried.
        corrected = max(1, round(candidate * target_layer_count / actual_count))
        if corrected not in attempted and corrected not in pending:
            pending.insert(0, corrected)

    if best_material is None:
        if thin_to_target:
            raise RuntimeError(
                f"no oversized secondary {material_name} realization was generated in {trials} trials; "
                f"increase --secondary-resolution-trials"
            )
        raise RuntimeError(f"no usable complete secondary {material_name} realization was generated")

    if thin_to_target:
        actual_count = len(best_material["u_jkg"])
        overshoot = (actual_count - target_layer_count) / target_layer_count
        if overshoot > overshoot_max:
            raise RuntimeError(
                f"smallest oversized secondary {material_name} realization overshoots the target by "
                f"{overshoot:.2%}, above --secondary-thin-overshoot-max {overshoot_max:.2%}; "
                f"increase --secondary-resolution-trials or relax the cap"
            )
        if actual_count > target_layer_count:
            if rng is None:
                rng = np.random.default_rng(0)
            keep = rng.choice(actual_count, target_layer_count, replace=False)
            keep.sort()
            deleted = actual_count - target_layer_count
            for key in ("position_m", "velocity_ms", "u_jkg", "density_kgm3", "temperature_k"):
                if key in best_material:
                    best_material[key] = best_material[key][keep]
            print(
                f"  thinned {material_name}: {actual_count:,} -> {target_layer_count:,} "
                f"(deleted {deleted:,} = {deleted / actual_count:.2%})"
            )

    print(
        f"  selected complete {material_name}: full realization {best_source_count:,}, "
        f"layer particles {len(best_material['u_jkg']):,}"
    )
    return best_material, best_source_count


def extract_material(realization, planet, woma_material_id: int) -> dict[str, np.ndarray]:
    """Center one realization, then extract a material without shifting layers."""
    material = np.asarray(realization.A1_mat_id)
    mask = material == woma_material_id
    if not np.any(mask):
        available = np.unique(material).tolist()
        raise ValueError(
            f"WoMa material ID {woma_material_id} is absent; available IDs are {available}"
        )

    position = np.asarray(realization.A2_pos, dtype=np.float64)
    velocity = np.asarray(realization.A2_vel, dtype=np.float64)
    mass = np.asarray(realization.A1_m, dtype=np.float64)
    if position.ndim != 2 or position.shape[1] != 3 or velocity.shape != position.shape:
        raise ValueError("WoMa particle positions/velocities do not have shape (N, 3)")
    if mass.shape != (position.shape[0],):
        raise ValueError("WoMa particle masses do not have shape (N,)")

    # Center the full realization, not each selected layer independently.
    center = np.average(position, axis=0, weights=mass)
    bulk_velocity = np.average(velocity, axis=0, weights=mass)
    position = position - center
    velocity = velocity - bulk_velocity

    if hasattr(realization, "A1_rho") and np.asarray(realization.A1_rho).shape == mass.shape:
        density = np.asarray(realization.A1_rho, dtype=np.float64)
    else:
        profile_r = np.asarray(planet.A1_r, dtype=np.float64)
        profile_density = np.asarray(planet.A1_rho, dtype=np.float64)
        if profile_r.ndim != 1 or profile_density.shape != profile_r.shape or profile_r.size < 2:
            raise ValueError("WoMa provides neither particle densities nor a valid radial density profile")
        density = np.interp(np.linalg.norm(position, axis=1), profile_r, profile_density)
    if density.shape != mass.shape or np.any(~np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("WoMa returned invalid particle densities")

    result = {
        "position_m": position[mask].copy(),
        "velocity_ms": velocity[mask].copy(),
        "u_jkg": np.asarray(realization.A1_u, dtype=np.float64)[mask].copy(),
        "density_kgm3": density[mask].copy(),
    }

    # Per-particle temperature is written to the ICs when WoMa provides it
    # (planetg reads Temperature at RestartFlag 0 under EOS_CARRIES_TEMPERATURE).
    if hasattr(realization, "A1_T") and np.asarray(realization.A1_T).shape == mass.shape:
        temperature = np.asarray(realization.A1_T, dtype=np.float64)
        if np.all(np.isfinite(temperature)) and np.all(temperature > 0.0):
            result["temperature_k"] = temperature[mask].copy()
    return result


def smoothing_lengths(mass_kg: np.ndarray, density_kgm3: np.ndarray, neighbors: int) -> np.ndarray:
    """Initial h for a kernel whose compact support radius is 2h."""
    eta = (3.0 * neighbors / (32.0 * math.pi)) ** (1.0 / 3.0)
    return eta * np.cbrt(mass_kg / density_kgm3)


def assemble_particles(
    core: dict[str, np.ndarray],
    mantle: dict[str, np.ndarray],
    core_mass_kg: float,
    mantle_mass_kg: float,
    core_eos_id: int,
    mantle_eos_id: int,
    core_imat: int,
    mantle_imat: int,
    neighbors: int,
    units: dict[str, float],
    fixed_core_particle_mass_kg: float | None = None,
    fixed_mantle_particle_mass_kg: float | None = None,
) -> dict[str, np.ndarray]:
    """Assign layer particle masses and convert WoMa SI arrays to code units."""
    n_core = len(core["u_jkg"])
    n_mantle = len(mantle["u_jkg"])
    if n_core == 0 or n_mantle == 0:
        raise ValueError("both core and mantle must contain particles")

    if (fixed_core_particle_mass_kg is None) != (fixed_mantle_particle_mass_kg is None):
        raise ValueError("fixed core and mantle particle masses must be supplied together")
    if fixed_core_particle_mass_kg is None:
        core_particle_mass = core_mass_kg / n_core
        mantle_particle_mass = mantle_mass_kg / n_mantle
    else:
        core_particle_mass = fixed_core_particle_mass_kg
        mantle_particle_mass = fixed_mantle_particle_mass_kg
    if core_particle_mass <= 0.0 or mantle_particle_mass <= 0.0:
        raise ValueError("core and mantle particle masses must be positive")
    physical_mass = np.concatenate(
        (
            np.full(n_core, core_particle_mass, dtype=np.float64),
            np.full(n_mantle, mantle_particle_mass, dtype=np.float64),
        )
    )
    physical_density = np.concatenate((core["density_kgm3"], mantle["density_kgm3"]))
    h_m = smoothing_lengths(physical_mass, physical_density, neighbors)

    position = np.concatenate((core["position_m"], mantle["position_m"]))
    velocity = np.concatenate((core["velocity_ms"], mantle["velocity_ms"]))
    u_jkg = np.concatenate((core["u_jkg"], mantle["u_jkg"]))
    material = np.concatenate(
        (
            np.full(n_core, core_eos_id, dtype=np.uint32),
            np.full(n_mantle, mantle_eos_id, dtype=np.uint32),
        )
    )
    imat = np.concatenate(
        (
            np.full(n_core, core_imat, dtype=np.uint32),
            np.full(n_mantle, mantle_imat, dtype=np.uint32),
        )
    )

    particles = {
        "x": position[:, 0] / units["length_m"],
        "y": position[:, 1] / units["length_m"],
        "z": position[:, 2] / units["length_m"],
        "vx": velocity[:, 0] / units["velocity_ms"],
        "vy": velocity[:, 1] / units["velocity_ms"],
        "vz": velocity[:, 2] / units["velocity_ms"],
        "u": u_jkg / units["specific_energy_jkg"],
        "m": physical_mass / units["mass_kg"],
        "h": h_m / units["length_m"],
        "materialId": material,
        "imat": imat,
        "id": np.arange(n_core + n_mantle, dtype=np.uint64),
        "density_si": physical_density,
        "n_core": n_core,
    }
    if "temperature_k" in core and "temperature_k" in mantle:
        particles["temperature_k"] = np.concatenate((core["temperature_k"], mantle["temperature_k"]))
    return particles


def validate_particles(
    particles: dict[str, np.ndarray], args, ratio: float, enforce_mass_ratio: bool = True
) -> dict[str, float]:
    """Reject malformed or poorly matched merged realizations before writing."""
    n = len(particles["id"])
    required = ("x", "y", "z", "vx", "vy", "vz", "u", "m", "h", "materialId", "imat")
    for name in required:
        values = particles[name]
        if len(values) != n:
            raise ValueError(f"field {name} has {len(values)} entries instead of {n}")
        if name not in ("materialId", "imat") and np.any(~np.isfinite(values)):
            raise ValueError(f"field {name} contains non-finite values")
    if np.any(particles["m"] <= 0.0) or np.any(particles["h"] <= 0.0):
        raise ValueError("particle mass and smoothing length must be positive")
    if np.unique(particles["id"]).size != n:
        raise ValueError("particle IDs are not unique")

    n_core = int(particles["n_core"])
    core_particle_mass = float(particles["m"][0])
    mantle_particle_mass = float(particles["m"][n_core])
    actual_ratio = core_particle_mass / mantle_particle_mass
    relative_error = abs(actual_ratio / ratio - 1.0)
    if enforce_mass_ratio and relative_error > args.mass_ratio_tolerance:
        raise ValueError(
            f"particle-mass ratio {actual_ratio:.6g} differs from CMB density ratio "
            f"{ratio:.6g} by {100.0 * relative_error:.2f}%; change particle counts or tolerance"
        )

    radius = np.sqrt(particles["x"] ** 2 + particles["y"] ** 2 + particles["z"] ** 2)
    core_outer = float(np.max(radius[:n_core]))
    mantle_inner = float(np.min(radius[n_core:]))
    boundary_count = max(16, min(512, n_core, n - n_core))
    core_h = float(np.median(particles["h"][:n_core][np.argsort(radius[:n_core])[-boundary_count:]]))
    mantle_h = float(
        np.median(particles["h"][n_core:][np.argsort(radius[n_core:])[:boundary_count]])
    )
    return {
        "core_particle_mass_code": core_particle_mass,
        "mantle_particle_mass_code": mantle_particle_mass,
        "particle_mass_ratio": actual_ratio,
        "mass_ratio_relative_error": relative_error,
        "core_outer_radius": core_outer,
        "mantle_inner_radius": mantle_inner,
        "cmb_gap_over_h": (mantle_inner - core_outer) / (0.5 * (core_h + mantle_h)),
        "cmb_h_ratio": core_h / mantle_h,
    }


def realized_masses(particles: dict[str, np.ndarray], units: dict[str, float]) -> tuple[float, float, float]:
    """Return realized total, core, and mantle masses in kg."""
    n_core = int(particles["n_core"])
    mass_kg = np.asarray(particles["m"], dtype=np.float64) * units["mass_kg"]
    core_mass_kg = float(np.sum(mass_kg[:n_core], dtype=np.float64))
    mantle_mass_kg = float(np.sum(mass_kg[n_core:], dtype=np.float64))
    return core_mass_kg + mantle_mass_kg, core_mass_kg, mantle_mass_kg


def model_metadata(
    role: str,
    requested_mass_earth: float,
    core_mass_fraction: float,
    particles: dict[str, np.ndarray],
    units: dict[str, float],
    particle_mass_reference: str,
) -> dict[str, float | str]:
    """Create traceable paired-model metadata for the HDF5 step."""
    n_core = int(particles["n_core"])
    realized_total_kg, _, _ = realized_masses(particles, units)
    return {
        "planetModelRole": role,
        "requestedMassEarth": requested_mass_earth,
        "realizedMassEarth": realized_total_kg / M_EARTH_KG,
        "requestedCoreMassFraction": core_mass_fraction,
        "ironParticleMassKg": float(particles["m"][0]) * units["mass_kg"],
        "rockParticleMassKg": float(particles["m"][n_core]) * units["mass_kg"],
        "particleMassReference": particle_mass_reference,
    }


def write_gizmo(
    path: Path,
    particles: dict[str, np.ndarray],
    args,
    units: dict[str, float],
    metadata: dict[str, float | str] | None = None,
) -> None:
    """Atomically write a GIZMO (GADGET-3 HDF5) initial-conditions file."""
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("writing requires h5py (install it in the WoMa Python environment)") from error

    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not args.force:
        raise FileExistsError(f"output already exists: {path}; pass --force to replace it")
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"temporary output already exists: {temporary}")

    n = len(particles["id"])
    coordinates = np.column_stack((particles["x"], particles["y"], particles["z"]))
    velocities = np.column_stack((particles["vx"], particles["vy"], particles["vz"]))

    try:
        with h5py.File(temporary, "w") as output:
            header = output.create_group("Header")
            header.attrs["NumPart_ThisFile"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.int32)
            header.attrs["NumPart_Total"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.uint32)
            header.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
            header.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
            header.attrs["Time"] = 0.0
            header.attrs["Redshift"] = 0.0
            header.attrs["BoxSize"] = float(args.box_size)
            header.attrs["NumFilesPerSnapshot"] = 1
            header.attrs["Omega0"] = 0.0
            header.attrs["OmegaLambda"] = 0.0
            header.attrs["HubbleParam"] = 1.0
            header.attrs["Flag_Sfr"] = 0
            header.attrs["Flag_Cooling"] = 0
            header.attrs["Flag_StellarAge"] = 0
            header.attrs["Flag_Metals"] = 0
            header.attrs["Flag_Feedback"] = 0
            header.attrs["Flag_DoublePrecision"] = 0
            header.attrs["Flag_IC_Info"] = 0

            # Unit system and provenance (informational; planetg reads its units
            # from the parameter file, which must match these values).
            header.attrs["UnitLength_in_cm"] = units["unit_length_cm"]
            header.attrs["UnitMass_in_g"] = units["unit_mass_g"]
            header.attrs["UnitVelocity_in_cm_per_s"] = units["unit_velocity_cms"]
            header.attrs["unitTimeS"] = units["time_s"]
            if metadata:
                for name, value in metadata.items():
                    header.attrs[name] = value

            group = output.create_group("PartType0")
            group.create_dataset("Coordinates", data=coordinates, dtype=np.float64)
            group.create_dataset("Velocities", data=velocities, dtype=np.float64)
            group.create_dataset("Masses", data=np.asarray(particles["m"]), dtype=np.float64)
            group.create_dataset("InternalEnergy", data=np.asarray(particles["u"]), dtype=np.float64)
            group.create_dataset("ParticleIDs", data=particles["id"], dtype=np.uint64)
            # Materials carries the imat index used by EosTableMatIds, not the raw
            # EOS material ID (the raw IDs are recorded in EosMaterialId).
            group.create_dataset("Materials", data=particles["imat"], dtype=np.uint32)
            group.create_dataset("EosMaterialId", data=particles["materialId"], dtype=np.uint32)
            group.create_dataset(
                "SmoothingLength", data=np.asarray(particles["h"]), dtype=np.float64
            )
            if "temperature_k" in particles:
                group.create_dataset(
                    "Temperature", data=np.asarray(particles["temperature_k"]), dtype=np.float64
                )

        with h5py.File(temporary, "r") as check:
            expected = {
                "Coordinates", "Velocities", "Masses", "InternalEnergy",
                "ParticleIDs", "Materials", "SmoothingLength",
            }
            group = check["PartType0"]
            missing = expected.difference(group.keys())
            wrong_size = [name for name in expected if name in group and group[name].shape[0] != n]
            header = check["Header"]
            header_ok = (
                int(header.attrs["NumPart_ThisFile"][0]) == n
                and int(header.attrs["NumPart_Total"][0]) == n
                and float(header.attrs["MassTable"][0]) == 0.0
            )
            if missing or wrong_size or not header_ok:
                raise RuntimeError(
                    "written GIZMO file failed validation; "
                    f"missing datasets={sorted(missing)}, wrong-size={sorted(wrong_size)}, "
                    f"header consistent={header_ok}"
                )
        os.replace(temporary, path)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def parse_arguments(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--primary-output", "--output", dest="output", type=Path, default=Path("earth_085_planetg.h5")
    )
    parser.add_argument("--secondary-output", type=Path, default=Path("secondary_planetg.h5"))
    parser.add_argument("--force", action="store_true", help="replace existing output file(s)")
    parser.add_argument(
        "--primary-mass-earth", "--planet-mass-earth", dest="planet_mass_earth", type=float, default=0.85
    )
    parser.add_argument("--secondary-mass-earth", type=float, help="enable paired mode with this secondary mass")
    parser.add_argument(
        "--primary-particle-count",
        type=int,
        help="approximate merged particle count for the primary; required in paired mode",
    )
    parser.add_argument("--core-mass-fraction", type=float, default=0.30)
    parser.add_argument("--surface-pressure-pa", type=float, default=1.0e5)
    parser.add_argument("--surface-temperature-k", type=float, default=2000.0)
    parser.add_argument(
        "--primary-radius-min-earth", "--radius-min-earth", dest="radius_min_earth",
        type=float,
        help="explicit minimum WoMa search radius in Earth radii (default: 0.95 M^0.27)",
    )
    parser.add_argument(
        "--primary-radius-max-earth", "--radius-max-earth", dest="radius_max_earth",
        type=float,
        help="explicit maximum WoMa search radius in Earth radii (default: 1.05 M^0.27)",
    )
    parser.add_argument(
        "--secondary-radius-min-earth",
        type=float,
        help="explicit secondary minimum WoMa search radius in Earth radii",
    )
    parser.add_argument(
        "--secondary-radius-max-earth",
        type=float,
        help="explicit secondary maximum WoMa search radius in Earth radii",
    )
    parser.add_argument("--mantle-particle-mass-code", type=float, default=2.65089e-5)
    parser.add_argument("--high-resolution-count", type=int)
    parser.add_argument("--low-resolution-count", type=int)
    parser.add_argument("--mass-ratio-tolerance", type=float, default=0.05)
    parser.add_argument(
        "--secondary-mass-tolerance",
        type=float,
        default=0.01,
        help="maximum relative layer or total mass error after complete-realization search (default: 0.01)",
    )
    parser.add_argument(
        "--secondary-resolution-trials",
        type=int,
        default=7,
        help="maximum complete WoMa realizations tested per secondary layer (default: 7)",
    )
    parser.add_argument(
        "--secondary-thin-overshoot-max",
        type=float,
        default=0.05,
        help="accept oversized secondary layer realizations (acceptance ladder 1%%, 2%%, ... per trial "
        "up to this cap) and randomly thin them to the exact target count (default: 0.05; 0 restores "
        "the legacy unmodified-realization selection)",
    )
    parser.add_argument("--woma-core-id", type=int, default=401)
    parser.add_argument("--woma-mantle-id", type=int, default=400)
    parser.add_argument("--eos-core-id", type=int, default=63)
    parser.add_argument("--eos-mantle-id", type=int, default=62)
    parser.add_argument(
        "--imat-core",
        type=int,
        default=1,
        help="Materials value written for iron-core particles (default: 1)",
    )
    parser.add_argument(
        "--imat-mantle",
        type=int,
        default=0,
        help="Materials value written for rock-mantle particles (default: 0)",
    )
    parser.add_argument(
        "--neighbors",
        type=int,
        default=100,
        help="neighbor count assumed for the informational SmoothingLength dataset (default: 100)",
    )
    parser.add_argument(
        "--box-size",
        type=float,
        default=8.0,
        help="BoxSize header attribute in code units (informational for non-periodic runs; default: 8)",
    )
    parser.add_argument("--unit-length-cm", type=float, default=DEFAULT_UNIT_LENGTH_CM)
    parser.add_argument("--unit-mass-g", type=float, default=DEFAULT_UNIT_MASS_G)
    parser.add_argument("--unit-velocity-cms", type=float, default=DEFAULT_UNIT_VELOCITY_CMS)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--dry-run", action="store_true", help="build only the radial profile and report counts")
    parser.add_argument("--show-units", action="store_true", help="report code units without importing WoMa")
    return parser.parse_args(argv)


def report_units(units: dict[str, float]) -> None:
    print("planetg code units:")
    for key in (
        "unit_length_cm", "unit_mass_g", "unit_velocity_cms",
        "length_m", "velocity_ms", "time_s", "mass_kg", "earth_mass_code",
        "density_kgm3", "specific_energy_jkg", "pressure_pa",
    ):
        print(f"  {key:28s} = {units[key]:.12g}")


def main(argv: list[str] | None = None) -> int:
    args = parse_arguments(argv)
    units = planetg_units(args.unit_length_cm, args.unit_mass_g, args.unit_velocity_cms)
    if args.show_units:
        report_units(units)
        return 0

    paired = args.secondary_mass_earth is not None
    if not (0.0 < args.core_mass_fraction < 1.0):
        raise ValueError("core mass fraction must lie between zero and one")
    if not math.isfinite(args.planet_mass_earth) or args.planet_mass_earth <= 0.0:
        raise ValueError("primary mass must be finite and positive")
    if args.primary_particle_count is not None and args.primary_particle_count < 2:
        raise ValueError("primary-particle-count must be at least two")
    if args.primary_particle_count is not None and (
        args.high_resolution_count is not None or args.low_resolution_count is not None
    ):
        raise ValueError("primary-particle-count cannot be combined with explicit high/low realization counts")
    if paired:
        if not math.isfinite(args.secondary_mass_earth) or args.secondary_mass_earth <= 0.0:
            raise ValueError("secondary mass must be finite and positive")
        if args.secondary_mass_earth > args.planet_mass_earth:
            raise ValueError("primary-mass-earth must be at least as large as secondary-mass-earth")
        if args.primary_particle_count is None:
            raise ValueError("paired mode requires --primary-particle-count")
        if not math.isfinite(args.secondary_mass_tolerance) or args.secondary_mass_tolerance <= 0.0:
            raise ValueError("secondary-mass-tolerance must be finite and positive")
        if args.secondary_resolution_trials < 1:
            raise ValueError("secondary-resolution-trials must be positive")
        if (
            not math.isfinite(args.secondary_thin_overshoot_max)
            or not 0.0 <= args.secondary_thin_overshoot_max <= 0.5
        ):
            raise ValueError("secondary-thin-overshoot-max must be finite and lie between 0 and 0.5")
        if args.output.resolve() == args.secondary_output.resolve():
            raise ValueError("primary and secondary output paths must be different")
    elif not math.isfinite(args.mantle_particle_mass_code) or args.mantle_particle_mass_code <= 0.0:
        raise ValueError("mantle particle mass must be finite and positive")
    if not math.isfinite(args.mass_ratio_tolerance) or args.mass_ratio_tolerance <= 0.0:
        raise ValueError("mass-ratio-tolerance must be finite and positive")
    if (args.radius_min_earth is None) != (args.radius_max_earth is None):
        raise ValueError("primary radius minimum and maximum must be supplied together")
    if args.radius_min_earth is not None and not (0.0 < args.radius_min_earth < args.radius_max_earth):
        raise ValueError("require 0 < primary-radius-min-earth < primary-radius-max-earth")
    if (args.secondary_radius_min_earth is None) != (args.secondary_radius_max_earth is None):
        raise ValueError("secondary radius minimum and maximum must be supplied together")
    if args.secondary_radius_min_earth is not None and not (
        0.0 < args.secondary_radius_min_earth < args.secondary_radius_max_earth
    ):
        raise ValueError("require 0 < secondary-radius-min-earth < secondary-radius-max-earth")
    if args.high_resolution_count is not None and args.high_resolution_count < 1:
        raise ValueError("high-resolution-count must be positive")
    if args.low_resolution_count is not None and args.low_resolution_count < 1:
        raise ValueError("low-resolution-count must be positive")
    if args.neighbors < 1:
        raise ValueError("neighbors must be positive")
    if not (0 <= args.imat_core <= 255 and 0 <= args.imat_mantle <= 255):
        raise ValueError("imat values must lie between 0 and 255")
    if args.imat_core == args.imat_mantle:
        raise ValueError("core and mantle imat values must differ")
    if not math.isfinite(args.box_size) or args.box_size <= 0.0:
        raise ValueError("box-size must be finite and positive")

    output_paths = [args.output, args.secondary_output] if paired else [args.output]
    if not args.force and not args.dry_run:
        existing = [path.resolve() for path in output_paths if path.exists()]
        if existing:
            raise FileExistsError(
                "output already exists: " + ", ".join(str(path) for path in existing) + "; pass --force to replace"
            )

    try:
        import woma
    except ImportError as error:
        raise RuntimeError("planet generation requires WoMa in this Python environment") from error

    np.random.seed(args.seed)
    primary = make_planet_profile(woma, "Primary", args.planet_mass_earth, args)
    primary_planet, total_mass_kg, core_mass_kg, mantle_mass_kg, ratio = primary

    secondary = None
    if paired:
        secondary = make_planet_profile(woma, "Secondary", args.secondary_mass_earth, args)

    if args.primary_particle_count is not None:
        low_count, high_count, target_core_particle_mass_kg, target_mantle_particle_mass_kg = primary_resolution(
            total_mass_kg,
            core_mass_kg,
            mantle_mass_kg,
            ratio,
            args.primary_particle_count,
        )
    else:
        target_mantle_particle_mass_kg = args.mantle_particle_mass_code * units["mass_kg"]
        target_core_particle_mass_kg = ratio * target_mantle_particle_mass_kg
        high_count = args.high_resolution_count or max(
            1, round(total_mass_kg / target_mantle_particle_mass_kg)
        )
        low_count = args.low_resolution_count or max(
            1, round(total_mass_kg / target_core_particle_mass_kg)
        )

    report_units(units)
    print("Primary model:")
    print(f"  requested mass              = {args.planet_mass_earth:.12g} Earth masses")
    if args.primary_particle_count is not None:
        print(f"  approximate merged count    = {args.primary_particle_count:,} particles")
    print(f"  CMB density ratio           = {ratio:.12g}")
    print(f"  high-resolution realization = {high_count:,} particles")
    print(f"  low-resolution realization  = {low_count:,} particles")
    print(f"  target iron particle mass   = {target_core_particle_mass_kg:.12g} kg")
    print(f"  target rock particle mass   = {target_mantle_particle_mass_kg:.12g} kg")

    if paired:
        _, _, secondary_core_mass_kg, secondary_mantle_mass_kg, secondary_ratio = secondary
        predicted_core_count, predicted_mantle_count = secondary_layer_counts(
            secondary_core_mass_kg,
            secondary_mantle_mass_kg,
            target_core_particle_mass_kg,
            target_mantle_particle_mass_kg,
        )
        print("Secondary model (predicted from primary target particle masses):")
        print(f"  requested mass              = {args.secondary_mass_earth:.12g} Earth masses")
        print(f"  CMB density ratio           = {secondary_ratio:.12g}")
        print(f"  target iron particles       = {predicted_core_count:,}")
        print(f"  target rock particles       = {predicted_mantle_count:,}")
        print(f"  approximate merged count    = {predicted_core_count + predicted_mantle_count:,}")
        print(f"  resolution trials per layer = {args.secondary_resolution_trials}")

    if args.dry_run:
        return 0

    print("Generating primary lower-resolution realization for the iron core...")
    low = particle_set(woma, primary_planet, low_count, args.neighbors)
    core = extract_material(low, primary_planet, args.woma_core_id)
    del low

    print("Generating primary higher-resolution realization for the rock mantle...")
    high = particle_set(woma, primary_planet, high_count, args.neighbors)
    mantle = extract_material(high, primary_planet, args.woma_mantle_id)
    del high

    particles = assemble_particles(
        core,
        mantle,
        core_mass_kg,
        mantle_mass_kg,
        args.eos_core_id,
        args.eos_mantle_id,
        args.imat_core,
        args.imat_mantle,
        args.neighbors,
        units,
    )
    del core, mantle
    diagnostics = validate_particles(particles, args, ratio)
    n_core = int(particles["n_core"])
    primary_core_particle_mass_kg = float(particles["m"][0]) * units["mass_kg"]
    primary_mantle_particle_mass_kg = float(particles["m"][n_core]) * units["mass_kg"]
    primary_metadata = model_metadata(
        "primary",
        args.planet_mass_earth,
        args.core_mass_fraction,
        particles,
        units,
        "self",
    )
    print(f"  selected primary core       = {n_core:,}")
    print(f"  selected primary mantle     = {len(particles['id']) - n_core:,}")
    print(f"  actual primary total        = {len(particles['id']):,}")
    print(f"  iron particle mass          = {primary_core_particle_mass_kg:.12g} kg")
    print(f"  rock particle mass          = {primary_mantle_particle_mass_kg:.12g} kg")
    for key, value in diagnostics.items():
        print(f"  {key:26s} = {value:.12g}")

    write_gizmo(args.output, particles, args, units, primary_metadata)
    print(f"Wrote primary with {len(particles['id']):,} particles to {args.output.resolve()}")
    del particles

    if not paired:
        print("planetg params: EosTableMatIds must list the mantle material ID first, e.g. '62,63'")
        return 0

    secondary_planet, secondary_total_mass_kg, secondary_core_mass_kg, secondary_mantle_mass_kg, secondary_ratio = (
        secondary
    )
    target_secondary_core_count, target_secondary_mantle_count = secondary_layer_counts(
        secondary_core_mass_kg,
        secondary_mantle_mass_kg,
        primary_core_particle_mass_kg,
        primary_mantle_particle_mass_kg,
    )
    print("Secondary model (using actual primary particle masses):")
    print(f"  target iron particles       = {target_secondary_core_count:,}")
    print(f"  target rock particles       = {target_secondary_mantle_count:,}")
    thin_to_target = args.secondary_thin_overshoot_max > 0.0
    print("Searching complete secondary lower-resolution realizations for the iron core...")
    secondary_core, secondary_low_count = best_complete_material_realization(
        woma,
        secondary_planet,
        target_secondary_core_count,
        args.core_mass_fraction,
        args.neighbors,
        args.woma_core_id,
        args.secondary_resolution_trials,
        "iron-core",
        thin_to_target=thin_to_target,
        overshoot_max=args.secondary_thin_overshoot_max,
        rng=np.random.default_rng(args.seed + 101),
    )

    print("Searching complete secondary higher-resolution realizations for the rock mantle...")
    secondary_mantle, secondary_high_count = best_complete_material_realization(
        woma,
        secondary_planet,
        target_secondary_mantle_count,
        1.0 - args.core_mass_fraction,
        args.neighbors,
        args.woma_mantle_id,
        args.secondary_resolution_trials,
        "rock-mantle",
        thin_to_target=thin_to_target,
        overshoot_max=args.secondary_thin_overshoot_max,
        rng=np.random.default_rng(args.seed + 202),
    )
    secondary_particles = assemble_particles(
        secondary_core,
        secondary_mantle,
        secondary_core_mass_kg,
        secondary_mantle_mass_kg,
        args.eos_core_id,
        args.eos_mantle_id,
        args.imat_core,
        args.imat_mantle,
        args.neighbors,
        units,
        fixed_core_particle_mass_kg=primary_core_particle_mass_kg,
        fixed_mantle_particle_mass_kg=primary_mantle_particle_mass_kg,
    )
    del secondary_core, secondary_mantle
    secondary_diagnostics = validate_particles(
        secondary_particles, args, secondary_ratio, enforce_mass_ratio=False
    )
    realized_total_kg, realized_core_kg, realized_mantle_kg = realized_masses(secondary_particles, units)
    mass_errors = {
        "total_mass_relative_error": abs(realized_total_kg / secondary_total_mass_kg - 1.0),
        "core_mass_relative_error": abs(realized_core_kg / secondary_core_mass_kg - 1.0),
        "mantle_mass_relative_error": abs(realized_mantle_kg / secondary_mantle_mass_kg - 1.0),
    }
    if max(mass_errors.values()) > args.secondary_mass_tolerance:
        raise ValueError(
            "best complete secondary realizations exceed secondary-mass-tolerance; increase "
            "--secondary-resolution-trials, increase --primary-particle-count, or relax the tolerance"
        )

    secondary_metadata = model_metadata(
        "secondary",
        args.secondary_mass_earth,
        args.core_mass_fraction,
        secondary_particles,
        units,
        f"primary:{args.output.name}",
    )
    secondary_metadata.update(
        {
            "secondaryCoreSourceResolution": secondary_low_count,
            "secondaryMantleSourceResolution": secondary_high_count,
            "secondaryResolutionTrials": args.secondary_resolution_trials,
            "secondarySelectionMethod": "overshoot-thin" if thin_to_target else "complete-resolution-search",
            "secondaryThinOvershootMax": args.secondary_thin_overshoot_max,
        }
    )
    print(f"  actual secondary total      = {len(secondary_particles['id']):,}")
    for key, value in secondary_diagnostics.items():
        print(f"  {key:26s} = {value:.12g}")
    for key, value in mass_errors.items():
        print(f"  {key:26s} = {value:.12g}")

    write_gizmo(args.secondary_output, secondary_particles, args, units, secondary_metadata)
    print(
        f"Wrote secondary with {len(secondary_particles['id']):,} particles to "
        f"{args.secondary_output.resolve()}"
    )
    print("planetg params: EosTableMatIds must list the mantle material ID first, e.g. '62,63'")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, FileExistsError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
