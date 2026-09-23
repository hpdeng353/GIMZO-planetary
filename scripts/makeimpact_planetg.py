#!/usr/bin/env python3
"""Combine two relaxed planetg planets into giant-impact initial conditions.

Ported from sphexa's makeimpact_sphexa.py. Inputs are relaxed planetg (GIZMO /
GADGET-3 HDF5) snapshots or makeplanet_planetg.py outputs; the output is a
single GIZMO-format IC file with regenerated particle IDs, plus optionally a
parameter file patched from a template.

The trajectory is parabolic: its relative velocity at infinity is zero. It has
a contact speed of one mutual escape velocity and a default impact angle of 45
degrees. The initial centre separation is four code-length units, matching the
canonical WoMa tutorial setup when one code length is one Earth radius.

Units: planetg snapshots carry no unit attributes, so the code units are given
through --unit-length-cm / --unit-mass-g / --unit-velocity-cms (defaults match
noon1.params) and MUST be identical to the parameter file used for the run.
The two-body gravitational parameter uses G_code = GRAVITY * M * T^2 / L^3 with
GRAVITY = 6.672e-8 cgs exactly as planetg's begrun computes it (= 1.000035 in
the default units); override with --g-code if the build differs.

Materials handling: the Materials dataset (imat indices, see EosTableMatIds) is
carried through unchanged; ParticleIDs are regenerated as 1..N; Temperature is
carried when both inputs provide it (read by EOS_CARRIES_TEMPERATURE);
SmoothingLength is carried when present but only informational (READ_HSML off);
Density/Pressure/Entropy are recomputed by planetg at startup and not written.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import re

import numpy as np


# planetg allvars.h GRAVITY in cgs.
PLANETG_GRAVITY_CGS = 6.672e-8

DEFAULT_UNIT_LENGTH_CM = 6.37869e8
DEFAULT_UNIT_MASS_G = 9.56072e25
DEFAULT_UNIT_VELOCITY_CMS = 1.0e5

REQUIRED_DATASETS = ("Coordinates", "Velocities", "Masses", "InternalEnergy", "Materials")
OPTIONAL_CARRIED = ("Temperature", "SmoothingLength")


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="relaxed target planetg snapshot/IC")
    parser.add_argument("impactor", type=Path, help="relaxed impactor planetg snapshot/IC")
    parser.add_argument("-o", "--output", type=Path, default=Path("impact-planetg-ics.hdf5"))
    parser.add_argument("--impact-angle-deg", type=float, default=45.0)
    parser.add_argument(
        "--contact-speed-vesc",
        type=float,
        default=1.0,
        help="contact speed in mutual escape-speed units; must be 1 for the enforced parabolic orbit",
    )
    parser.add_argument(
        "--separation-code",
        type=float,
        default=4.0,
        help="initial centre separation in code-length units (default: 4)",
    )
    parser.add_argument("--target-radius-code", type=float, help="override the estimated target radius")
    parser.add_argument("--impactor-radius-code", type=float, help="override the estimated impactor radius")
    parser.add_argument(
        "--radius-mass-fraction",
        type=float,
        default=0.999,
        help="enclosed-mass fraction used to estimate each radius (default: 0.999)",
    )
    parser.add_argument(
        "--box-padding",
        type=float,
        default=1.25,
        help="BoxSize = 2 * padding * (particle extent + 2 h_max); informational for open boundaries (default: 1.25)",
    )
    parser.add_argument("--unit-length-cm", type=float, default=DEFAULT_UNIT_LENGTH_CM)
    parser.add_argument("--unit-mass-g", type=float, default=DEFAULT_UNIT_MASS_G)
    parser.add_argument("--unit-velocity-cms", type=float, default=DEFAULT_UNIT_VELOCITY_CMS)
    parser.add_argument(
        "--g-code",
        type=float,
        help="gravitational constant in code units (default: recomputed from the units with GRAVITY=6.672e-8 cgs)",
    )
    parser.add_argument(
        "--params-template",
        type=Path,
        help="optional planetg parameter file to clone and patch for the impact run",
    )
    parser.add_argument(
        "--params-output",
        type=Path,
        help="patched parameter file path (default: <output without extension>.params)",
    )
    parser.add_argument(
        "--time-max",
        type=float,
        default=20.0,
        help="TimeMax written into the patched parameter file, in code time units (default: 20)",
    )
    parser.add_argument("--force", action="store_true", help="atomically replace an existing output file")
    args = parser.parse_args(argv)

    if not 0.0 <= args.impact_angle_deg < 90.0:
        parser.error("--impact-angle-deg must satisfy 0 <= angle < 90")
    if not math.isfinite(args.contact_speed_vesc) or not math.isclose(
        args.contact_speed_vesc, 1.0, rel_tol=0.0, abs_tol=1.0e-12
    ):
        parser.error("--contact-speed-vesc must be 1 so that velocity at infinity is zero")
    args.contact_speed_vesc = 1.0
    if not math.isfinite(args.separation_code) or args.separation_code <= 0.0:
        parser.error("--separation-code must be finite and positive")
    if not 0.0 < args.radius_mass_fraction <= 1.0:
        parser.error("--radius-mass-fraction must satisfy 0 < fraction <= 1")
    if not math.isfinite(args.box_padding) or args.box_padding <= 1.0:
        parser.error("--box-padding must be finite and greater than one")
    for name in ("target_radius_code", "impactor_radius_code"):
        value = getattr(args, name)
        if value is not None and (not math.isfinite(value) or value <= 0.0):
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    for name in ("unit_length_cm", "unit_mass_g", "unit_velocity_cms"):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if args.g_code is not None and (not math.isfinite(args.g_code) or args.g_code <= 0.0):
        parser.error("--g-code must be finite and positive")
    if not math.isfinite(args.time_max) or args.time_max <= 0.0:
        parser.error("--time-max must be finite and positive")
    return args


def planetg_units(args) -> dict[str, float]:
    """Code-unit conversions plus the code-unit gravitational constant."""
    length_cm = args.unit_length_cm
    mass_g = args.unit_mass_g
    velocity_cms = args.unit_velocity_cms
    length_m = length_cm * 1.0e-2
    velocity_ms = velocity_cms * 1.0e-2
    time_s = length_cm / velocity_cms
    g_code = (
        args.g_code
        if args.g_code is not None
        else PLANETG_GRAVITY_CGS * mass_g * time_s * time_s / length_cm**3
    )
    return {
        "length_m": length_m,
        "velocity_ms": velocity_ms,
        "time_s": time_s,
        "mass_kg": mass_g * 1.0e-3,
        "g_code": g_code,
        "unit_length_cm": length_cm,
        "unit_mass_g": mass_g,
        "unit_velocity_cms": velocity_cms,
    }


def load_body(path, h5py):
    """Read one relaxed planetg GIZMO-format snapshot."""
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"snapshot not found: {path}")

    with h5py.File(path, "r") as source:
        if "PartType0" not in source:
            raise ValueError(f"{path} has no /PartType0 group (not a GIZMO/GADGET HDF5 file?)")
        group = source["PartType0"]
        missing = [name for name in REQUIRED_DATASETS if name not in group]
        if missing:
            raise ValueError(f"{path} is missing required datasets: {', '.join(missing)}")

        particle_count = group["Coordinates"].shape[0]
        if particle_count == 0:
            raise ValueError(f"{path} contains no particles")
        data = {}
        for name in REQUIRED_DATASETS:
            dataset = group[name]
            if dataset.shape[0] != particle_count:
                raise ValueError(
                    f"{path}:/PartType0/{name} has {dataset.shape[0]} rows but Coordinates has "
                    f"{particle_count}; truncated per-particle arrays are not accepted"
                )
            data[name] = np.asarray(dataset, dtype=np.float64)
        data["Materials"] = np.asarray(group["Materials"], dtype=np.uint32)
        for name in OPTIONAL_CARRIED:
            if name in group:
                if group[name].shape[0] != particle_count:
                    raise ValueError(
                        f"{path}:/PartType0/{name} has {group[name].shape[0]} rows but Coordinates has "
                        f"{particle_count}; truncated per-particle arrays are not accepted"
                    )
                data[name] = np.asarray(group[name], dtype=np.float64)

        header = dict(source["Header"].attrs) if "Header" in source else {}
        return {"path": path, "data": data, "header": header}


def check_units_compatible(target, impactor):
    """Cross-check unit attributes when both inputs carry them (makeplanet_planetg ICs do)."""
    for name in ("UnitLength_in_cm", "UnitMass_in_g", "UnitVelocity_in_cm_per_s"):
        values = [float(body["header"][name]) for body in (target, impactor) if name in body["header"]]
        if len(values) == 2 and not np.isclose(values[0], values[1], rtol=1e-12, atol=0.0):
            raise ValueError(f"input unit mismatch for {name}: {values[0]:.12g} vs {values[1]:.12g}")
    return None


def body_frame(body, radius_mass_fraction, radius_override):
    data = body["data"]
    mass = np.asarray(data["Masses"], dtype=np.float64)
    if not np.all(np.isfinite(mass)) or np.any(mass <= 0.0):
        raise ValueError(f"{body['path']} contains invalid particle masses")
    total_mass = float(np.sum(mass, dtype=np.float64))

    position = np.asarray(data["Coordinates"], dtype=np.float64)
    velocity = np.asarray(data["Velocities"], dtype=np.float64)
    if position.ndim != 2 or position.shape[1] != 3 or velocity.shape != position.shape:
        raise ValueError(f"{body['path']} positions/velocities do not have shape (N, 3)")
    if not np.all(np.isfinite(position)) or not np.all(np.isfinite(velocity)):
        raise ValueError(f"{body['path']} contains non-finite positions or velocities")

    center = np.sum(position * mass[:, None], axis=0, dtype=np.float64) / total_mass
    bulk_velocity = np.sum(velocity * mass[:, None], axis=0, dtype=np.float64) / total_mass
    position = position - center
    velocity = velocity - bulk_velocity

    radius_values = np.linalg.norm(position, axis=1)
    order = np.argsort(radius_values)
    enclosed_mass = np.cumsum(mass[order], dtype=np.float64)
    radius_index = min(np.searchsorted(enclosed_mass, radius_mass_fraction * total_mass), len(order) - 1)
    estimated_radius = float(radius_values[order[radius_index]])
    radius = radius_override if radius_override is not None else estimated_radius
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError(f"could not determine a positive radius for {body['path']}")

    return {
        "position": position,
        "velocity": velocity,
        "mass": total_mass,
        "radius": float(radius),
        "estimated_radius": estimated_radius,
        "maximum_radius": float(np.max(radius_values)),
        "center_offset": center,
        "bulk_velocity_offset": bulk_velocity,
    }


def woma_trajectory(args, target_frame, impactor_frame, units):
    try:
        import woma
    except ImportError as error:
        raise RuntimeError("impact generation requires WoMa in this Python environment") from error

    trajectory = getattr(woma, "impact_pos_vel_b_v_c_r", None)
    if trajectory is None:
        try:
            from woma.misc.utils import impact_pos_vel_b_v_c_r as trajectory
        except ImportError as error:
            raise RuntimeError("installed WoMa does not provide impact_pos_vel_b_v_c_r") from error

    length_m = units["length_m"]
    velocity_ms = units["velocity_ms"]
    mass_kg = units["mass_kg"]
    separation_m = args.separation_code * length_m
    target_radius_m = target_frame["radius"] * length_m
    impactor_radius_m = impactor_frame["radius"] * length_m
    if separation_m <= target_radius_m + impactor_radius_m:
        raise ValueError(
            f"initial separation {args.separation_code:.6g} is not larger than the contact separation "
            f"{target_frame['radius'] + impactor_frame['radius']:.6g} code lengths"
        )

    result = trajectory(
        b=args.impact_angle_deg,
        v_c=args.contact_speed_vesc,
        r=separation_m,
        R_t=target_radius_m,
        R_i=impactor_radius_m,
        M_t=target_frame["mass"] * mass_kg,
        M_i=impactor_frame["mass"] * mass_kg,
        units_b="B",
        units_v_c="v_esc",
        return_t=True,
    )
    relative_position = np.asarray(result[0], dtype=np.float64) / length_m
    time_to_contact = float(result[2]) / units["time_s"]
    if not np.all(np.isfinite(relative_position)) or relative_position.shape != (3,):
        raise ValueError("WoMa returned a non-finite impact trajectory")

    # Reconstruct the velocity in code units. WoMa supplies the orbit orientation, while this
    # calculation makes the two-body specific orbital energy exactly zero under the same G_code
    # that planetg computes from the parameter-file units. This avoids a small energy offset when
    # WoMa and planetg use slightly different numerical values of G.
    separation = float(np.linalg.norm(relative_position))
    if not math.isclose(separation, args.separation_code, rel_tol=1.0e-10, abs_tol=1.0e-12):
        raise ValueError(
            f"WoMa returned separation {separation:.12g}, expected {args.separation_code:.12g} code lengths"
        )

    total_mass = target_frame["mass"] + impactor_frame["mass"]
    mu = units["g_code"] * total_mass
    contact_distance = target_frame["radius"] + impactor_frame["radius"]
    if not math.isfinite(mu) or mu <= 0.0:
        raise ValueError("the two-body gravitational parameter must be finite and positive")

    radial_direction = relative_position / separation
    impact_parameter = math.sin(math.radians(args.impact_angle_deg))
    contact_speed = math.sqrt(2.0 * mu / contact_distance)
    specific_angular_momentum = contact_distance * contact_speed * impact_parameter
    tangential_speed = specific_angular_momentum / separation
    parabolic_speed = math.sqrt(2.0 * mu / separation)
    radial_speed_squared = parabolic_speed * parabolic_speed - tangential_speed * tangential_speed
    if radial_speed_squared < -1.0e-12 * parabolic_speed * parabolic_speed:
        raise ValueError("the requested impact geometry has no real inward parabolic trajectory")

    if tangential_speed == 0.0:
        tangential_direction = np.zeros(3, dtype=np.float64)
    else:
        reference_axis = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(reference_axis, radial_direction))) > 0.9:
            reference_axis = np.array([0.0, 1.0, 0.0])
        tangential_direction = np.cross(reference_axis, radial_direction)
        tangential_direction /= np.linalg.norm(tangential_direction)

    relative_velocity = (
        -math.sqrt(max(0.0, radial_speed_squared)) * radial_direction
        + tangential_speed * tangential_direction
    )
    specific_orbital_energy = 0.5 * float(np.dot(relative_velocity, relative_velocity)) - mu / separation
    energy_tolerance = 32.0 * np.finfo(np.float64).eps * max(1.0, mu / separation)
    if abs(specific_orbital_energy) > energy_tolerance:
        raise RuntimeError(
            "failed to construct a zero-energy parabolic orbit: "
            f"specific orbital energy is {specific_orbital_energy:.12g}"
        )
    return relative_position, relative_velocity, time_to_contact


def combine_particles(target, impactor, target_frame, impactor_frame, relative_position, relative_velocity):
    total_mass = target_frame["mass"] + impactor_frame["mass"]
    target_fraction = target_frame["mass"] / total_mass
    impactor_fraction = impactor_frame["mass"] / total_mass

    target_center = -impactor_fraction * relative_position
    impactor_center = target_fraction * relative_position
    target_bulk_velocity = -impactor_fraction * relative_velocity
    impactor_bulk_velocity = target_fraction * relative_velocity

    target_position = target_frame["position"] + target_center
    impactor_position = impactor_frame["position"] + impactor_center
    target_velocity = target_frame["velocity"] + target_bulk_velocity
    impactor_velocity = impactor_frame["velocity"] + impactor_bulk_velocity
    position = np.concatenate([target_position, impactor_position], axis=0)
    velocity = np.concatenate([target_velocity, impactor_velocity], axis=0)

    particles = {
        "Coordinates": position,
        "Velocities": velocity,
    }
    for name in ("Masses", "InternalEnergy", "Materials"):
        particles[name] = np.concatenate([target["data"][name], impactor["data"][name]])
    for name in OPTIONAL_CARRIED:
        if name in target["data"] and name in impactor["data"]:
            particles[name] = np.concatenate([target["data"][name], impactor["data"][name]])

    mass = np.asarray(particles["Masses"], dtype=np.float64)
    center = np.sum(position * mass[:, None], axis=0, dtype=np.float64) / np.sum(mass, dtype=np.float64)
    momentum = np.sum(velocity * mass[:, None], axis=0, dtype=np.float64)
    if np.linalg.norm(center) > 1e-10 * max(1.0, np.max(np.linalg.norm(position, axis=1))):
        raise RuntimeError(f"combined center of mass is not zero: {center}")
    if np.linalg.norm(momentum) > 1e-10 * max(1.0, np.sum(mass) * np.max(np.linalg.norm(velocity, axis=1))):
        raise RuntimeError(f"combined momentum is not zero: {momentum}")
    return particles, target_center, impactor_center, target_bulk_velocity, impactor_bulk_velocity


def write_gizmo_ic(path, particles, header_metadata, units, box_size, force, h5py):
    """Atomically write the combined GIZMO-format IC file and re-read it for validation."""
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        raise FileExistsError(f"output already exists: {path}; pass --force to replace it")
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"temporary output already exists: {temporary}")

    n = len(particles["Masses"])
    try:
        with h5py.File(temporary, "w") as output:
            header = output.create_group("Header")
            header.attrs["NumPart_ThisFile"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.int32)
            header.attrs["NumPart_Total"] = np.array([n, 0, 0, 0, 0, 0], dtype=np.uint32)
            header.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.uint32)
            header.attrs["MassTable"] = np.zeros(6, dtype=np.float64)
            header.attrs["Time"] = 0.0
            header.attrs["Redshift"] = 0.0
            header.attrs["BoxSize"] = float(box_size)
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

            header.attrs["UnitLength_in_cm"] = units["unit_length_cm"]
            header.attrs["UnitMass_in_g"] = units["unit_mass_g"]
            header.attrs["UnitVelocity_in_cm_per_s"] = units["unit_velocity_cms"]
            header.attrs["gravConstantCode"] = units["g_code"]
            for name, value in header_metadata.items():
                header.attrs[name] = value

            group = output.create_group("PartType0")
            group.create_dataset("Coordinates", data=particles["Coordinates"], dtype=np.float64)
            group.create_dataset("Velocities", data=particles["Velocities"], dtype=np.float64)
            group.create_dataset("Masses", data=np.asarray(particles["Masses"]), dtype=np.float64)
            group.create_dataset(
                "InternalEnergy", data=np.asarray(particles["InternalEnergy"]), dtype=np.float64
            )
            group.create_dataset("ParticleIDs", data=np.arange(1, n + 1, dtype=np.uint64))
            group.create_dataset("Materials", data=np.asarray(particles["Materials"]), dtype=np.uint32)
            for name in OPTIONAL_CARRIED:
                if name in particles:
                    group.create_dataset(name, data=np.asarray(particles[name]), dtype=np.float64)

        with h5py.File(temporary, "r") as check:
            group = check["PartType0"]
            expected = {"Coordinates", "Velocities", "Masses", "InternalEnergy", "ParticleIDs", "Materials"}
            missing = expected.difference(group.keys())
            wrong_size = [name for name in expected if name in group and group[name].shape[0] != n]
            header = check["Header"]
            header_ok = (
                int(header.attrs["NumPart_ThisFile"][0]) == n
                and int(header.attrs["NumPart_Total"][0]) == n
                and float(header.attrs["MassTable"][0]) == 0.0
                and float(header.attrs["Time"]) == 0.0
            )
            if missing or wrong_size or not header_ok:
                raise RuntimeError(
                    "written impact file failed validation; "
                    f"missing datasets={sorted(missing)}, wrong-size={sorted(wrong_size)}, "
                    f"header consistent={header_ok}"
                )
        os.replace(temporary, path)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise
    return path


def write_params(template: Path, params_path: Path, ic_path: Path, time_max: float, force: bool) -> Path:
    """Clone a planetg parameter file, patching the entries an impact run needs."""
    params_path = params_path.resolve()
    if params_path.exists() and not force:
        raise FileExistsError(f"parameter file already exists: {params_path}; pass --force to replace it")

    # InitCondFile wants the path without the .hdf5/.h5 suffix (planetg appends it).
    ic_reference = str(ic_path)
    for suffix in (".hdf5", ".h5"):
        if ic_reference.lower().endswith(suffix):
            ic_reference = ic_reference[: -len(suffix)]
            break

    overrides = {
        "InitCondFile": ic_reference,
        "TimeBegin": "0",
        "TimeMax": f"{time_max:.12g}",
        "RelaxTimescale": "0",
        "RelaxUntil": "0",
        "SphericalRelaxUntil": "0",
        "SphericalRelaxReleaseDuration": "0",
        "RelaxIsentropic": "0",
    }
    lines = template.read_text(encoding="utf-8").splitlines()
    written = set()
    output_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith("%"):
            parts = stripped.split()
            if parts and parts[0] in overrides:
                line = f"{parts[0]:35s}{overrides[parts[0]]}"
                written.add(parts[0])
        output_lines.append(line)
    missing = [key for key in overrides if key not in written]
    if missing:
        output_lines.append("")
        output_lines.append("% entries appended by makeimpact_planetg.py")
        for key in missing:
            output_lines.append(f"{key:35s}{overrides[key]}")
    params_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return params_path


def main(argv=None):
    args = parse_arguments(argv)
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("impact generation requires h5py in this Python environment") from error

    if args.target.resolve() == args.impactor.resolve():
        raise ValueError("target and impactor snapshots must be different files")
    if args.output.resolve() in {args.target.resolve(), args.impactor.resolve()}:
        raise ValueError("output must not overwrite either input snapshot")

    units = planetg_units(args)
    target = load_body(args.target, h5py)
    impactor = load_body(args.impactor, h5py)
    check_units_compatible(target, impactor)

    target_frame = body_frame(target, args.radius_mass_fraction, args.target_radius_code)
    impactor_frame = body_frame(impactor, args.radius_mass_fraction, args.impactor_radius_code)
    relative_position, relative_velocity, time_to_contact = woma_trajectory(
        args, target_frame, impactor_frame, units
    )
    particles, target_center, impactor_center, target_velocity, impactor_velocity = combine_particles(
        target, impactor, target_frame, impactor_frame, relative_position, relative_velocity
    )

    max_extent = float(np.max(np.abs(particles["Coordinates"])))
    max_h = float(np.max(particles["SmoothingLength"])) if "SmoothingLength" in particles else 0.0
    box_size = 2.0 * args.box_padding * (max_extent + 2.0 * max_h)

    total_mass = target_frame["mass"] + impactor_frame["mass"]
    contact_distance = target_frame["radius"] + impactor_frame["radius"]
    mutual_escape_speed = math.sqrt(2.0 * units["g_code"] * total_mass / contact_distance)
    initial_separation = float(np.linalg.norm(relative_position))
    initial_relative_speed = float(np.linalg.norm(relative_velocity))
    specific_orbital_energy = (
        0.5 * initial_relative_speed * initial_relative_speed
        - units["g_code"] * total_mass / initial_separation
    )
    header_metadata = {
        "planetModelRole": "impact",
        "impactAngleDeg": args.impact_angle_deg,
        "impactParameter": math.sin(math.radians(args.impact_angle_deg)),
        "impactOrbitType": "parabolic",
        "impactVelocityAtInfinityCode": 0.0,
        "impactSpecificOrbitalEnergyCode": specific_orbital_energy,
        "mutualEscapeVelocityCode": mutual_escape_speed,
        "mutualEscapeVelocityMPerS": mutual_escape_speed * units["velocity_ms"],
        "initialRelativeSpeedCode": initial_relative_speed,
        "initialParabolicSpeedCode": math.sqrt(2.0 * units["g_code"] * total_mass / initial_separation),
        "initialSeparationCode": initial_separation,
        "estimatedTimeToContactCode": time_to_contact,
        "targetMassCode": target_frame["mass"],
        "impactorMassCode": impactor_frame["mass"],
        "targetRadiusCode": target_frame["radius"],
        "impactorRadiusCode": impactor_frame["radius"],
        "targetParticleCount": np.uint64(len(target["data"]["Masses"])),
        "impactorParticleCount": np.uint64(len(impactor["data"]["Masses"])),
        "sourceTarget": str(target["path"]),
        "sourceImpactor": str(impactor["path"]),
    }

    output = write_gizmo_ic(args.output, particles, header_metadata, units, box_size, args.force, h5py)

    params_path = None
    if args.params_template is not None:
        if not args.params_template.is_file():
            raise FileNotFoundError(f"parameter template not found: {args.params_template}")
        params_path = args.params_output or output.with_suffix("").with_suffix(".params")
        params_path = write_params(args.params_template, params_path, output, args.time_max, args.force,)

    print("planetg giant-impact initial conditions")
    print(f"  target:                   {target['path']}")
    print(f"  impactor:                 {impactor['path']}")
    print(f"  particles:                {len(particles['Masses']):,}")
    print(f"  radii (target, impactor): {target_frame['radius']:.8g}, {impactor_frame['radius']:.8g} code length")
    print(f"  target COM drift removed: {target_frame['center_offset']}, bulk {target_frame['bulk_velocity_offset']}")
    print(f"  impactor COM drift removed: {impactor_frame['center_offset']}, bulk {impactor_frame['bulk_velocity_offset']}")
    print(f"  initial separation:       {initial_separation:.8g} code length")
    print(f"  impact angle:             {args.impact_angle_deg:.8g} deg")
    print("  orbit:                    parabolic (velocity at infinity = 0)")
    print(f"  G (code units):           {units['g_code']:.8g}")
    print(f"  contact mutual v_escape:  {mutual_escape_speed:.8g} code velocity")
    print(f"  initial relative speed:   {initial_relative_speed:.8g} code velocity")
    print(f"  estimated contact time:   {time_to_contact:.8g} code time")
    print(f"  target center/velocity:   {target_center}, {target_velocity}")
    print(f"  impactor center/velocity: {impactor_center}, {impactor_velocity}")
    print(f"Wrote {output}")
    if params_path is not None:
        print(f"Wrote {params_path}")
    print("Reminder: run with all four MOONRELAX parameters at 0; Materials indices map")
    print("through EosTableMatIds (e.g. '62,63' = mantle,core) and must match the inputs.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=__import__("sys").stderr)
        raise SystemExit(2)
