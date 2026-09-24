#!/usr/bin/env python3

"""Build a planetg/SPH-EXA .spheos EOS table from EOSlib MANEOStable inputs.

Deterministic: the same MANEOStable inputs always produce a byte-identical
table (fixed field order, FNV-1a checksum in the header, self-validated before
the atomic rename). See docs/impact_pipeline.md section 0 for the provenance
checksums of the tables used by this project.

Usage:
    python scripts/make_eos_table.py         --forsterite MANEOStable_fosterite.in --iron MANEOStable_iron.in         -o impact-out/eos/rock_planet_aneos_62_63.spheos
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import struct
import tempfile

import numpy as np


MAGIC = b"SPXEOST1"
FORMAT_VERSION = 3
ENDIAN_MARKER = 0x01020304
MATERIALS = ((62, "forsterite"), (63, "iron"))


def fnv1a64(data: bytes) -> int:
    value = 14695981039346656037
    for byte in data:
        value ^= byte
        value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def read_maneos(path: Path) -> dict[str, object]:
    size = path.stat().st_size
    with path.open("rb") as stream:
        header = stream.read(16)
        if len(header) != 16:
            raise ValueError(f"{path}: truncated MANEOStable header")
        rho0, nrho, ntemp = struct.unpack("<dii", header)
        if nrho < 2 or ntemp < 2:
            raise ValueError(f"{path}: invalid grid dimensions {nrho} x {ntemp}")
        rho = np.fromfile(stream, dtype="<f8", count=nrho)
        temperature = np.fromfile(stream, dtype="<f8", count=ntemp)
    if rho.size != nrho or temperature.size != ntemp:
        raise ValueError(f"{path}: truncated grid axes")
    if not np.all(np.isfinite(rho)) or not np.all(np.diff(rho) > 0) or np.any(rho <= 0):
        raise ValueError(f"{path}: density axis must be finite, positive, and strictly increasing")
    if not np.all(np.isfinite(temperature)) or not np.all(np.diff(temperature) > 0) or np.any(temperature <= 0):
        raise ValueError(f"{path}: temperature axis must be finite, positive, and strictly increasing")

    axis_end = 16 + 8 * (nrho + ntemp)
    values_per_field = nrho * ntemp
    required_size = axis_end + 4 * values_per_field * 8
    if size not in (required_size, required_size + 1024):
        raise ValueError(
            f"{path}: size {size} does not match a {nrho} x {ntemp} MANEOStable table "
            f"({required_size} bytes, optionally plus a 1024-byte material string)"
        )

    fields = {}
    for field_index, name in enumerate(("pressure", "energy", "entropy", "sound_speed")):
        field = np.memmap(
            path,
            dtype="<f8",
            mode="r",
            offset=axis_end + field_index * values_per_field * 8,
            shape=(ntemp, nrho),
        )
        if not np.all(np.isfinite(field)):
            raise ValueError(f"{path}: {name} contains non-finite values")
        fields[name] = field
    if np.any(fields["pressure"] <= 0) or np.any(fields["sound_speed"] <= 0):
        raise ValueError(f"{path}: pressure and sound speed must be strictly positive")

    return {
        "path": path,
        "rho0": rho0,
        "rho": rho,
        "temperature": temperature,
        **fields,
    }


def append_array(payload: bytearray, values: np.ndarray, dtype: str) -> None:
    contiguous = np.ascontiguousarray(values, dtype=dtype)
    payload.extend(contiguous.tobytes(order="C"))


def append_material(payload: bytearray, material_id: int, table: dict[str, object]) -> None:
    rho = table["rho"]
    temperature = table["temperature"]
    nrho = int(rho.size)
    ntemp = int(temperature.size)
    # Flags: logarithmic positive outputs, row-dependent energy, native-T grid, entropy present.
    payload.extend(struct.pack("<IIII", material_id, nrho, ntemp, 15))
    append_array(payload, np.log(rho), "<f8")

    # EOSlib arrays are [temperature][density]; SPH-EXA stores [density][temperature].
    append_array(payload, np.asarray(table["energy"]).T, "<f8")
    append_array(payload, np.log(np.asarray(table["pressure"]).T), "<f4")
    append_array(payload, np.log(np.asarray(table["sound_speed"]).T), "<f4")
    log_temperature = np.broadcast_to(np.log(temperature), (nrho, ntemp))
    append_array(payload, log_temperature, "<f4")
    append_array(payload, np.asarray(table["entropy"]).T, "<f4")


def validate_spheos(path: Path, expected: list[tuple[int, int, int]]) -> None:
    data = path.read_bytes()
    header_size = 32
    if len(data) < header_size:
        raise ValueError("generated SPH-EXA table has a truncated header")
    magic, version, endian, payload_size, checksum = struct.unpack("<8sIIQQ", data[:header_size])
    payload = data[header_size:]
    if magic != MAGIC or version != FORMAT_VERSION or endian != ENDIAN_MARKER:
        raise ValueError("generated SPH-EXA table header is invalid")
    if payload_size != len(payload) or checksum != fnv1a64(payload):
        raise ValueError("generated SPH-EXA table size or checksum is invalid")

    offset = 0
    (count,) = struct.unpack_from("<I", payload, offset)
    offset += 4
    if count != len(expected):
        raise ValueError("generated SPH-EXA table has the wrong material count")
    for expected_id, expected_nrho, expected_ntemp in expected:
        material_id, nrho, ntemp, flags = struct.unpack_from("<IIII", payload, offset)
        offset += 16
        if (material_id, nrho, ntemp, flags) != (expected_id, expected_nrho, expected_ntemp, 15):
            raise ValueError("generated SPH-EXA material metadata is invalid")
        offset += nrho * 8
        offset += nrho * ntemp * 8
        offset += 4 * nrho * ntemp * 4
    if offset != len(payload):
        raise ValueError("generated SPH-EXA table has trailing or missing payload data")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--forsterite", required=True, type=Path, help="EOSlib MANEOStable forsterite input")
    parser.add_argument("--iron", required=True, type=Path, help="EOSlib MANEOStable iron input")
    parser.add_argument("-o", "--output", required=True, type=Path, help="output .spheos table")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = {"forsterite": args.forsterite.resolve(), "iron": args.iron.resolve()}
    tables = {name: read_maneos(inputs[name]) for _, name in MATERIALS}

    payload = bytearray(struct.pack("<I", len(MATERIALS)))
    expected = []
    for material_id, name in MATERIALS:
        table = tables[name]
        append_material(payload, material_id, table)
        expected.append((material_id, int(table["rho"].size), int(table["temperature"].size)))

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    header = struct.pack(
        "<8sIIQQ", MAGIC, FORMAT_VERSION, ENDIAN_MARKER, len(payload), fnv1a64(payload)
    )
    descriptor, temporary_name = tempfile.mkstemp(prefix=output.name + ".", suffix=".tmp", dir=output.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(header)
            stream.write(payload)
        validate_spheos(temporary, expected)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    print(f"Wrote {output} ({output.stat().st_size} bytes)")
    print(f"SHA-256 {hashlib.sha256(output.read_bytes()).hexdigest()}")
    for material_id, name in MATERIALS:
        table = tables[name]
        print(
            f"material {material_id} {name}: {table['rho'].size} x {table['temperature'].size}, "
            f"rho={table['rho'][0]:.8g}..{table['rho'][-1]:.8g} g/cm^3, "
            f"T={table['temperature'][0]:.8g}..{table['temperature'][-1]:.8g} K"
        )


if __name__ == "__main__":
    main()
