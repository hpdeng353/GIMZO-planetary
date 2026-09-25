# Giant-Impact Pipeline (planetg)

End-to-end recipe for giant-impact simulations with the tabulated ANEOS EOS,
ported from the sphexa toolchain. All steps except the planetg runs themselves
are Python scripts in `scripts/`; compile and run on the cluster (H100:
`bash build_h100.sh`, SYSTYPE=h100).

## 0. Prerequisites

- EOS table: a binary `.spheos` file, e.g.
  `impact-out/eos/rock_planet_aneos_62_63.spheos` (materials 62 = forsterite
  mantle, 63 = iron core). The bulk float fields are shared read-only by all
  ranks on a node through mmap and the page cache. The unaligned density and
  energy axes are copied privately, about 34.3 MiB per rank for this table.

### Generating the EOS table (the `.spheos` file is not in the git repo)

The runtime table (103 MB) exceeds git hosting limits and is therefore not
tracked. Regenerate it deterministically from the EOSlib M-ANEOS tables
(`MANEOStable_*.in`, the binary tables distributed with M-ANEOS / EOSlib —
they are not part of this repository either). EOSlib is public on GitHub:
https://github.com/Halbarath/EOSlib.git

When using these EOS tables, please cite:

- Deng et al. 2019
- Thomas Meier, & Christian Reinhardt. (2021). Halbarath/EOSlib: fixed the
  licensing (Version v.1.0.1) [Computer software]. Zenodo.
  https://doi.org/10.5281/zenodo.4704950

Build the runtime table with:

```bash
python scripts/make_eos_table.py \
    --forsterite MANEOStable_fosterite.in --iron MANEOStable_iron.in \
    -o impact-out/eos/rock_planet_aneos_62_63.spheos
```

The build is byte-deterministic. Verify your inputs and output against the
provenance checksums of this project (SHA256):

| file | bytes | SHA256 |
|---|---|---|
| MANEOStable_fosterite.in | 71852328 | 2354ea8348bd6b0e8c5806256cc0b6a8cf3c083ef78b4380a1ce463297681d30 |
| MANEOStable_iron.in | 71852328 | 694e76a725e1e59c337ad8d3f9edd151273ab920a904b26bb37d325a3eeb10b8 |
| rock_planet_aneos_62_63.spheos | 107763396 | 8b73314c275b946497f2cb791ac1386d1472cc03c7f219661006c05415db0bb3 |

Then validate the table independently of the C reader:

```bash
python scripts/test_eos_table.py impact-out/eos/rock_planet_aneos_62_63.spheos
```

and, after compiling planetg, cross-check the production C interpolation path
round-trip with the `aneos/test_eos_table.c` unit test (see `scripts/h100_phase6.sh`
for how it is built and run on the cluster).
- Runtime parameters (see `impact-out/init/noon1.params`):
  - `EosTable        eos/rock_planet_aneos_62_63.spheos` — relative paths
    resolve against `OutputDir` (which must keep its trailing `/`); absolute
    paths are used as-is.
  - `EosTableMatIds  62,63` — maps IC material index `imat` 0 → table material
  62, `imat` 1 → 63. At startup each rank prints the mapped material bounds;
    a particle whose `imat` has no mapping is a fatal error. A table missing
    entropy is also rejected because the impact build writes entropy and may
    use it for fixed-entropy relaxation.
- Compile flag: `EOS_ANEOS`. It is the giant-impact master switch and enables
  temperature, entropy, material input, clipping, and `ORELAX`; `READ_HSML`
  stays **off**.
- Riemann face states use slope-limited reconstructed density, pressure,
  internal energy and sound speed by default. Invalid/non-finite reconstructed
  states fall back as a complete side to the EOS-consistent particle-centred
  state. `EOS_REEVALUATE_FACE_STATES` can instead re-query ANEOS at every face
  for verification, but should not be used for production: it performs two
  table interpolations per interacting pair and is typically several times
  slower.

EOS range violations are clamped to the nearest table boundary and reported
with a rate-limited diagnostic. Non-finite input uses the finite table safety
fallback. An invalid table state is fatal. Valid low-pressure states retain
their table pressure and sound speed; the legacy hard-coded pressure/sound
speed floor is not applied to the `.spheos` branch.

## 1. Build the planets — `scripts/makeplanet_planetg.py`

WoMa adiabatic profiles integrated on the `.spheos` tables, realised as
equal-mass-per-layer GIZMO HDF5 particles.

```
python scripts/makeplanet_planetg.py \
    --primary-mass-earth 0.85 --core-mass-fraction 0.30 \
    --primary-output target.h5 --force
```

Key options: `--high-resolution-count` / `--low-resolution-count` (particle
counts), `--mantle-particle-mass-code`, `--surface-temperature-k`,
`--surface-pressure-pa`, `--woma-core-id 401 / --woma-mantle-id 400` (WoMa
materials), `--eos-core-id 63 / --eos-mantle-id 62` (table materials).
Paired mode (`--secondary-mass-earth`) builds the impactor in the same
resolution search.

Output datasets (PartType0): `Coordinates`, `Velocities` (zero), `Masses`,
`InternalEnergy`, `ParticleIDs` (uint64), `Materials` (uint32 imat),
`Temperature`, `SmoothingLength`. planetg reads up to `InternalEnergy` plus
`Materials` and `Temperature`; `SmoothingLength` is for inspection only
(`READ_HSML` off). Masses are centred on the origin.

Note: very small planets (~0.01 M_earth) can fail WoMa convergence; that is a
WoMa limitation, not a bug in the writer.

## 2. Relax each planet — ORELAX

Run planetg on the single-planet IC with the staged relaxation enabled via
four parameters (all zero = disabled, e.g. for the impact run):

- `RelaxTimescale` — linear drag time `dv/dt = -v/tau` while
  `Time < RelaxUntil`. The drag removes kinetic energy (not thermalised).
- `RelaxUntil` — stop all damping.
- `SphericalRelaxUntil` — while `Time` is below this, the particle velocity and
  each momentum kick are projected onto the radius vector about the origin
  (kinematically radial motion).
- `SphericalRelaxReleaseDuration` — over
  `[SphericalRelaxUntil, SphericalRelaxUntil + duration]` the tangential kick
  is restored with the cubic smoothstep `w = 3s^2 - 2s^3`. Velocity is not
  reprojected and is never scaled by `w`.
- `RelaxIsentropic` (0/1) — fixed-entropy relaxation, a port of SPH-EXA's
  `--relax-isentropic`. While `Time < RelaxUntil`, every particle's internal
  energy is reset to `u(rho, s0)` from the EOS table before the pressure
  evaluation (stored in `SphP.RelaxEntropy0`). The anchor entropy `s0` comes
  from the IC's optional `PartType0/AnchorEntropy` block when present —
  `makeplanet_planetg.py --s0-eos-table` writes it from the *smooth WoMa*
  `(rho, u)` state on the runtime's own `.spheos` table (same algorithm as
  SPH-EXA's `makeplanet_sphexa.py` `s0` field), so the pin is the identity
  on the WoMa profile and does not inherit the noisy t=0 SPH density
  evaluation (critical for thinned secondaries). ICs without the block fall
  back to adopting `s0` from the t=0 `(rho, u)` table entropy on the first
  force evaluation. This holds the constructed entropy profile exactly while
  positions/densities settle; energy conservation is intentionally violated,
  so it must stay 0 for production/impact runs. Requires an EOS table with
  entropy (both materials in `rock_planet_aneos_62_63.spheos` have it).
  Unlike SPH-EXA, `s0` is not written to snapshots: a restarted relaxation
  re-reads the IC block (RestartFlag 0) or re-adopts from the current state,
  which is equivalent as long as the state was pinned. Implementation:
  `orelax_isentropic_pin()` in `eos/eos.c`, called in the density loop
  just before `get_pressure()` (`hydro/density.c`); IC block read in
  `read_anchor_entropy()` (`read_ic.c`); entropy inversion
  `eos_table_invert_energy()` is a C port of SPH-EXA's
  `TabulatedEos::invertEnergy`.

- Compile-time flags: `EOS_ANEOS` is the single master switch. It
  auto-enables `EOS_TABULATED`, `EOS_CARRIES_TEMPERATURE`,
  `EOS_CARRIES_ENTROPY`, `MOON`, `ORELAX`, `READ_IMAT`, `CLIPPING` and
  `PREVENT_PARTICLE_MERGE_SPLIT` (see `eos/eos.h`), so `impact.conf` only
  needs `EOS_ANEOS` plus the hydro/kernel/I/O choices.

Implementation: `orelax_modify_kick()` in `run.c`, called from
`do_the_kick` (`kicks.c`) after the momentum kick `dp` is computed and before
`Vel` is updated, scaled by each particle's own kick interval. Single,
non-rotating, origin-centred planets only — never enable the spherical stage
for impact runs. Typical use: spherical stage for the first ~10–20 code time
units, smooth release over a few units, 3D damping until the planet is quiet;
then continue unrelaxed to confirm the profile holds.

## 3. Build the impact IC — `scripts/makeimpact_planetg.py`

```
python scripts/makeimpact_planetg.py relaxed_target.h5 relaxed_impactor.h5 \
    -o impact-ics.hdf5 --impact-angle-deg 45 --contact-speed-vesc 1.0 \
    --separation-code 4.0 \
    --params-template impact-out/init/noon1.params --time-max 20 --force
```

Reads two relaxed planetg snapshots, recentres each body, places them on a
parabolic (zero-energy) orbit with the requested contact speed in units of the
mutual escape speed (`--contact-speed-vesc` must be 1), impact angle measured
from the centre-to-centre line, and initial centre separation
(`--separation-code`). Radii are estimated from an enclosed-mass fraction
(`--radius-mass-fraction 0.999`) unless overridden. The system is shifted to
the centre of mass and total momentum is zeroed; `BoxSize` is set from the
particle extent (`--box-padding`). Particle IDs are renumbered to stay unique;
`Materials`, `InternalEnergy`, `Temperature` are carried over.

With `--params-template`, a patched parameter file is written
(`--params-output`, default `<output>.params`): `InitCondFile`, `OutputDir`,
`TimeMax`, unit system and `G` in code units are filled in, and the four
ORELAX parameters are forced to 0.

## 4. Run the impact

Run planetg on the impact IC with the patched params. Confirm the startup log
shows the `EOS_ANEOS: loaded .spheos table ...` summary with the expected
material bounds, and that no ORELAX damping line appears (all four
parameters are 0).

## 5. Notes and gotchas

- Code units: G_code = 6.672e-8 * UnitMass * UnitTime^2 / UnitLength^3
  (≈ 1.000035 for the noon1 units); the scripts recompute it from the units.
- The EOS table covers negative specific energies; the old `u < 2e-3` clamp is
  gone. Out-of-range lookups are counted per rank and rate-limited in the log.
- The legacy ASCII tables (rho-T interpolation, `RHOT` switch) were removed in
  Phase 5b; only `.spheos` tables are supported.
- Snapshot/IC format is plain GIZMO HDF5; the old truncated
  `noon1-ics.hdf5` (Entropy/Temperature blocks shorter than NumPart) is
  rejected by the new writers.
