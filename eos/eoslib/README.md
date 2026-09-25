# M-ANEOS source tables

This directory contains the binary source tables used to build planetg and
SPH-EXA's shared `rock_planet_aneos_62_63.spheos` runtime table. Run

```bash
python scripts/make_eos_table.py
```

from the repository root. The generated `.spheos` file remains untracked
because it is larger than GitHub's ordinary 100 MB per-file limit.

## Provenance

- `MANEOStable_fosterite.in`: material 62, Forsterite-ANEOS-SLVTv1.0G1,
  71,852,328 bytes, SHA-256
  `2354ea8348bd6b0e8c5806256cc0b6a8cf3c083ef78b4380a1ce463297681d30`.
  Source release: <https://doi.org/10.5281/zenodo.3478631>.
- `MANEOStable_iron.in`: material 63, Iron-ANEOS-SLVTv0.2G1,
  71,852,328 bytes, SHA-256
  `694e76a725e1e59c337ad8d3f9edd151273ab920a904b26bb37d325a3eeb10b8`.
  Source release: <https://doi.org/10.5281/zenodo.3866507>.

Both are expanded 21 May 2025 M-ANEOS tables with 1,402 density nodes and
1,601 temperature nodes. The generated runtime table covers
`1e-25`--`1000 g/cm^3` and `1`--`1e8 K` without cropping. Its expected size is
107,763,396 bytes and its SHA-256 is
`8b73314c275b946497f2cb791ac1386d1472cc03c7f219661006c05415db0bb3`.

The iron source release includes the MIT terms reproduced in
`iron_LICENSE.md`. The forsterite Zenodo record labels its license
`other-open`; its archive does not contain a separate license file. Preserve
this README and the accompanying upstream metadata when redistributing these
tables.

When using these assets, cite the Stewart material releases above, Deng et al.
2019, and the EOSlib release used in the surrounding workflow:

Thomas Meier & Christian Reinhardt (2021), Halbarath/EOSlib v1.0.1,
<https://doi.org/10.5281/zenodo.4704950>.
