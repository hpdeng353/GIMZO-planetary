/*
 * eos_table.h -- reader + interpolator for SPH-EXA binary EOS tables (.spheos)
 *
 * This is a C port of sphexa's sph/src/tabulated_eos.cpp (format version 1-3).
 * Tables are cgs-backed: densities in g/cm^3, specific energies in erg/g,
 * pressures in dyn/cm^2 (=erg/cm^3), sound speeds in cm/s, temperatures in K,
 * specific entropies in erg/g/K.
 *
 * The table file is memory-mapped read-only (MAP_SHARED), so one 100+ MB
 * table is shared by all MPI ranks on the same node through the page cache;
 * only a small per-material index is allocated per process. If mmap is not
 * available the file is read into a single malloc'd buffer instead.
 *
 * Host must be little-endian (x86_64, AArch64-LE).
 */

#ifndef EOS_TABLE_H
#define EOS_TABLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Evaluation status codes, same names and order as sphexa's EosTableStatus
 * (invalidVolumeElementState is sphexa-GPU specific and omitted). */
typedef enum
{
    EOS_TABLE_SUCCESS = 0,
    EOS_TABLE_DENSITY_BELOW_RANGE,
    EOS_TABLE_DENSITY_ABOVE_RANGE,
    EOS_TABLE_ENERGY_BELOW_RANGE,
    EOS_TABLE_ENERGY_ABOVE_RANGE,
    EOS_TABLE_UNKNOWN_MATERIAL,
    EOS_TABLE_INVALID_INPUT,
    EOS_TABLE_ENTROPY_UNAVAILABLE,
    EOS_TABLE_INVALID_TABLE_STATE
} EosTableStatus;

/* Result of one EOS evaluation. All values are in table units (cgs) for
 * eos_table_evaluate(), or in code units for eos_table_evaluate_code(). */
typedef struct
{
    double pressure;
    double soundSpeed;
    double temperature;
    double entropy;   /* 0.0 when hasEntropy == 0 */
    int    status;    /* EosTableStatus */
    int    hasEntropy;
} EosTableState;

/* Runtime unit conversion between code units and the cgs-backed table.
 * Multiply a code-unit quantity by *ToTable before lookup; multiply a table
 * quantity by *FromTable to obtain code units. */
typedef struct
{
    double densityToTable;
    double energyToTable;
    double pressureFromTable;
    double soundSpeedFromTable;
    double entropyFromTable;
} EosTableUnits;

/* Opaque table handle. */
typedef struct EosTable EosTable;

/* Build conversion factors for a cgs-backed table, given the code unit of
 * density and of specific energy expressed in cgs. Mirrors sphexa's
 * eosTableUnitConversion():
 *   pressure unit = densityUnit * energyUnit,  speed unit = sqrt(energyUnit).
 * Returns 0 on success, -1 on invalid (non-positive/non-finite) input. */
int eos_table_units_cgs(EosTableUnits *units, double densityUnitCgs, double energyUnitCgs);

/* Load and validate a .spheos table. Returns NULL on failure and prints the
 * reason to stderr. Only one table should be loaded per process. */
EosTable *eos_table_load(const char *path);

void eos_table_free(EosTable *table);

/* Number of materials stored in the table. */
uint32_t eos_table_num_materials(const EosTable *table);

/* Material ID of the i-th stored material (i < eos_table_num_materials). */
uint32_t eos_table_material_id(const EosTable *table, uint32_t index);

/* Non-zero if the material uses a native (temperature-grid, row-dependent,
 * linear-energy) axis. Returns -1 if materialId is unknown. */
int eos_table_is_native(const EosTable *table, uint32_t materialId);

/* Non-zero if the material stores an entropy array, -1 if unknown. */
int eos_table_has_entropy(const EosTable *table, uint32_t materialId);

/* Density axis range [g/cm^3] and energy axis global range [erg/g] for a
 * material (energy range is over all rows; native grids are row-dependent).
 * Returns 0 on success, -1 if materialId unknown. */
int eos_table_bounds(const EosTable *table, uint32_t materialId,
                     double *rhoMin, double *rhoMax, double *uMin, double *uMax);

/* Evaluate P, cs, T (and S when requestEntropy != 0 and available) at
 * (rho, u) given in table units (cgs). Never fails catastrophically: on
 * out-of-range input the state is clamped and status reflects the clamp;
 * on invalid input/table state a finite fallback state is returned. */
EosTableState eos_table_evaluate(const EosTable *table, double rho, double u,
                                 uint32_t materialId, int requestEntropy);

/* Same, but takes and returns code-unit quantities using the given
 * conversion factors. Temperature is never rescaled (always K). */
EosTableState eos_table_evaluate_code(const EosTable *table, double rho, double u,
                                      uint32_t materialId, int requestEntropy,
                                      EosTableUnits units);

/* Entropy inversion at fixed density (port of sphexa's
 * TabulatedEos::invertEnergy): find the specific energy whose bilinear
 * entropy interpolant at (rho) equals entropyTarget. Density outside the
 * axis is clamped (status reflects it); entropy outside the tabulated
 * column clamps to the matching column end, mirroring the energy clamping
 * in eos_table_evaluate(). All quantities in table units (cgs); the
 * resulting u is returned via *uOut. Returns an EosTableStatus. */
int eos_table_invert_energy(const EosTable *table, double rho, double entropyTarget,
                            uint32_t materialId, double *uOut);

/* Code-unit wrapper: rho and entropyTargetCode in code units, result via
 * *uCodeOut in code specific energy. */
int eos_table_invert_energy_code(const EosTable *table, double rho, double entropyTargetCode,
                                 uint32_t materialId, EosTableUnits units, double *uCodeOut);

const char *eos_table_status_string(int status);

#ifdef __cplusplus
}
#endif

#endif /* EOS_TABLE_H */
