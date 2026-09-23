/*
 * test_eos_table.c -- standalone smoke test for the .spheos reader.
 *
 * Usage:  test_eos_table <path-to-table.spheos> [materialId]
 *
 * Loads the table, prints format/grid metadata, evaluates a set of sample
 * (rho, u) points, then scans the whole grid of every material checking
 * that all interpolated states are finite and positive. Exit status is 0
 * when all checks pass.
 *
 * Build:  cc -O2 -o test_eos_table eos_table.c test_eos_table.c -lm
 */

#include "eos_table.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_point(const EosTable *t, uint32_t materialId, double rho, double u)
{
    EosTableState s = eos_table_evaluate(t, rho, u, materialId, 1);
    printf("  rho=%12.5e g/cm^3  u=%12.5e erg/g -> P=%12.5e  cs=%12.5e  T=%12.5e  S=%12.5e  status=%s%s\n",
           rho, u, s.pressure, s.soundSpeed, s.temperature, s.entropy,
           eos_table_status_string(s.status), s.hasEntropy ? "" : " (no entropy)");
    int ok = isfinite(s.pressure) && isfinite(s.soundSpeed) && isfinite(s.temperature) &&
             isfinite(s.entropy) && s.soundSpeed > 0.0 && s.temperature > 0.0;
    if (!ok) printf("  ^^ FAIL: non-finite or non-positive state\n");
    return ok;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: %s <table.spheos> [materialId]\n", argv[0]);
        return 2;
    }

    EosTable *table = eos_table_load(argv[1]);
    if (!table)
    {
        fprintf(stderr, "failed to load %s\n", argv[1]);
        return 1;
    }

    uint32_t numMaterials = eos_table_num_materials(table);
    printf("loaded %s: %u material(s)\n", argv[1], numMaterials);

    int failures = 0;
    for (uint32_t i = 0; i < numMaterials; ++i)
    {
        uint32_t id = eos_table_material_id(table, i);
        double rhoMin, rhoMax, uMin, uMax;
        eos_table_bounds(table, id, &rhoMin, &rhoMax, &uMin, &uMax);
        printf("material %u: native=%d entropy=%d  rho in [%e, %e] g/cm^3  u in [%e, %e] erg/g\n",
               id, eos_table_is_native(table, id), eos_table_has_entropy(table, id),
               rhoMin, rhoMax, uMin, uMax);
    }

    /* Sample points: mid-range and near the corners for the requested or
     * first material. */
    uint32_t id = (argc > 2) ? (uint32_t)strtoul(argv[2], NULL, 10)
                             : eos_table_material_id(table, 0);
    double rhoMin, rhoMax, uMin, uMax;
    if (eos_table_bounds(table, id, &rhoMin, &rhoMax, &uMin, &uMax) != 0)
    {
        fprintf(stderr, "material %u not present in table\n", id);
        eos_table_free(table);
        return 1;
    }
    printf("sample evaluations for material %u:\n", id);
    double rhoMid = exp(0.5 * (log(rhoMin) + log(rhoMax)));
    double uMid   = (uMin > 0.0) ? exp(0.5 * (log(uMin) + log(uMax))) : 0.5 * (uMin + uMax);
    failures += !check_point(table, id, rhoMid, uMid);
    failures += !check_point(table, id, rhoMin, uMin);
    failures += !check_point(table, id, rhoMax, uMax);
    failures += !check_point(table, id, rhoMid, uMax);
    failures += !check_point(table, id, rhoMid, uMin);
    failures += !check_point(table, id, 0.1 * rhoMin, uMid);   /* below range: expect clamp */
    failures += !check_point(table, id, rhoMid, 100.0 * uMax); /* above range: expect clamp */

    /* Full-grid scan: evaluate at every stored node (plus a mid-cell point)
     * and require finite, positive P, cs, T. */
    printf("full-grid scan over all materials ...\n");
    fflush(stdout);
    for (uint32_t i = 0; i < numMaterials; ++i)
    {
        uint32_t mid2 = eos_table_material_id(table, i);
        double rMin, rMax, eMin, eMax;
        eos_table_bounds(table, mid2, &rMin, &rMax, &eMin, &eMax);
        const int NR = 40, NU = 40;
        for (int ir = 0; ir < NR; ++ir)
        {
            double rho = exp(log(rMin) + (log(rMax) - log(rMin)) * ir / (NR - 1));
            for (int iu = 0; iu < NU; ++iu)
            {
                double u = (eMin > 0.0)
                               ? exp(log(eMin) + (log(eMax) - log(eMin)) * iu / (NU - 1))
                               : eMin + (eMax - eMin) * iu / (NU - 1);
                EosTableState s = eos_table_evaluate(table, rho, u, mid2, 1);
                if (!isfinite(s.pressure) || !(s.pressure > 0.0) ||
                    !isfinite(s.soundSpeed) || !(s.soundSpeed > 0.0) ||
                    !isfinite(s.temperature) || !(s.temperature > 0.0) ||
                    !isfinite(s.entropy))
                {
                    printf("  FAIL material %u at rho=%e u=%e: P=%e cs=%e T=%e S=%e status=%s\n",
                           mid2, rho, u, s.pressure, s.soundSpeed, s.temperature, s.entropy,
                           eos_table_status_string(s.status));
                    ++failures;
                    if (failures > 20) { printf("too many failures, aborting scan\n"); goto done; }
                }
            }
        }
        printf("  material %u: %d x %d grid OK\n", mid2, NR, NU);
    }

    /* Entropy-inversion round trip: sample (rho,u) strictly inside the grid,
     * take the interpolated entropy, and require eos_table_invert_energy to
     * recover u. Both paths interpolate the same stored float32 nodes in
     * double, so the round trip should be near machine precision; the loose
     * tolerance guards against logic errors (wrong row/stride/clamp), not
     * float32 noise. */
    printf("entropy inversion round trip ...\n");
    for (uint32_t i = 0; i < numMaterials; ++i)
    {
        uint32_t mid2 = eos_table_material_id(table, i);
        if (!eos_table_has_entropy(table, mid2))
        {
            double uDummy = 0.0;
            int st = eos_table_invert_energy(table, 1.0, 1.0, mid2, &uDummy);
            if (st != EOS_TABLE_ENTROPY_UNAVAILABLE)
            {
                printf("  FAIL material %u: inversion without entropy returned %s\n",
                       mid2, eos_table_status_string(st));
                ++failures;
            }
            continue;
        }
        double rMin, rMax, eMin, eMax;
        eos_table_bounds(table, mid2, &rMin, &rMax, &eMin, &eMax);
        int native = eos_table_is_native(table, mid2);
        const int NR = 20, NU = 20;
        double worstRel = 0.0;
        int skipped = 0;
        for (int ir = 0; ir < NR; ++ir)
        {
            double rho = exp(log(rMin) + (log(rMax) - log(rMin)) * (ir + 0.5) / NR);
            /* The hot low-density vapor corner of these tables has a
             * non-monotonic entropy column S(u): the same entropy occurs at
             * several energies, so the inverse (which, like sphexa's
             * invertEnergy, returns the lowest bracketing u) is genuinely
             * ambiguous there. S(u) at fixed rho is piecewise linear between
             * table nodes, so a dense profile pinpoints the wiggles; only
             * require the round trip where S(u) exceeds every value below. */
            const int ND = 1024;
            double uDense[1024], sDense[1024];
            for (int id = 0; id < ND; ++id)
            {
                uDense[id] = native ? eMin + (eMax - eMin) * id / (ND - 1)
                                    : exp(log(eMin) + (log(eMax) - log(eMin)) * id / (ND - 1));
                EosTableState sd = eos_table_evaluate(table, rho, uDense[id], mid2, 1);
                sDense[id] = (sd.hasEntropy && isfinite(sd.entropy)) ? sd.entropy : -DBL_MAX;
            }
            for (int iu = 0; iu < NU; ++iu)
            {
                double u = native ? eMin + (eMax - eMin) * (iu + 0.5) / NU
                                  : exp(log(eMin) + (log(eMax) - log(eMin)) * (iu + 0.5) / NU);
                EosTableState s = eos_table_evaluate(table, rho, u, mid2, 1);
                if (!s.hasEntropy || s.status != EOS_TABLE_SUCCESS)
                    continue;
                double sMaxBelow = -DBL_MAX;
                for (int id = 0; id < ND && uDense[id] < u; ++id)
                    if (sDense[id] > sMaxBelow) sMaxBelow = sDense[id];
                if (!(s.entropy > sMaxBelow + 1e-9 * fabs(s.entropy)))
                {
                    ++skipped;
                    continue;
                }
                double uBack = 0.0;
                int st = eos_table_invert_energy(table, rho, s.entropy, mid2, &uBack);
                if (st != EOS_TABLE_SUCCESS)
                {
                    printf("  FAIL material %u at rho=%e u=%e: inversion status=%s\n",
                           mid2, rho, u, eos_table_status_string(st));
                    ++failures;
                    continue;
                }
                double scale = fabs(u);
                double absFloor = 1e-9 * (fabs(eMin) + fabs(eMax));
                if (scale < absFloor) scale = absFloor;
                double rel = fabs(uBack - u) / scale;
                if (rel > worstRel) worstRel = rel;
                if (rel > 1e-4)
                {
                    printf("  FAIL material %u at rho=%e u=%e: round trip u_back=%e (rel err %g)\n",
                           mid2, rho, u, uBack, rel);
                    ++failures;
                    if (failures > 20) { printf("too many failures, aborting scan\n"); goto done; }
                }
            }
        }
        printf("  material %u: round trip worst relative error %.3e (%d non-unique-entropy points skipped)\n",
               mid2, worstRel, skipped);
    }

    /* code-unit wrapper: evaluate_code -> invert_energy_code round trip */
    {
        EosTableUnits units;
        if (eos_table_units_cgs(&units, 0.36838, 1.0e10) == 0)
        {
            for (uint32_t i = 0; i < numMaterials; ++i)
            {
                uint32_t mid2 = eos_table_material_id(table, i);
                if (!eos_table_has_entropy(table, mid2)) continue;
                double rMin, rMax, eMin, eMax;
                eos_table_bounds(table, mid2, &rMin, &rMax, &eMin, &eMax);
                double rhoC = exp(0.5 * (log(rMin) + log(rMax))) / units.densityToTable;
                double uC = (eos_table_is_native(table, mid2) ? 0.5 * (eMin + eMax)
                            : exp(0.5 * (log(eMin) + log(eMax)))) / units.energyToTable;
                EosTableState s = eos_table_evaluate_code(table, rhoC, uC, mid2, 1, units);
                double uBack = 0.0;
                int st = eos_table_invert_energy_code(table, rhoC, s.entropy, mid2, units, &uBack);
                double rel = fabs(uBack - uC) / fabs(uC);
                printf("  material %u code-units wrapper: rel err %.3e (status %s)\n",
                       mid2, rel, eos_table_status_string(st));
                if (st != EOS_TABLE_SUCCESS || !(rel < 1e-6) || !isfinite(uBack))
                    ++failures;
            }
        }
    }

done:
    eos_table_free(table);
    if (failures)
    {
        printf("RESULT: %d failure(s)\n", failures);
        return 1;
    }
    printf("RESULT: all checks passed\n");
    return 0;
}
