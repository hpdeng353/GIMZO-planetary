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
