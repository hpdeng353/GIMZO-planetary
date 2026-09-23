/*
 * eos_table.c -- reader + interpolator for SPH-EXA binary EOS tables (.spheos)
 *
 * Faithful C port of sphexa's sph/src/tabulated_eos.cpp (format v1-3),
 * including the native (temperature-grid) energy-axis location algorithm
 * with all of its fallback branches. Do not "simplify" the native lookups:
 * they are deliberately defensive against non-monotonic T(rho,u) columns
 * and plateaus in ANEOS tables.
 *
 * Binary layout (all little-endian):
 *   header (32 bytes):
 *     char     magic[8]      = "SPXEOST1"
 *     uint32_t version       (1..3 supported)
 *     uint32_t endianMarker  = 0x01020304
 *     uint64_t payloadSize
 *     uint64_t checksum      = FNV-1a-64 over payload
 *   payload:
 *     uint32_t numMaterials
 *     per material:
 *       uint32_t materialId, nRho, nU, flags
 *         flags bit0 = 1 (required)
 *         flags bit1 = rowDependentU       (v2+)
 *         flags bit2 = nativeTemperatureGrid (v2+)
 *         flags bit3 = hasEntropy          (v3+)
 *       double logRho[nRho]
 *       double energyAxis[ rowDependentU ? nRho*nU : nU ]
 *         (log-u values, or LINEAR u values -- possibly signed -- when native)
 *       float  logPressure[nRho*nU]     (log-encoded)
 *       float  logSoundSpeed[nRho*nU]   (log-encoded)
 *       float  logTemperature[nRho*nU]  (log-encoded)
 *       float  entropy[nRho*nU]         (linear, may be negative; v3+ only)
 */

#include "eos_table.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#include <io.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#define EOS_TABLE_MAGIC "SPXEOST1"
#define EOS_TABLE_MAGIC_SIZE 8
#define EOS_TABLE_HEADER_SIZE 32 /* 8 + 2*4 + 2*8 */
#define EOS_TABLE_ENDIAN_MARKER 0x01020304u
#define EOS_TABLE_FORMAT_VERSION 3u
#define EOS_TABLE_MAX_MATERIALS 65536u

#define FLAG_ROW_DEPENDENT_U 0x2u
#define FLAG_NATIVE_TGRID    0x4u
#define FLAG_HAS_ENTROPY     0x8u

/* ------------------------------------------------------------------ */
/* internal material view                                              */
/* ------------------------------------------------------------------ */

typedef struct
{
    uint32_t materialId;
    uint32_t nRho;
    uint32_t nU;
    int      rowDependentU;
    int      nativeGrid;
    int      hasEntropy;

    double  *logRho;     /* owned copy, nRho entries */
    double  *energyAxis; /* owned copy, nU or nRho*nU entries */
    uint32_t uStride;    /* nU when rowDependentU, else 0 */

    const float *logPressure;    /* points into mapped image, nRho*nU */
    const float *logSoundSpeed;  /* ditto */
    const float *logTemperature; /* ditto */
    const float *entropy;        /* ditto, NULL when absent */

    double minSoundSpeed; /* linear values, for the invalid-state fallback */
    double minTemperature;
} MaterialView;

struct EosTable
{
    void    *image;     /* mmap base or malloc'd buffer */
    size_t   imageSize;
    int      mapped;    /* non-zero when image is an mmap */
    uint32_t version;

    uint32_t       numMaterials;
    MaterialView  *materials;   /* numMaterials entries */
    int           *materialIndex; /* materialId -> index, -1 when absent */
    uint32_t       indexSize;     /* maxMaterialId + 1 */
};

/* ------------------------------------------------------------------ */
/* small utilities                                                     */
/* ------------------------------------------------------------------ */

static uint64_t fnv1a64(const unsigned char *bytes, size_t count)
{
    uint64_t hash = 14695981039346656037ull;
    for (size_t i = 0; i < count; ++i)
    {
        hash ^= bytes[i];
        hash *= 1099511628211ull;
    }
    return hash;
}

static double lerpd(double a, double b, double t)
{
    return a + t * (b - a);
}

static double clampd(double v, double lo, double hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

static int host_is_little_endian(void)
{
    const uint16_t one = 1;
    return *((const unsigned char *)&one) == 1;
}

/* ------------------------------------------------------------------ */
/* loading                                                             */
/* ------------------------------------------------------------------ */

static void eos_table_free_members(EosTable *table)
{
    if (!table) return;
    if (table->materials)
    {
        for (uint32_t i = 0; i < table->numMaterials; ++i)
        {
            free(table->materials[i].logRho);
            free(table->materials[i].energyAxis);
        }
        free(table->materials);
    }
    free(table->materialIndex);
    if (table->image)
    {
#if !defined(_WIN32)
        if (table->mapped) munmap(table->image, table->imageSize);
        else
#endif
            free(table->image);
    }
    free(table);
}

void eos_table_free(EosTable *table)
{
    eos_table_free_members(table);
}

/* Map (or read) the whole file. Returns 0 on success. */
static int load_image(const char *path, void **imageOut, size_t *sizeOut, int *mappedOut)
{
    *imageOut  = NULL;
    *sizeOut   = 0;
    *mappedOut = 0;

#if !defined(_WIN32)
    int fd = open(path, O_RDONLY);
    if (fd >= 0)
    {
        struct stat st;
        if (fstat(fd, &st) == 0 && st.st_size > 0)
        {
            void *map = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_SHARED, fd, 0);
            if (map != MAP_FAILED)
            {
                close(fd);
                *imageOut  = map;
                *sizeOut   = (size_t)st.st_size;
                *mappedOut = 1;
                return 0;
            }
        }
        close(fd);
        /* fall through to fread */
    }
#endif

    FILE *file = fopen(path, "rb");
    if (!file) return -1;
    if (fseek(file, 0, SEEK_END) != 0) { fclose(file); return -1; }
    long length = ftell(file);
    if (length <= 0) { fclose(file); return -1; }
    rewind(file);
    void *buffer = malloc((size_t)length);
    if (!buffer) { fclose(file); return -1; }
    if (fread(buffer, 1, (size_t)length, file) != (size_t)length)
    {
        free(buffer);
        fclose(file);
        return -1;
    }
    fclose(file);
    *imageOut  = buffer;
    *sizeOut   = (size_t)length;
    return 0;
}

/* Bounds-checked little-endian cursor over the payload. Host is required to
 * be little-endian, so scalar reads are plain loads. */
typedef struct
{
    const unsigned char *data;
    size_t               size;
    size_t               pos;
    int                  error;
} ByteCursor;

static uint32_t cursor_u32(ByteCursor *c)
{
    if (c->pos + 4 > c->size) { c->error = 1; return 0; }
    uint32_t value;
    memcpy(&value, c->data + c->pos, 4);
    c->pos += 4;
    return value;
}

static uint64_t cursor_u64(ByteCursor *c)
{
    if (c->pos + 8 > c->size) { c->error = 1; return 0; }
    uint64_t value;
    memcpy(&value, c->data + c->pos, 8);
    c->pos += 8;
    return value;
}

/* Copy n doubles out of the cursor into a fresh malloc'd array. */
static double *cursor_doubles(ByteCursor *c, size_t n)
{
    if (n == 0 || c->pos + n * sizeof(double) > c->size || n > SIZE_MAX / sizeof(double))
    {
        c->error = 1;
        return NULL;
    }
    double *out = (double *)malloc(n * sizeof(double));
    if (!out) { c->error = 1; return NULL; }
    memcpy(out, c->data + c->pos, n * sizeof(double));
    c->pos += n * sizeof(double);
    return out;
}

/* Return a pointer to n floats inside the mapped image (zero-copy).
 * 4-byte alignment is guaranteed by the layout: payload starts at offset 32,
 * descriptors are 16 bytes, axes are 8-byte doubles. */
static const float *cursor_floats(ByteCursor *c, size_t n)
{
    if (n == 0 || c->pos + n * sizeof(float) > c->size || n > SIZE_MAX / sizeof(float))
    {
        c->error = 1;
        return NULL;
    }
    const float *out = (const float *)(const void *)(c->data + c->pos);
    c->pos += n * sizeof(float);
    return out;
}

static int axis_is_strictly_increasing(const double *axis, uint32_t n)
{
    for (uint32_t i = 1; i < n; ++i)
    {
        if (!(axis[i] > axis[i - 1])) return 0;
    }
    return 1;
}

EosTable *eos_table_load(const char *path)
{
    if (!host_is_little_endian())
    {
        fprintf(stderr, "eos_table: .spheos tables require a little-endian host\n");
        return NULL;
    }

    EosTable *table = (EosTable *)calloc(1, sizeof(EosTable));
    if (!table) return NULL;

    if (load_image(path, &table->image, &table->imageSize, &table->mapped) != 0)
    {
        fprintf(stderr, "eos_table: cannot open or read %s\n", path);
        eos_table_free_members(table);
        return NULL;
    }

    if (table->imageSize < EOS_TABLE_HEADER_SIZE)
    {
        fprintf(stderr, "eos_table: %s header is truncated\n", path);
        goto fail;
    }

    {
        ByteCursor header = {(const unsigned char *)table->image, EOS_TABLE_HEADER_SIZE, 0, 0};
        char magic[EOS_TABLE_MAGIC_SIZE];
        memcpy(magic, header.data, EOS_TABLE_MAGIC_SIZE);
        header.pos = EOS_TABLE_MAGIC_SIZE;
        if (memcmp(magic, EOS_TABLE_MAGIC, EOS_TABLE_MAGIC_SIZE) != 0)
        {
            fprintf(stderr, "eos_table: %s has invalid magic\n", path);
            goto fail;
        }
        uint32_t version     = cursor_u32(&header);
        uint32_t endian      = cursor_u32(&header);
        uint64_t payloadSize = cursor_u64(&header);
        uint64_t checksum    = cursor_u64(&header);
        if (header.error || version < 1u || version > EOS_TABLE_FORMAT_VERSION)
        {
            fprintf(stderr, "eos_table: %s has unsupported format version %u\n", path, version);
            goto fail;
        }
        if (endian != EOS_TABLE_ENDIAN_MARKER)
        {
            fprintf(stderr, "eos_table: %s byte order is incompatible\n", path);
            goto fail;
        }
        if (payloadSize != (uint64_t)table->imageSize - EOS_TABLE_HEADER_SIZE)
        {
            fprintf(stderr, "eos_table: %s payload size mismatch\n", path);
            goto fail;
        }
        const unsigned char *payload = (const unsigned char *)table->image + EOS_TABLE_HEADER_SIZE;
        if (fnv1a64(payload, (size_t)payloadSize) != checksum)
        {
            fprintf(stderr, "eos_table: %s checksum mismatch\n", path);
            goto fail;
        }
        table->version = version;

        ByteCursor cursor = {payload, (size_t)payloadSize, 0, 0};
        uint32_t numMaterials = cursor_u32(&cursor);
        if (cursor.error || numMaterials == 0 || numMaterials > EOS_TABLE_MAX_MATERIALS)
        {
            fprintf(stderr, "eos_table: %s material count is invalid\n", path);
            goto fail;
        }
        table->numMaterials = numMaterials;
        table->materials    = (MaterialView *)calloc(numMaterials, sizeof(MaterialView));
        if (!table->materials) goto fail;

        uint32_t maxMaterialId = 0;
        for (uint32_t i = 0; i < numMaterials; ++i)
        {
            MaterialView *m = &table->materials[i];
            m->materialId = cursor_u32(&cursor);
            m->nRho       = cursor_u32(&cursor);
            m->nU         = cursor_u32(&cursor);
            uint32_t flags = cursor_u32(&cursor);
            if (cursor.error) goto truncated;

            m->rowDependentU = (version >= 2u && (flags & FLAG_ROW_DEPENDENT_U)) ? 1 : 0;
            m->nativeGrid    = (version >= 2u && (flags & FLAG_NATIVE_TGRID)) ? 1 : 0;
            m->hasEntropy    = (version >= 3u && (flags & FLAG_HAS_ENTROPY)) ? 1 : 0;
            uint32_t supportedFlags = (version >= 3u ? 15u : 7u);
            if ((flags & ~supportedFlags) != 0u || (flags & 1u) == 0u)
            {
                fprintf(stderr, "eos_table: %s material %u uses an unsupported value encoding\n",
                        path, m->materialId);
                goto fail;
            }
            if (m->nRho < 2 || m->nU < 2)
            {
                fprintf(stderr, "eos_table: %s material %u grid dimensions are invalid\n",
                        path, m->materialId);
                goto fail;
            }
            uint64_t valueCount64     = (uint64_t)m->nRho * m->nU;
            uint64_t energyAxisCount  = m->rowDependentU ? valueCount64 : (uint64_t)m->nU;
            if (valueCount64 > SIZE_MAX / sizeof(float))
            {
                fprintf(stderr, "eos_table: %s material %u grid is too large\n", path, m->materialId);
                goto fail;
            }

            m->logRho     = cursor_doubles(&cursor, m->nRho);
            m->energyAxis = cursor_doubles(&cursor, (size_t)energyAxisCount);
            size_t valueCount = (size_t)valueCount64;
            m->logPressure    = cursor_floats(&cursor, valueCount);
            m->logSoundSpeed  = cursor_floats(&cursor, valueCount);
            m->logTemperature = cursor_floats(&cursor, valueCount);
            if (m->hasEntropy) m->entropy = cursor_floats(&cursor, valueCount);
            if (cursor.error) goto truncated;

            if (!axis_is_strictly_increasing(m->logRho, m->nRho))
            {
                fprintf(stderr, "eos_table: %s material %u density axis is not increasing\n",
                        path, m->materialId);
                goto fail;
            }
            if (!m->rowDependentU && !axis_is_strictly_increasing(m->energyAxis, m->nU))
            {
                fprintf(stderr, "eos_table: %s material %u energy axis is not increasing\n",
                        path, m->materialId);
                goto fail;
            }
            m->uStride = m->rowDependentU ? m->nU : 0u;

            /* Fallback sound speed / temperature are minima of the decoded
             * (linear) values; stored values are logs, so take exp of the min. */
            double minCs = DBL_MAX, minT = DBL_MAX;
            for (size_t k = 0; k < valueCount; ++k)
            {
                if ((double)m->logSoundSpeed[k] < minCs) minCs = m->logSoundSpeed[k];
                if ((double)m->logTemperature[k] < minT) minT = m->logTemperature[k];
            }
            m->minSoundSpeed  = exp(minCs);
            m->minTemperature = exp(minT);

            if (m->materialId > maxMaterialId) maxMaterialId = m->materialId;
        }
        if (cursor.pos != cursor.size)
        {
            fprintf(stderr, "eos_table: %s contains trailing payload data\n", path);
            goto fail;
        }

        table->indexSize     = maxMaterialId + 1;
        table->materialIndex = (int *)malloc((size_t)table->indexSize * sizeof(int));
        if (!table->materialIndex) goto fail;
        for (uint32_t i = 0; i < table->indexSize; ++i) table->materialIndex[i] = -1;
        for (uint32_t i = 0; i < numMaterials; ++i)
        {
            uint32_t id = table->materials[i].materialId;
            if (table->materialIndex[id] >= 0)
            {
                fprintf(stderr, "eos_table: %s has duplicate material ID %u\n", path, id);
                goto fail;
            }
            table->materialIndex[id] = (int)i;
        }
    }
    return table;

truncated:
    fprintf(stderr, "eos_table: %s material payload is truncated\n", path);
fail:
    eos_table_free_members(table);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* metadata queries                                                    */
/* ------------------------------------------------------------------ */

uint32_t eos_table_num_materials(const EosTable *table)
{
    return table ? table->numMaterials : 0;
}

uint32_t eos_table_material_id(const EosTable *table, uint32_t index)
{
    if (!table || index >= table->numMaterials) return 0;
    return table->materials[index].materialId;
}

static const MaterialView *find_material(const EosTable *table, uint32_t materialId)
{
    if (!table || materialId >= table->indexSize) return NULL;
    int index = table->materialIndex[materialId];
    if (index < 0) return NULL;
    return &table->materials[index];
}

int eos_table_is_native(const EosTable *table, uint32_t materialId)
{
    const MaterialView *m = find_material(table, materialId);
    return m ? m->nativeGrid : -1;
}

int eos_table_has_entropy(const EosTable *table, uint32_t materialId)
{
    const MaterialView *m = find_material(table, materialId);
    return m ? m->hasEntropy : -1;
}

int eos_table_bounds(const EosTable *table, uint32_t materialId,
                     double *rhoMin, double *rhoMax, double *uMin, double *uMax)
{
    const MaterialView *m = find_material(table, materialId);
    if (!m) return -1;
    if (rhoMin) *rhoMin = exp(m->logRho[0]);
    if (rhoMax) *rhoMax = exp(m->logRho[m->nRho - 1]);
    if (uMin || uMax)
    {
        uint64_t count = m->rowDependentU ? (uint64_t)m->nRho * m->nU : (uint64_t)m->nU;
        double lo = m->energyAxis[0], hi = m->energyAxis[0];
        for (uint64_t i = 1; i < count; ++i)
        {
            if (m->energyAxis[i] < lo) lo = m->energyAxis[i];
            if (m->energyAxis[i] > hi) hi = m->energyAxis[i];
        }
        if (!m->nativeGrid) { lo = exp(lo); hi = exp(hi); }
        if (uMin) *uMin = lo;
        if (uMax) *uMax = hi;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* interpolation (ported from tabulated_eos.cpp, keep logic identical) */
/* ------------------------------------------------------------------ */

static size_t lower_cell(const double *axis, size_t n, double value)
{
    /* upper_bound(axis, value) - 1, clamped to [0, n-2] */
    size_t lo = 0, hi = n; /* first index with axis[i] > value in [lo, hi) */
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (value < axis[mid]) hi = mid;
        else                   lo = mid + 1;
    }
    if (lo == 0) return 0;
    if (lo == n) return n - 2;
    return lo - 1;
}

static double bilinear(const float *values, unsigned nU, size_t ir, size_t iu, double fr, double fu)
{
    const float *row0 = values + ir * (size_t)nU;
    const float *row1 = values + (ir + 1) * (size_t)nU;
    double lo = lerpd((double)row0[iu], (double)row0[iu + 1], fu);
    double hi = lerpd((double)row1[iu], (double)row1[iu + 1], fu);
    return lerpd(lo, hi, fr);
}

/* Decode (exp) the four nodes first, then interpolate linearly. Used for
 * native grids, where the u coordinate is linear and a log-space
 * interpolation can produce negative values across the huge dynamic range. */
static double bilinear_decoded(const float *values, unsigned nU, size_t ir, size_t iu, double fr, double fu)
{
    const float *row0 = values + ir * (size_t)nU;
    const float *row1 = values + (ir + 1) * (size_t)nU;
    double lo = lerpd(exp((double)row0[iu]), exp((double)row0[iu + 1]), fu);
    double hi = lerpd(exp((double)row1[iu]), exp((double)row1[iu + 1]), fu);
    return lerpd(lo, hi, fr);
}

/* Virtual energy value at fractional density row (ir, fr) and energy index
 * iu. With a shared axis (stride 0) the row indices are irrelevant. */
static double virtual_energy(const double *axis, unsigned stride, size_t ir, size_t iu, double fr)
{
    if (stride == 0) return axis[iu];
    double lo = axis[ir * stride + iu];
    double hi = axis[(ir + 1) * stride + iu];
    return lerpd(lo, hi, fr);
}

/* Same, but the axis is always row-dependent (native grids). */
static double virtual_energy_linear(const double *axis, unsigned stride, size_t ir, size_t iu, double fr)
{
    double lo = axis[ir * stride + iu];
    double hi = axis[(ir + 1) * stride + iu];
    return lerpd(lo, hi, fr);
}

static void native_energy_bounds(const double *axis, unsigned stride, unsigned nU, size_t ir, double fr,
                                 double *minOut, double *maxOut)
{
    double minimum = virtual_energy_linear(axis, stride, ir, 0, fr);
    double maximum = minimum;
    for (size_t iu = 1; iu < nU; ++iu)
    {
        double value = virtual_energy_linear(axis, stride, ir, iu, fr);
        if (value < minimum) minimum = value;
        if (value > maximum) maximum = value;
    }
    *minOut = minimum;
    *maxOut = maximum;
}

static size_t lower_virtual_energy_cell(const double *axis, unsigned stride, unsigned nU,
                                        size_t ir, double fr, double value)
{
    size_t lo = 0, hi = nU - 1;
    while (hi - lo > 1)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (virtual_energy(axis, stride, ir, mid, fr) <= value) lo = mid;
        else hi = mid;
    }
    return lo;
}

static int usable_native_energy_segment(double z0, double z1, double value)
{
    const double relativeTolerance = 64.0 * DBL_EPSILON;
    double scale = 1.0;
    double a0 = fabs(z0), a1 = fabs(z1), av = fabs(value);
    if (a0 > scale) scale = a0;
    if (a1 > scale) scale = a1;
    if (av > scale) scale = av;
    return isfinite(z0) && isfinite(z1) && fabs(z1 - z0) > relativeTolerance * scale;
}

static int contains_native_energy(double z0, double z1, double value)
{
    double zLow  = lerpd(z0, z1, -0.0001);
    double zHigh = lerpd(z0, z1, 1.0001);
    return (zLow <= zHigh && value >= zLow && value <= zHigh) ||
           (zHigh < zLow && value >= zHigh && value <= zLow);
}

typedef struct
{
    size_t cell;
    double fraction;
} NativeEnergyLocation;

static NativeEnergyLocation native_energy_location(const double *axis, unsigned stride, unsigned nU,
                                                   size_t ir, double fr, double value)
{
    size_t a = 0, b = nU - 2, c = (a + b) / 2;
    while (1)
    {
        double z0 = virtual_energy_linear(axis, stride, ir, c, fr);
        double z1 = virtual_energy_linear(axis, stride, ir, c + 1, fr);
        if (usable_native_energy_segment(z0, z1, value) && contains_native_energy(z0, z1, value))
        {
            NativeEnergyLocation loc = {c, clampd((value - z0) / (z1 - z0), 0.0, 1.0)};
            return loc;
        }
        if (b - a < 2) break;
        double zLow  = z0 < z1 ? z0 : z1;
        double zHigh = z0 < z1 ? z1 : z0;
        if (value > zHigh) a = c;
        else if (value < zLow) b = c;
        else break;
        c = (a + b) / 2;
    }
    for (size_t j = 0; j + 1 < nU; ++j)
    {
        double z0 = virtual_energy_linear(axis, stride, ir, j, fr);
        double z1 = virtual_energy_linear(axis, stride, ir, j + 1, fr);
        if (usable_native_energy_segment(z0, z1, value) && contains_native_energy(z0, z1, value))
        {
            NativeEnergyLocation loc = {j, clampd((value - z0) / (z1 - z0), 0.0, 1.0)};
            return loc;
        }
    }

    /* A native ANEOS temperature column is not guaranteed to be monotonic in
     * energy and may contain plateaus. If no segment brackets the requested
     * energy, use the closest endpoint of a non-degenerate segment. This is
     * deterministic and always produces a finite thermodynamic state. */
    NativeEnergyLocation closest = {0, 0.0};
    double closestDistance = HUGE_VAL;
    for (size_t j = 0; j + 1 < nU; ++j)
    {
        double z0 = virtual_energy_linear(axis, stride, ir, j, fr);
        double z1 = virtual_energy_linear(axis, stride, ir, j + 1, fr);
        if (!usable_native_energy_segment(z0, z1, value)) continue;
        double d0 = fabs(value - z0);
        double d1 = fabs(value - z1);
        double distance = d0 < d1 ? d0 : d1;
        if (distance < closestDistance)
        {
            closestDistance = distance;
            closest.cell    = j;
            closest.fraction = d0 <= d1 ? 0.0 : 1.0;
        }
    }
    if (isfinite(closestDistance)) return closest;

    /* Entirely flat columns carry no invertible T(rho,u) information. Select
     * the closest stored node without interpolation. */
    size_t closestNode = 0;
    closestDistance    = HUGE_VAL;
    for (size_t j = 0; j < nU; ++j)
    {
        double distance = fabs(value - virtual_energy_linear(axis, stride, ir, j, fr));
        if (distance < closestDistance)
        {
            closestDistance = distance;
            closestNode     = j;
        }
    }
    NativeEnergyLocation loc;
    if (closestNode + 1 < nU) { loc.cell = closestNode;     loc.fraction = 0.0; }
    else                      { loc.cell = closestNode - 1; loc.fraction = 1.0; }
    return loc;
}

/* ------------------------------------------------------------------ */
/* evaluation                                                          */
/* ------------------------------------------------------------------ */

EosTableState eos_table_evaluate(const EosTable *table, double rho, double u,
                                 uint32_t materialId, int requestEntropy)
{
    EosTableState state = {0.0, 0.0, 0.0, 0.0, EOS_TABLE_UNKNOWN_MATERIAL, 0};
    const MaterialView *m = find_material(table, materialId);
    if (!m)
    {
        state.status = EOS_TABLE_UNKNOWN_MATERIAL;
        return state;
    }

    int hasEntropy = requestEntropy && m->hasEntropy;
    state.hasEntropy = hasEntropy;

    if (!(rho > 0.0) || !isfinite(rho) || !isfinite(u))
    {
        state.pressure    = 0.0;
        state.soundSpeed  = m->minSoundSpeed;
        state.temperature = m->minTemperature;
        state.status      = EOS_TABLE_INVALID_INPUT;
        return state;
    }
    if (!m->nativeGrid && !(u > 0.0))
    {
        state.pressure    = 0.0;
        state.soundSpeed  = m->minSoundSpeed;
        state.temperature = m->minTemperature;
        state.status      = EOS_TABLE_INVALID_INPUT;
        return state;
    }

    const double *rhoAxis = m->logRho;
    const double *uAxis   = m->energyAxis;
    unsigned nU           = m->nU;
    unsigned uStride      = m->uStride;

    double logRho = log(rho);
    double logU   = m->nativeGrid ? 0.0 : log(u);
    int status    = EOS_TABLE_SUCCESS;
    if (logRho < rhoAxis[0]) status = EOS_TABLE_DENSITY_BELOW_RANGE;
    else if (logRho > rhoAxis[m->nRho - 1]) status = EOS_TABLE_DENSITY_ABOVE_RANGE;
    logRho = clampd(logRho, rhoAxis[0], rhoAxis[m->nRho - 1]);

    size_t ir = lower_cell(rhoAxis, m->nRho, logRho);
    double fr = (logRho - rhoAxis[ir]) / (rhoAxis[ir + 1] - rhoAxis[ir]);
    if (m->nativeGrid)
    {
        /* native grids interpolate the axes linearly in rho, not in log rho.
           Use the CLAMPED density (exp of the clamped logRho) so out-of-range
           states sit exactly on the boundary row -- using the raw rho would
           silently extrapolate (fr < 0 or fr > 1), giving negative pressures
           below range and runaway pressures above range. Matches sphexa's
           tabulated_eos.cpp (rhoInterpolated) and tabulated_eos_gpu.cu. */
        double rho0 = exp(rhoAxis[ir]), rho1 = exp(rhoAxis[ir + 1]);
        double rhoInterpolated = exp(logRho);
        fr = (rhoInterpolated - rho0) / (rho1 - rho0);
    }

    double energyCoordinate = m->nativeGrid ? u : logU;
    if (m->nativeGrid)
    {
        double nativeMin, nativeMax;
        native_energy_bounds(uAxis, uStride, nU, ir, fr, &nativeMin, &nativeMax);
        if (status == EOS_TABLE_SUCCESS && energyCoordinate < nativeMin)
            status = EOS_TABLE_ENERGY_BELOW_RANGE;
        else if (status == EOS_TABLE_SUCCESS && energyCoordinate > nativeMax)
            status = EOS_TABLE_ENERGY_ABOVE_RANGE;
        energyCoordinate = clampd(energyCoordinate, nativeMin, nativeMax);
    }
    else
    {
        double uMin = virtual_energy(uAxis, uStride, ir, 0, fr);
        double uMax = virtual_energy(uAxis, uStride, ir, nU - 1, fr);
        if (status == EOS_TABLE_SUCCESS && energyCoordinate < uMin)
            status = EOS_TABLE_ENERGY_BELOW_RANGE;
        else if (status == EOS_TABLE_SUCCESS && energyCoordinate > uMax)
            status = EOS_TABLE_ENERGY_ABOVE_RANGE;
        energyCoordinate = clampd(energyCoordinate, uMin, uMax);
        logU = energyCoordinate;
    }

    NativeEnergyLocation nativeLocation = {0, 0.0};
    size_t iu;
    if (m->nativeGrid)
    {
        nativeLocation = native_energy_location(uAxis, uStride, nU, ir, fr, energyCoordinate);
        iu = nativeLocation.cell;
    }
    else
    {
        iu = lower_virtual_energy_cell(uAxis, uStride, nU, ir, fr, logU);
    }

    double fu;
    if (m->nativeGrid)
    {
        fu = nativeLocation.fraction;
    }
    else
    {
        double u0 = virtual_energy(uAxis, uStride, ir, iu, fr);
        double u1 = virtual_energy(uAxis, uStride, ir, iu + 1, fr);
        fu = (logU - u0) / (u1 - u0);
    }
    if (!isfinite(fu))
    {
        state.pressure    = 0.0;
        state.soundSpeed  = m->minSoundSpeed;
        state.temperature = m->minTemperature;
        state.hasEntropy  = 0;
        state.status      = EOS_TABLE_INVALID_TABLE_STATE;
        return state;
    }
    /* The segment-matching tolerance must not turn into thermodynamic
     * extrapolation: even a 1e-4 overshoot can make a linear interpolation
     * of positive decoded values negative across the ANEOS dynamic range. */
    fu = clampd(fu, 0.0, 1.0);

    double entropy = hasEntropy ? bilinear(m->entropy, nU, ir, iu, fr, fu) : 0.0;
    if (m->nativeGrid)
    {
        state.pressure    = bilinear_decoded(m->logPressure, nU, ir, iu, fr, fu);
        state.soundSpeed  = bilinear_decoded(m->logSoundSpeed, nU, ir, iu, fr, fu);
        state.temperature = bilinear_decoded(m->logTemperature, nU, ir, iu, fr, fu);
    }
    else
    {
        state.pressure    = exp(bilinear(m->logPressure, nU, ir, iu, fr, fu));
        state.soundSpeed  = exp(bilinear(m->logSoundSpeed, nU, ir, iu, fr, fu));
        state.temperature = exp(bilinear(m->logTemperature, nU, ir, iu, fr, fu));
    }
    state.entropy = entropy;
    state.status  = status;

    if (!(state.pressure > 0.0) || !isfinite(state.pressure) ||
        !(state.soundSpeed > 0.0) || !isfinite(state.soundSpeed) ||
        !(state.temperature > 0.0) || !isfinite(state.temperature))
    {
        state.pressure    = 0.0;
        state.soundSpeed  = m->minSoundSpeed;
        state.temperature = m->minTemperature;
        state.entropy     = 0.0;
        state.status      = EOS_TABLE_INVALID_TABLE_STATE;
        return state;
    }
    return state;
}

/* ------------------------------------------------------------------ */
/* unit conversion                                                     */
/* ------------------------------------------------------------------ */

int eos_table_units_cgs(EosTableUnits *units, double densityUnitCgs, double energyUnitCgs)
{
    if (!units) return -1;
    if (!(densityUnitCgs > 0.0) || !(energyUnitCgs > 0.0) ||
        !isfinite(densityUnitCgs) || !isfinite(energyUnitCgs))
    {
        return -1;
    }
    double pressureUnitCgs = densityUnitCgs * energyUnitCgs;
    double speedUnitCgs    = sqrt(energyUnitCgs);
    if (!isfinite(pressureUnitCgs) || !isfinite(speedUnitCgs) ||
        pressureUnitCgs <= 0.0 || speedUnitCgs <= 0.0)
    {
        return -1;
    }
    units->densityToTable      = densityUnitCgs;
    units->energyToTable       = energyUnitCgs;
    units->pressureFromTable   = 1.0 / pressureUnitCgs;
    units->soundSpeedFromTable = 1.0 / speedUnitCgs;
    units->entropyFromTable    = 1.0 / energyUnitCgs;
    return 0;
}

EosTableState eos_table_evaluate_code(const EosTable *table, double rho, double u,
                                      uint32_t materialId, int requestEntropy,
                                      EosTableUnits units)
{
    EosTableState state = eos_table_evaluate(table, rho * units.densityToTable,
                                             u * units.energyToTable, materialId, requestEntropy);
    state.pressure   *= units.pressureFromTable;
    state.soundSpeed *= units.soundSpeedFromTable;
    state.entropy    *= units.entropyFromTable;
    return state;
}

const char *eos_table_status_string(int status)
{
    switch (status)
    {
    case EOS_TABLE_SUCCESS:             return "success";
    case EOS_TABLE_DENSITY_BELOW_RANGE: return "densityBelowRange";
    case EOS_TABLE_DENSITY_ABOVE_RANGE: return "densityAboveRange";
    case EOS_TABLE_ENERGY_BELOW_RANGE:  return "energyBelowRange";
    case EOS_TABLE_ENERGY_ABOVE_RANGE:  return "energyAboveRange";
    case EOS_TABLE_UNKNOWN_MATERIAL:    return "unknownMaterial";
    case EOS_TABLE_INVALID_INPUT:       return "invalidInput";
    case EOS_TABLE_ENTROPY_UNAVAILABLE: return "entropyUnavailable";
    case EOS_TABLE_INVALID_TABLE_STATE: return "invalidTableState";
    default:                            return "unknownStatus";
    }
}
