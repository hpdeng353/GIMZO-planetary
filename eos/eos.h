#ifndef EOS_INTERFACE_H
#define EOS_INTERFACE_H

#include "../GIZMO_config.h"


#if (defined(EOS_HELMHOLTZ) || defined(EOS_TILLOTSON) || defined(EOS_ANEOS) || defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)) && !defined(EOS_GENERAL)
#define EOS_GENERAL
#endif

#ifdef EOS_HELMHOLTZ
#define EOS_TABULATED
#define EOS_USES_CGS
#define EOS_CARRIES_YE
#define EOS_CARRIES_ABAR
#define EOS_CARRIES_TEMPERATURE
#define EOS_PROVIDES_ENTROPY
#define EOS_PROVIDES_CV
#endif

#ifdef EOS_ANEOS
/* EOS_ANEOS is the giant-impact master switch: pull in every flag the impact
 * pipeline needs so users only enable EOS_ANEOS in Config.sh.
 *   EOS_TABULATED            - table-based EOS machinery
 *   EOS_CARRIES_TEMPERATURE  - IC/snapshot Temperature block
 *   EOS_CARRIES_ENTROPY      - snapshot Entropy block
 *   MOON                     - giant-impact master (imat field, IO_IMAT, core/mantle handling)
 *   ORELAX                - relaxation (drag, spherical fix, isentropic pin)
 *   READ_IMAT                - read Materials from ICs
 *   CLIPPING                 - domain/merge-split clipping used by impact runs
 *   PREVENT_PARTICLE_MERGE_SPLIT - impact runs never split/merge particles
 */
#define EOS_TABULATED
#define EOS_CARRIES_TEMPERATURE
#define EOS_CARRIES_ENTROPY
#ifndef MOON
#define MOON
#endif
#define ORELAX
#define READ_IMAT
#define CLIPPING
#define PREVENT_PARTICLE_MERGE_SPLIT
#endif

struct eos_input
{
  double rho;         /* Density */
  double eps;         /* Specific internal energy */
#ifdef EOS_CARRIES_YE
  double Ye;          /* Electron fraction */
#endif
#ifdef EOS_CARRIES_ABAR
  double Abar;        /* Mean atomic weight (in atomic mass units) */
#endif
#ifdef EOS_CARRIES_TEMPERATURE
  double temp;        /* Temperature initial guess */
#endif
};

struct eos_output
{
  double press;       /* Pressure */
  double csound;      /* Sound speed */
#ifdef EOS_CARRIES_TEMPERATURE
  double temp;        /* Temperature (in Kelvin) */
#endif
#ifdef EOS_PROVIDES_ENTROPY
  double entropy;     /* Entropy (in CGS) */
#endif
#ifdef EOS_PROVIDES_CV
  double cv;          /* Specific heat at constant volume (in CGS) */
#endif
};


#ifdef EOS_TABULATED
int eos_init(char const * eos_table_fname);
int eos_cleanup();
#endif

int eos_compute(struct eos_input const * in, struct eos_output * out);

#endif

