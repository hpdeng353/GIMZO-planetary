/*
 * The header file for aneostable.c.
 */
#ifndef TESTLIBMANEOS_HINCLUDED
#define TESTLIBMANEOS_HINCLUDED

/*
 * Each entry has (u, P, cs, s) as a function of rho and T. rhoarr and tarr are external global variables
 */
typedef struct ANEOSLOOKUPENTRY
{
  double P;
  double S;		// specific entropy
  double cs;
  double u;
} ANEOSLookupEntry;

typedef struct ANEOSTABLE
{
	int iMat;
	int nTableRho;
	int nTableT;
  //    double rho_min,rho_max,drho;
  //    double u_min, u_max, du;
	ANEOSLookupEntry **Lookup;
} ANEOSTable;

// Functions to initialize and finalize the table.
//ANEOSTable *ANEOSInitTable(int iMat, double rho_min, double rho_max, double drho, double u_min, double u_max, double du);
ANEOSTable *ANEOSInitTable(int iMat, int nrho, int ntemp);
//ANEOSTable *ANEOSInitTable(int iMat, double rho_min, double rho_max, double drho, double u_min, double u_max, double du);
void ANEOSFinalizeTable(ANEOSTable *table);
int ANEOSTableAlloc(ANEOSTable *table, int nRho, int nT);
ANEOSLookupEntry **ANEOSMatrixAlloc(int nRow, int nCol);

// Functions to read the lookup table from a file.
int ANEOSTableRead(ANEOSTable *table, char *file, int nRho, int nT);
void ANEOSTableConvertCodeunits(ANEOSTable *table, double Lunit, double Munit, double Tunit);
void ANEOSTablePrint(ANEOSTable *table, char *file);

// Functions to do a lookup in rho and u.
int ANEOSLookupRhoIndex(double *rhoarr, double rho, unsigned int nRow);
int ANEOSLookupUIndex(ANEOSLookupEntry **Lookup, double u, unsigned int iRow, unsigned int nCol);
int ANEOSLookupTIndex(double *Tarr, double T, unsigned int nCol);
void ANEOSInterpolateRhoU(ANEOSTable *table, double *rhoarr, double *tarr, double rho, double u, double *pP, double *pcs, double *pT,double *ps);
void ANEOSInterpolateRhoT(ANEOSTable *table, double *rhoarr, double *tarr, double rho, double T, double *pu, double *pP, double *pcs, double *ps);

// Check the table by comparing it to the values obtained from callaneos().


#endif
