/*
 * The header file for testlibmaneos.c.
 */
#ifndef TESTLIBMANEOS_HINCLUDED
#define TESTLIBMANEOS_HINCLUDED

/*
 * Each entry has (T, P, S, cv, cs) as a function of rho and u.
 */
typedef struct ANEOSLOOKUPENTRY
{
	double T;
	double P;
	double S;		// specific entropy
	double cv;      // specific heat capacity
	double cs;
    double rho;
    double u;
} ANEOSLookupEntry;

typedef struct ANEOSTABLE
{
	int iMat;
	int nTableRho;
	int nTableU;
    double rho_min,rho_max,drho;
    double u_min, u_max, du;
	ANEOSLookupEntry **Lookup;
} ANEOSTable;

#if 0
/*
 * The final version will only have (T, P, S, cv, cs)
 */
struct ANOESLookupEntry
{
	double T;
	double rho;
	int iMat;
	double P;
	double u;
	double S;		// specific entropy
	double cv;
	double dPdT;
	double dPdrho;
	double fkros;
	double cs;
	int iPhase;
	double rhoL;
	double rhoH;
	double ion;
}
#endif

ANEOSTable *ANEOSInitTable(int iMat, double rho_min, double rho_max, double drho, double u_min, double u_max, double du);
void ANEOSFinalizeTable(ANEOSTable *table);
int ANEOSTableAlloc(ANEOSTable *table, int nRho, int nU);
ANEOSLookupEntry **ANEOSMatrixAlloc(int nRow, int nCol);
int ANEOSTableRead(ANEOSTable *table, char *file, int nRho, int nU);
void ANEOSTablePrint(ANEOSTable *table, char *file);
int ANEOSLookupRhoIndex(ANEOSLookupEntry **Lookup, double rho, unsigned int nRow, unsigned int iCol);
void ANEOSTableComparetoEOS(ANEOSTable *table);
#endif
