#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__) && !defined(__MACH__)
#include <malloc.h>
#endif
#include <assert.h>
#include "maneos.h"
int main(int argc, char *argv[])
{
#if 0
	double T, rho, u, S;
	double P, cv, dPdT, dPdrho, fkros, cs, rhoL, rhoH, ion;
	int iMat, iPhase;
#endif
	ANEOSTable *table;
    // Load Federico's table that was rewritten in terms of rho and u
	char inFile[256] = "MANEOS-water/MANEOS-water.dat";
    double P, cs;
//    int i,j;
	
    /*
     * The min. and max. values of rho and u in the table (readme-ANEOS.txt). Keep in mind, that the spacing is logarithmic!
     */
	table = ANEOSInitTable(1, 1.0e-2, 1.0e10, 0.02, 1.0e6, 1.0e16, 0.02);

    fprintf(stderr,"Reading %s.\n",inFile);
	if (ANEOSTableRead(table,inFile,601,501))
	{
			fprintf(stderr,"Failed reading %s.\n", inFile);
			assert(0);
	}

    fprintf(stderr,"Done.\n");
#if 0
    fprintf(stderr,"Writing output to %s (nRho= %i nU= %i).\n",outFile,table->nTableRho,table->nTableU);
	ANEOSTablePrint(table,outFile);
#endif

#if 0
	fprintf(stderr, "Initializing material...\n");
	initaneos();
    
    ANEOSTableComparetoEOS(table);

    ANEOSTableRhoIsMonotonic(table);
    ANEOSTableUIsMonotonic(table);
#endif

    ANEOSInterpolateRhoU(table, table->Lookup[table->nTableRho-1][0].rho, table->Lookup[table->nTableRho-1][table->nTableU-1].u, &P, &cs);



    printf("\n");
    printf("rho_min= %15.7E", table->rho_min);
    printf("  rho_max= %15.7E", table->rho_max);
    printf("  drho= %15.7E", table->drho);
    printf("  nRho= %3d", table->nTableRho);
    printf("\n");

    printf("u_min= %15.7E", table->u_min);
    printf("  u_max= %15.7E", table->u_max);
    printf("  du= %15.7E", table->du);
    printf("  nU= %3d", table->nTableU);
    printf("\n");

    printf("rho[0][0]= %15.7E", table->Lookup[0][0].rho);
    printf("  rho[%i][0]= %15.7E", table->nTableRho, table->Lookup[table->nTableRho-1][0].rho);
    printf("\n");

    printf("u[0][0]  = %15.7E", table->Lookup[0][0].u);
    printf("  u[0][%i]  = %15.7E", table->nTableU, table->Lookup[0][table->nTableU-1].u);
    printf("\n");


    printf("%15.7E", table->u_min*pow(10.0,table->nTableU*table->du));
    printf("%15.7E", table->Lookup[0][table->nTableU-1].u-table->u_max);
    printf("\n");
#if 0
    /*
     * Initialize to some values.
     */
	T = 100.0;
	rho = 3.2;
	iMat = 1;
	P = 0.0;
	u = 0.0;
	S = 0.0;
	cv = 0.0;
	dPdT = 0.0;
	dPdrho = 0.0;
	fkros = 0.0;
	cs = 0.0;
	iPhase = 0;
	rhoL = 0.0;
	rhoH = 0.0;
	ion = 0.0;

    fprintf(stderr, "Calling MANEOS.\n");
	callaneos(T, rho, iMat, &P, &u, &S, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH, &ion);

	fprintf(stderr,"%15.7E", T);
	fprintf(stderr,"%15.7E", rho);
	fprintf(stderr,"%2i", iMat);
	fprintf(stderr,"%15.7E", P);
	fprintf(stderr,"%15.7E", u);
	fprintf(stderr,"%15.7E", S);
	fprintf(stderr,"%15.7E", cv);
	fprintf(stderr,"%15.7E", dPdT);
	fprintf(stderr,"%15.7E", dPdrho);
	fprintf(stderr,"%15.7E", fkros);
	fprintf(stderr,"%15.7E", cs);
	fprintf(stderr,"%2i", iPhase);
	fprintf(stderr,"%15.7E", rhoL);
	fprintf(stderr,"%15.7E", rhoH);
	fprintf(stderr,"%15.7E", ion);
	fprintf(stderr,"\n");
#endif 
}
