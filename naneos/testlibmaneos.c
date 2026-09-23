/*
 ** This is a simple program to test the C wrapper to libmaneos.f.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__) && !defined(__MACH__)
#include <malloc.h>
#endif
#include <assert.h>
#include "maneos.h"
#include "testlibmaneos.h"

/*
 * Initialize the lookup table.
 */
ANEOSTable *ANEOSInitTable(int iMat, double rho_min, double rho_max, double drho, double u_min, double u_max, double du)
{
	ANEOSTable *table;
	int i;

	table = (ANEOSTable*) malloc(sizeof(ANEOSTable));
	assert(table != NULL);
	
	table->iMat = iMat;
	table->nTableRho = 0;
	table->nTableU = 0;

    table->rho_min = rho_min;
    table->rho_max = rho_max;
    table->drho = drho;

    table->u_min = u_min;
    table->u_max = u_max;
    table->du = du;

#if 0
	/*
	 * We allocate the array as described in the Numerical Recipes.
	 */	
	table->Lookup = (ANEOSLookupEntry**) malloc(table->nTableRho*sizeof(ANEOSLookupEntry*));
	table->Lookup[0] = (ANEOSLookupEntry*) malloc(table->nTableRho*table->nTableU*sizeof(ANEOSLookupEntry));

	for (i=1; i<table->nTableRho; i++)
	{
		table->Lookup[i] = table->Lookup[i-1]+table->nTableU;
	}
#endif
	return (table);
}

/*
 * Free memory.
 */
void ANEOSFinalizeTable(ANEOSTable *table)
{
    int i;
	if (table != NULL)
	{
        if (table->Lookup != NULL)
        {
            /* Set a pointer to each row. */
            for (i=1; i<table->nTableRho; i++)
            {
                free(table->Lookup[i]);
            }
            
            free (table->Lookup);
        }
        free(table);
	}
}

/*
 * Allocate memory for the look up table.
 */
int ANEOSTableAlloc(ANEOSTable *table, int nRho, int nU)
{
	int i;

	assert(table != NULL);

	/*
	 * nRho: number of rows
	 * nU:   number of columns
	 */
	table->nTableRho = nRho;
	table->nTableU = nU;
	
	table->Lookup = ANEOSMatrixAlloc(nRho, nU);
    assert(table->Lookup != NULL);
}

/*
 * We allocate a lookup table as an array as described in the Numerical Recipes.
 */	
ANEOSLookupEntry **ANEOSMatrixAlloc(int nRow, int nCol)
{
    ANEOSLookupEntry **Lookup;
	int i;

	Lookup = (ANEOSLookupEntry**) malloc(nRow*sizeof(ANEOSLookupEntry*));
	Lookup[0] = (ANEOSLookupEntry*) malloc(nRow*nCol*sizeof(ANEOSLookupEntry));

	assert(Lookup != NULL);

	/* Set a pointer to each row. */
	for (i=1; i<nRow; i++)
	{
		Lookup[i] = Lookup[i-1]+nCol;
	}

	return (Lookup);
}

/*
 * Read the lookup table from an ASCII file. The density and internal energy of each entry are obtained from rho_min and u_min
 * assuming a logarithmic spacing.
 */
int ANEOSTableRead(ANEOSTable *table, char *file, int nRho, int nU)
{
	FILE *fp;
	int iRet, i, j;

	assert(table != NULL);
	assert(table->Lookup == NULL);

    /* Allocate memory. */
    ANEOSTableAlloc(table, nRho, nU);

	/* Open the file. */
	fp = fopen(file,"r");
	assert(fp != NULL);

	/* Read the table. */
	for (i=0; i<table->nTableRho; i++)
	{
		for (j=0; j<table->nTableU; j++)
		{
			// Read one tuple (T, P, S, cv, cs)
			iRet = fscanf(fp, "%lf %lf %lf %lf %lf", &table->Lookup[i][j].T,&table->Lookup[i][j].P,&table->Lookup[i][j].S,&table->Lookup[i][j].cv,&table->Lookup[i][j].cs);
			if (iRet <= 0 && feof(fp))
			{
				// Error occured
				fprintf(stderr,"i=%i j=%i iRet=%i\n",i,j,iRet);
				iRet = 1;
				return (iRet);
//				break;
			}
			assert(iRet > 0);
            /*
             * Check that the values are sensible.
             */
#if 0
            assert(table->Lookup[i][j].T >= 0.0);
            assert(table->Lookup[i][j].P >= 0.0);
            assert(table->Lookup[i][j].S >= 0.0);
            assert(table->Lookup[i][j].cv >= 0.0);
            assert(table->Lookup[i][j].cs >= 0.0);
#endif
            // Logarithmic spacing
            table->Lookup[i][j].rho = table->rho_min * pow(10.0, i*table->drho);
            table->Lookup[i][j].u = table->u_min * pow(10.0, j*table->du);
        }
	}
	
	// File read without problems	
	iRet = 0;
	return (iRet);
}

/*
 * Print the lookup table to an ASCII file.
 */
void ANEOSTablePrint(ANEOSTable *table, char *file)
{
	FILE *fp;
	int iRet, i, j;

	assert(table != NULL);
	assert(table->Lookup != NULL);

	/* Open the file. */
	fp = fopen(file,"w");
	assert(fp != NULL);

	/* Print the table. */
	for (i=0; i<table->nTableRho; i++)
	{
		for (j=0; j<table->nTableU; j++)
		{
            // Print one tuple (T, P, S, cv, cs)
			fprintf(fp, "%24.15E%24.15E%24.15E%24.15E%24.15E", table->Lookup[i][j].T,table->Lookup[i][j].P,table->Lookup[i][j].S,table->Lookup[i][j].cv,table->Lookup[i][j].cs);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

/*
 * Do bisection to find the values rho_i and rho_i+1 that
 * bracket rho in one column of the the lookup table assuming
 * it is sorted. The function returns the index i.
 */
int ANEOSLookupRhoIndex(ANEOSLookupEntry **Lookup, double rho, unsigned int nRow, unsigned int iCol)
{
	unsigned int iLower,iUpper,i;
	
	iLower = 0;
	iUpper = nRow-1;

	/*
	 * Make sure that rho is in the lookup table.
	 * If this is not the case return -1 or
	 * nRow. */
	if (rho < Lookup[iLower][iCol].rho)
	{
		return (-1);
	} else if (rho > Lookup[iUpper][iCol].rho)
	{
		return (nRow);
	}

	assert(Lookup[iLower][iCol].rho < Lookup[iUpper][iCol].rho);
	
	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (rho >= Lookup[i][iCol].rho)
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (rho == Lookup[0][iCol].rho) return 0;
	/* 
	 * Return rhomax-1 as the lower index, so the
	 * desired value is always bracketed by (i,i+1).
	 */
	if (rho == Lookup[nRow-1][iCol].rho)
	{
			return (nRow-2);
	}

	printf("rho= %g rho(imax,j)= %g iLower= %i\n",rho, Lookup[nRow-1][iCol].rho,iLower);
	return iLower;
}

/*
 * Do bisection to find the values u_j and u_j+1 that
 * bracket u in one row of the the lookup table assuming
 * it is sorted. The function returns the index i.
 */
int ANEOSLookupUIndex(ANEOSLookupEntry **Lookup, double u, unsigned int iRow, unsigned int nCol)
{
	unsigned int iLower,iUpper,i;
	
	iLower = 0;
	iUpper = nCol-1;

	/*
	 * Make sure that T is in the lookup table.
	 * If this is not the case return -1 or
	 * nCol. */
	if (u < Lookup[iRow][iLower].u)
	{
		return (-1);
	} else if (u > Lookup[iRow][iUpper].u)
    {
		return (nCol);
	}

	assert(Lookup[iRow][iLower].u < Lookup[iRow][iUpper].u);
	
	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (u >= Lookup[iRow][i].u)
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (u == Lookup[iRow][0].u) return 0;
	/* 
	 * Return umax-1 as the lower index, so the
	 * desired value is always bracketed by (i,i+1).
	 */
	if (u == Lookup[iRow][nCol-1].u)
	{
			return (nCol-2);
	}

	printf("u= %g u(i,jmax)= %g iLower= %i\n",u, Lookup[iRow][nCol-1].u,iLower);
	return iLower;
}

void ANEOSInterpolateRhoU(ANEOSTable *table, double rho, double u, double *pP, double *pcs)
{
    double x,y;
    int i,j;

    /*
     * Since rho(i,j) = rho(i) and u(i,j) = u(j) we can do two 1D bisection to determine i and j.
     */
    i = ANEOSLookupRhoIndex(table->Lookup, rho, table->nTableRho, 0);
    printf("i=%i ",i);
    
    j = ANEOSLookupUIndex(table->Lookup, u, 0, table->nTableU);
    printf("j=%i ",j);

    /*
     * Check if the data is in the lookup table.
     */
    if (i < 0 || i >= table->nTableRho)
    {
        printf("rho= %15.7E i= %i is outside of the table.\n", rho, i);
        assert(0);
    }

    if (j < 0 || j >= table->nTableU)
    {
        printf("u= %15.7E j= %i is outside of the table.\n", u, j);
        assert(0);
    }

    /*
     * We assume, that the grid is evenly spaced in rho and u.
     */
    x = (table->Lookup[i+1][j].rho - rho)/(table->Lookup[i+1][j].rho - table->Lookup[i][j].rho);
    y = (table->Lookup[i][j+1].u - u)/(table->Lookup[i][j+1].u-table->Lookup[i][j].u);

    *pP = x*y*table->Lookup[i][j].P + x*(1.0-y)*table->Lookup[i][j+1].P +
        (1.0-x)*y*table->Lookup[i+1][j].P + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].P;

    *pcs = x*y*table->Lookup[i][j].cs + x*(1.0-y)*table->Lookup[i][j+1].cs +
        (1.0-x)*y*table->Lookup[i+1][j].cs + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].cs;

    /*
     * Some code for debugging the interpolator.
     */
    double P = *pP;
    fprintf(stderr,"rho= %g u= %g i=%i, j=%i\n",rho,u,i,j);
    fprintf(stderr,"P[%i][%i]   = %15.7E\n",i,j,table->Lookup[i][j].P);
    fprintf(stderr,"P[%i][%i]   = %15.7E\n",i+1,j,table->Lookup[i+1][j].P);
    fprintf(stderr,"P[%i][%i]   = %15.7E\n",i,j+1,table->Lookup[i][j+1].P);
    fprintf(stderr,"P[%i][%i]   = %15.7E\n",i+1,j+1,table->Lookup[i+1][j+1].P);
    fprintf(stderr,"P           = %15.7E\n",P);
}

/*
 * Check the lookup table by comparing it to callaneos().
 */
void ANEOSTableComparetoEOS(ANEOSTable *table)
{
	double T, rho, u, S;
	double P, cv, dPdT, dPdrho, fkros, cs, rhoL, rhoH, ion;
	int iMat, iPhase;
    int i,j;
    
    assert (table->Lookup != NULL);

    iMat = 2;

    fprintf(stderr, "Comparing the lookup table to the MANEOS routine.\n");
    for (i=600; i<800; i++)
    {
        for (j=600; j<800; j++)
        {
            rho = table->Lookup[i][j].rho;
            T = table->Lookup[i][j].T;
            
//            fprintf(stderr,"i= %i j= %i\n", i, j);
            // Call the ANEOS routine to determine P, u, S, cv and cs.
	        callaneos(T, rho, iMat, &P, &u, &S, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH, &ion);
            if (fabs(table->Lookup[i][j].P-P) > 1e-30)
            {
                fprintf(stderr,"P= %15.7E %15.7E\n",table->Lookup[i][j].P,P);
                assert(fabs(table->Lookup[i][j].P-P) < 1e-30);
            }

            if (fabs(table->Lookup[i][j].u-u) < 1e-30)
            {
                assert(fabs(table->Lookup[i][j].u-u) < 1e-30);
            }
            
            if (fabs(table->Lookup[i][j].S-S) < 1e-30)
            {
                assert(fabs(table->Lookup[i][j].S-S) < 1e-30);
            }
            
            if (fabs(table->Lookup[i][j].cv-cv) < 1e-30)
            {
                assert(fabs(table->Lookup[i][j].cv-cv) < 1e-30);
            }
            
            if (fabs(table->Lookup[i][j].cs-cs) < 1e-30)
            {
                assert(fabs(table->Lookup[i][j].cs-cs) < 1e-30);
            }
        }
    }
}

/*
 * Check if the table is monotonic in rho.
 */
int ANEOSTableRhoIsMonotonic(ANEOSTable *table)
{
    int i,j,iRet;
    
    assert (table->Lookup != NULL);

    // Set iRet to true
    iRet = 1;

    for (j=0; j<table->nTableU; j++)
    {
        for (i=0; i<table->nTableRho-1; i++)
        {
            if (table->Lookup[i][j].rho  <= table->Lookup[i+1][j].rho)
            {
            } else {
                fprintf(stderr,"i= %i j= %i rho[i][j]= %15.7E rho[i+1][j]= %15.7E.\n",i,j,table->Lookup[i][j].rho, table->Lookup[i+1][j].rho);
                iRet = 0;
            }
        }
    }
    return (iRet);
}
/*
 * Check if the table is monotonic in u.
 */
int ANEOSTableUIsMonotonic(ANEOSTable *table)
{
    int i,j,iRet;
    
    assert (table->Lookup != NULL);

    // Set iRet to true
    iRet = 1;

    for (i=0; i<table->nTableRho; i++)
    {
        for (j=0; j<table->nTableU-1; j++)
        {
            if (table->Lookup[i][j].u  <= table->Lookup[i][j+1].u)
            {
            } else {
                fprintf(stderr,"i= %i j= %i u[i][j]= %15.7E u[i][j+1]= %15.7E.\n",i,j,table->Lookup[i][j].u, table->Lookup[i][j+1].u);
                iRet = 0;
            }
        }
    }
    return (iRet);
}

int main(int argc, char *argv[])
{
#if 0
	double T, rho, u, S;
	double P, cv, dPdT, dPdrho, fkros, cs, rhoL, rhoH, ion;
	int iMat, iPhase;
#endif
	ANEOSTable *table;
    // Load Federico's table that was rewritten in terms of rho and u
	char inFile[256] = "aneos_Fe.dat";
    double P, cs;
//    int i,j;
	
    /*
     * The min. and max. values of rho and u in the table (readme-ANEOS.txt). Keep in mind, that the spacing is logarithmic!
     */
    //	table = ANEOSInitTable(1, 1.0e-2, 1.0e10, 0.02, 1.0e6, 1.0e16, 0.02);
    table=ANEOSInitTable(2, 1.0e-6, 1.0e6, 0.01, 1.0e6, 1.0e14, 0.005);
    
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


	fprintf(stderr, "Initializing material...\n");
	initaneos("aneos.input");
    
    ANEOSTableComparetoEOS(table);

    //    ANEOSTableRhoIsMonotonic(table);
    //    ANEOSTableUIsMonotonic(table);


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
