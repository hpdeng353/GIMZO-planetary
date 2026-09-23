/*
 ** This file provides a set of routines to handle ANEOS lookup tables.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__) && !defined(__MACH__)
#include <malloc.h>
#endif
#include <assert.h>
#include "maneos.h"
#include "aneostable.h"

/*
 * Initialize the lookup table.
 */
ANEOSTable *ANEOSInitTable(int iMat, int nrho, int ntemp)
{
	ANEOSTable *table;

	table = (ANEOSTable*) malloc(sizeof(ANEOSTable));
	assert(table != NULL);
	
	table->iMat = iMat;
	table->nTableRho = nrho;
	table->nTableT = ntemp;

 
        //        table->Lookup = NULL;
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
            /* Free each row of the matrix. */
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
int ANEOSTableAlloc(ANEOSTable *table, int nRho, int nT)
{
	int i;

	assert(table != NULL);

	/*
	 * nRho: number of rows
	 * nU:   number of columns
	 */
	table->nTableRho = nRho;
	table->nTableT = nT;
	
	table->Lookup = ANEOSMatrixAlloc(nRho, nT);
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
int ANEOSTableRead(ANEOSTable *table, char *file, int nRho, int nT)
{
	FILE *fp;
	int iRet, i, j;

	assert(table != NULL);
        //	assert(table->Lookup == NULL);

        /* Allocate memory. */
        ANEOSTableAlloc(table, nRho, nT);

	/* Open the file. */
	fp = fopen(file,"r");
	assert(fp != NULL);

	/* Read the table. */
	for (j=0; j<table->nTableT; j++)
	  {
	    for (i=0; i<table->nTableRho; i++)
	      {

		iRet = fscanf(fp, "%lf %lf %lf %lf", &table->Lookup[i][j].u,&table->Lookup[i][j].P,&table->Lookup[i][j].cs,&table->Lookup[i][j].S);

		//		if (&table->iMat ==0)"i=%i j=%i u=%lf\n",i,j,&table->Lookup[i][j].u);
		//		fprintf(stderr,"%lf %lf %lf %lf\n", &table->Lookup[i][j].u,&table->Lookup[i][j].P,&table->Lookup[i][j].cs,&table->Lookup[i][j].S);
		if (iRet <= 0 && feof(fp))
		  {
		    // Error occured
		    fprintf(stderr,"i=%i j=%i iRet=%i\n",i,j,iRet);
		    iRet = 1;
		    return (iRet);
		    //				break;
		  }
		assert(iRet > 0);

	      }
	  }
	// File read without problems	
	iRet = 0;
	return (iRet);
}

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
	for (j=0; j<table->nTableT; j++)
	  {
	    for (i=0; i<table->nTableRho; i++)
	      {
		fprintf(fp, "%24.15E%24.15E%24.15E%24.15E\n", table->Lookup[i][j].u,table->Lookup[i][j].P,table->Lookup[i][j].cs,table->Lookup[i][j].S);
	      }
	    
	}
	fclose(fp);
}

/*
 * Convert the table from cgs (and eV for the temperature) to code units.
 */

void ANEOSTableConvertCodeunits(ANEOSTable *table, double Lunit, double Munit, double Tunit)
{
  int i,j;
  double rhounit,vunit,eunit;
  
  rhounit = Munit/Lunit/Lunit/Lunit;
  vunit = Lunit/Tunit;
  eunit = vunit*vunit;
  for (i=0; i<table->nTableRho; i++)
    {
      for (j=0; j<table->nTableT; j++)
      {
        table->Lookup[i][j].P /= (rhounit*eunit);
        table->Lookup[i][j].S /= eunit; 
        table->Lookup[i][j].cs /= vunit;          
        table->Lookup[i][j].u /= eunit;
      }
    }
}

/*
 * Print the lookup table to an ASCII file.
 */


/*
 * Do bisection to find the values rho_i and rho_i+1 that
 * bracket rho in one column of the the lookup table assuming
 * it is sorted. The function returns the index i.
 */
int ANEOSLookupRhoIndex(double *rhoarr, double rho, unsigned int nRow)
{
	unsigned int iLower,iUpper,i;
	
	iLower = 0;
	iUpper = nRow-1;

	/*
	 * Make sure that rho is in the lookup table.
	 * If this is not the case return -1 or
	 * nRow. */
	if (rho < rhoarr[iLower])
	{
		return (-1);
	} else if (rho > rhoarr[iUpper])
	{
		return (nRow);
	}

	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (rho >= rhoarr[i])
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (rho == rhoarr[0]) return 0;
	/* 
	 * Return rhomax-1 as the lower index, so the
	 * desired value is always bracketed by (i,i+1).
	 */
	if (rho == rhoarr[nRow-1])
	{
			return (nRow-2);
	}

	return iLower;
}


int ANEOSLookupTIndex(double *Tarr, double T, unsigned int nCol)
{
	unsigned int iLower,iUpper,i;
	
	iLower = 0;
	iUpper = nCol-1;

	/*
	 * Make sure that rho is in the lookup table.
	 * If this is not the case return -1 or
	 * nCol. */
	if (T < Tarr[iLower])
	{
		return (-1);
	} else if (T > Tarr[iUpper])
	{
		return (nCol);
	}

	/* Do bisection. */
	while (iUpper-iLower > 1)
	{
		/* Compute the midpoint. */
		i = (iUpper + iLower) >> 1;

		if (T >= Tarr[i])
		{
			iLower = i;
		} else {
			iUpper = i;
		}
	}
	
	if (T == Tarr[0]) return 0;
	/* 
	 * Return Tmax-1 as the lower index, so the
	 * desired value is always bracketed by (i,i+1).
	 */
	if (T == Tarr[nCol-1])
	{
			return (nCol-2);
	}

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

	//printf("u= %g u(i,jmax)= %g iLower= %i\n",u, Lookup[iRow][nCol-1].u,iLower);
	return iLower;
}





void ANEOSInterpolateRhoU(ANEOSTable *table, double *rhoarr, double *tarr, double rho, double u, double *pP, double *pcs, double *pT, double *ps)
{
    double x,y;
    int i,j1,j2;

    double rho1,rho2;
    double frho,fu1,fu2;
    double u1,u2,u3,u4;
    double tp;
    /*
     * Since rho(i,j) = rho(i) and u(i,j) = u(j) we can do two 1D bisection to determine i and j.
     */
    i = ANEOSLookupRhoIndex(rhoarr, rho, table->nTableRho);
    //    printf("i=%i ",i);
    j1= ANEOSLookupUIndex(table->Lookup, u, i, table->nTableT);
      
    j2 = ANEOSLookupUIndex(table->Lookup, u, i+1, table->nTableT);

    /*assume bilinear interpolation of parallelogram,  https://math.stackexchange.com/questions/2007116/quadrilateral-interpolation*/


    rho1=rhoarr[i];
    rho2=rhoarr[i+1];

    u1=table->Lookup[i][j1].u;
    u2=table->Lookup[i][j1+1].u;
    u3=table->Lookup[i+1][j2].u;
    u4=table->Lookup[i+1][j2+1].u;

    rho1=logf(rho1);
    rho2=logf(rho2);
    u1=logf(u1);
    u2=logf(u2);
    u3=logf(u3);
    u4=logf(u4);
    
    frho=(logf(rho)-rho1)/(rho2-rho1);
    fu1=(logf(u)-u1)/(u2-u1);
    fu2=(logf(u)-u3)/(u4-u3);
 
    
    if (i < 0 || i >= table->nTableRho)
    {
      fprintf(stderr,"rho= %15.7E i= %i is outside of the table.\n", rho, i);
      //      assert(0);

    }



    tp = (1 - frho)*((1 - fu1) * logf(table->Lookup[i][j1].P) + fu1*logf(table->Lookup[i][j1+1].P))
      + frho*((1 - fu2) * logf(table->Lookup[i][j2].P) + fu2*logf(table->Lookup[i][j2+1].P));

    *pP = expf(tp);

    tp = (1 - frho)*((1 - fu1) * logf(table->Lookup[i][j1].cs) + fu1*logf(table->Lookup[i][j1+1].cs))
      + frho*((1 - fu2) * logf(table->Lookup[i][j2].cs) + fu2*logf(table->Lookup[i][j2+1].cs));

    *pcs = expf(tp);

    tp = (1 - frho)*((1 - fu1) * logf(table->Lookup[i][j1].S) + fu1*logf(table->Lookup[i][j1+1].S))
      + frho*((1 - fu2) * logf(table->Lookup[i][j2].S) + fu2*logf(table->Lookup[i][j2+1].S));

    *ps = expf(tp);

    tp = (1 - frho)*((1 - fu1) * logf(tarr[j1]) + fu1*logf(tarr[j1+1]))
      + frho*((1 - fu2) * logf(tarr[j2]) + fu2*logf(tarr[j2+1]));

    *pT = expf(tp);
	
}
/* old version in rho u space
void ANEOSInterpolateRhoU(ANEOSTable *table, double *rhoarr, double *tarr, double rho, double u, double *pP, double *pcs, double *pT, double *ps)
{
    double x,y;
    int i,j1,j2;

    double rho00,rho01,rho10,u00,u10,u01;
    double varA,varu,varv;
    double du1,du2;

     // Since rho(i,j) = rho(i) and u(i,j) = u(j) we can do two 1D bisection to determine i and j.
    
    i = ANEOSLookupRhoIndex(rhoarr, rho, table->nTableRho);
    //    printf("i=%i ",i);
    j1= ANEOSLookupUIndex(table->Lookup, u, i, table->nTableT);
      
    j2 = ANEOSLookupUIndex(table->Lookup, u, i+1, table->nTableT);

    //assume bilinear interpolation of parallelogram,  https://math.stackexchange.com/questions/2007116/quadrilateral-interpolation


    rho00=rhoarr[i];
    rho10=rhoarr[i];
    rho01=rhoarr[i+1];

    u00=table->Lookup[i][j1].u;
    u10=table->Lookup[i][j1+1].u;
    //    du1=u10-u00;
    u01=table->Lookup[i+1][j2].u;
    //    du2=table->Lookup[i+1][j2+1].u-u01;

    //    if(fabs(du1/du2-1)>0.2){fprintf(stderr,"Warning far from parallelogram, abs(du1/du2-1)='%lf' with rho= '%lf', u='%lf'\n", fabs(du1/du2-1), rho00, u00);}
    //fprintf(stderr,"%lf\n", fabs(du1/du2-1));

    varA=rho00*(u01-u10)+rho01*(u10-u00)+rho10*(u00-u01);
    varu=((rho01-rho00)*u-(u01-u00)*rho+rho00*u01-u00*rho01)/varA;
    varv=((rho00-rho10)*u-(u00-u10)*rho-rho00*u10+u00*rho10)/varA;
      
    
    if (i < 0 || i >= table->nTableRho)
    {
      fprintf(stderr,"rho= %15.7E i= %i is outside of the table.\n", rho, i);
      //      assert(0);

    }




     // We assume, that the grid is evenly spaced in rho and u.
     
    

    *pP = (1-varu)*(1-varv)*table->Lookup[i][j1].P + varu*(1.0-varv)*table->Lookup[i+1][j2].P +
        (1.0-varu)*varv*table->Lookup[i][j1+1].P + varu*varv*table->Lookup[i+1][j2+1].P;

    *pcs = (1-varu)*(1-varv)*table->Lookup[i][j1].cs + varu*(1.0-varv)*table->Lookup[i+1][j2].cs +
        (1.0-varu)*varv*table->Lookup[i][j1+1].cs + varu*varv*table->Lookup[i+1][j2+1].cs;

    *ps = (1-varu)*(1-varv)*table->Lookup[i][j1].S + varu*(1.0-varv)*table->Lookup[i+1][j2].S +
        (1.0-varu)*varv*table->Lookup[i][j1+1].S + varu*varv*table->Lookup[i+1][j2+1].S;

    *pT = (1-varu)*(1-varv)*tarr[j1] + varu*(1.0-varv)*tarr[j2] +
        (1.0-varu)*varv*tarr[j1+1] + varu*varv*tarr[j2+1];
	
}
*/
void ANEOSInterpolateRhoT(ANEOSTable *table, double *rhoarr, double *tarr, double rho, double T, double *pu, double *pP, double *pcs, double *ps)
{
    double x,y;
    int i,j;


    i = ANEOSLookupRhoIndex(rhoarr, rho, table->nTableRho);
    j = ANEOSLookupTIndex(tarr, T, table->nTableT);

    //    printf("i=%i ",i);
    
    //    printf("j=%i ",j);

    /*
     * Check if the data is in the lookup table.
     */
    if (i < 0 || i >= table->nTableRho)
    {
      fprintf(stderr,"rho= %15.7E i= %i is outside of the table.\n", rho, i);
      //      assert(0);

    }

    /*    if (j < 0 || j >= table->nTableT)
    {
      fprintf(stderr,"u= %15.7E j= %i is outside of the table.\n", u, j);
      //      assert(0);
      }*/

    /*
     * We assume, that the grid is evenly spaced in rho and u.
     */
    x = (rhoarr[i+1] - rho)/(rhoarr[i+1] - rhoarr[i]);
    y = (tarr[j+1] - T)/(tarr[j+1]-tarr[j]);

    *pP = x*y*table->Lookup[i][j].P + x*(1.0-y)*table->Lookup[i][j+1].P +
        (1.0-x)*y*table->Lookup[i+1][j].P + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].P;

    *pcs = x*y*table->Lookup[i][j].cs + x*(1.0-y)*table->Lookup[i][j+1].cs +
        (1.0-x)*y*table->Lookup[i+1][j].cs + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].cs;

    *pu = x*y*table->Lookup[i][j].u + x*(1.0-y)*table->Lookup[i][j+1].u +
        (1.0-x)*y*table->Lookup[i+1][j].u + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].u;

    *ps = x*y*table->Lookup[i][j].S + x*(1.0-y)*table->Lookup[i][j+1].S +
        (1.0-x)*y*table->Lookup[i+1][j].S + (1.0-x)*(1.0-y)*table->Lookup[i+1][j+1].S;
    
}

/*
 * Check the lookup table by comparing it to callaneos().
 */

