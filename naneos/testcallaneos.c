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



int main(int argc, char *argv[])
{
//	char inFile[256] = "maneos.in";
    //    int i,j;
	
#if 0
    fprintf(stderr,"Reading %s.\n",inFile);

    if (ANEOSTableRead(table,inFile,601,501))
	{
			fprintf(stderr,"Failed reading %s.\n", inFile);
			assert(0);
	}

    fprintf(stderr,"Done.\n");
#endif

    fprintf(stderr, "Initializing material...\n");
	initaneos(inFile);
    
    /*
     * Initialize to some values.
     */



    fprintf(stderr, "Calling MANEOS.\n");
	callaneos(T, rho, iMat, &P, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH, &ion);
    fprintf(stderr, "\n");
	fprintf(stderr,"callaneos():\n");
	fprintf(stderr,"T=      %15.7E\n", T);
	fprintf(stderr,"rho =   %15.7E\n", rho);
	fprintf(stderr,"iMat=   %2i\n", iMat);
	fprintf(stderr,"P=      %15.7E\n", P);
	fprintf(stderr,"u=      %15.7E\n", u);
	fprintf(stderr,"s=      %15.7E\n", s);
	fprintf(stderr,"cv=     %15.7E\n", cv);
	fprintf(stderr,"dPdT=   %15.7E\n", dPdT);
	fprintf(stderr,"dPdrho= %15.7E\n", dPdrho);
	fprintf(stderr,"fkros=  %15.7E\n", fkros);
	fprintf(stderr,"cs=     %15.7E\n", cs);
	fprintf(stderr,"iPhase= %2i\n", iPhase);
	fprintf(stderr,"rhoL=   %15.7E\n", rhoL);
	fprintf(stderr,"rhoH=   %15.7E\n", rhoH);
	fprintf(stderr,"ion=    %15.7E\n", ion);
	fprintf(stderr,"\n");
	fprintf(stderr,"All quantities are in CGS.\n");

        u1=6.0634504E+10;



        u=0.1;
        while(fabs(u-u1)>0.0005E+10)
        {
          temp=(T1+T0)/2;
          callaneos(temp, rho, iMat, &P, &u, &s, &cv, &dPdT, &dPdrho, &fkros, &cs, &iPhase, &rhoL, &rhoH, &ion);          
          if(u>u1)
          {
            T1=temp;
            fprintf(stderr,"T1=      %15.7E\n", T1);
          }else{
            T0=temp;
            fprintf(stderr,"T0=      %15.7E\n", T0);
          }
        }
	fprintf(stderr,"T1=      %15.7E\n", T1);
        fprintf(stderr,"T0=      %15.7E\n", T0);
        fprintf(stderr,"Temp=      %15.7E\n", temp);
	fprintf(stderr,"rho =   %15.7E\n", rho);
	fprintf(stderr,"iMat=   %2i\n", iMat);
	fprintf(stderr,"P=      %15.7E\n", P);
	fprintf(stderr,"u=      %15.7E\n", u);
	fprintf(stderr,"s=      %15.7E\n", s);
	fprintf(stderr,"cv=     %15.7E\n", cv);
	fprintf(stderr,"dPdT=   %15.7E\n", dPdT);
	fprintf(stderr,"dPdrho= %15.7E\n", dPdrho);
	fprintf(stderr,"fkros=  %15.7E\n", fkros);
	fprintf(stderr,"cs=     %15.7E\n", cs);
	fprintf(stderr,"iPhase= %2i\n", iPhase);
	fprintf(stderr,"rhoL=   %15.7E\n", rhoL);
	fprintf(stderr,"rhoH=   %15.7E\n", rhoH);
	fprintf(stderr,"ion=    %15.7E\n", ion);
    return 0;
}
