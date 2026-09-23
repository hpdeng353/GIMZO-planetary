/*
 * A C wrapper to the Fortran functions of the M-ANEOS equation of state (Melosh 2007).
 * It is based on a C++ code that was provided by F. Benitez.
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "maneos.h"

/**
 * Initialize the ANEOS library.
 *
 * This is a wrapper arround aneosinit Fortran subroutine.
 */
void initaneos(char *matFilename) {
	static int done = FALSE;
//	char matFilename[256];

	if ( done ) {
		return;
	}

	done = TRUE;
	fprintf(stderr,"Initializing ANEOS with input file %s.\n", matFilename);

	// Use maneos.input for MANEOS or aneos.input for ANEOS
	//sprintf(matFilename, "maneos.in");
    assert(matFilename != NULL);
    assert(strlen(matFilename) > 0);
    
    // At some point it would be smart to check, if the desired file exists.
    //    aneosinit_( matFilename, strlen(matFilename) );
}

/**
 * Compute the ANEOS equation of state for given ( T, rho, mat ) values.
 *
 * This is a wrapper around ANEOSV Fortran subroutine.
 * The library will be initialized is not already done, so there is no
 * requirement to call initaneos() before this function.
 * All units are CGS-eV.
 *
 * @param T         Temperature (input)
 * @param rho       Density (input)
 * @param mat       Material number, as defined in the input file (input)
 * @param p         Pressure (output)
 * @param u         Specific internal energy (output)
 * @param S         Specific entropy (output)
 * @param cv        Specific heat capacity at constant volume (output)
 * @param dpdt      Temperature derivative of the pressure (output)
 * @param dpdrho    Density derivative of the pressure (output)
 * @param fkros     Rossland mean opacity (output)
 * @param cs        Speed of sound (output)
 * @param kpa       Phase (output)
 * @param rhoL      Density of the lower phase, for multi-phase states (output)
 * @param rhoH      Density of the higher phase, for multi-phase states (output)
 * @param ion       Ionisation number (output)
 */
void callaneos (
	const double T, const double rho, const int mat,
	double* p, double* u, double* S, double* cv,
	double* dpdt, double* dpdrho, double* fkros, double* cs,
	int* kpa, double* rhoL, double* rhoH, double* ion
) {
	const int n = 1;

//	initaneos();
    // Make sure that the library was initalized.
//    assert(done);
//	aneosv_ (
//		&n, &T, &rho, &mat, p, u, S, cv, dpdt, dpdrho, fkros,
//		cs, kpa, rhoL, rhoH, ion
//	);
}


