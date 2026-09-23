#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#ifdef EOS_MANEOS
#include "../maneos/maneos.h"
#endif

/* Routines for gas equation-of-state terms (collects things like calculation of gas pressure)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* return the pressure of particle i */
double get_pressure(int i)
{
  MyFloat press = 0;
  press = GAMMA_MINUS1 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); /* ideal gas EOS (will get over-written it more complex EOS assumed) */
  
  //  press = 0.4 * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); 
#ifdef ISOTHERM_EQS 
  //  press = All.InitGasU * Particle_density_for_energy_i(i);
#endif // ISOTHERM_EQS 
  
  //  SphP[i].SoundSpeed=sqrt(1.4 * SphP[i].Pressure / Particle_density_for_energy_i(i));

  //    SphP[i].eosgamma= 0.4;
  //    SphP[i].eospsi=1.4*Particle_density_for_energy_i(i)*SphP[i].InternalEnergyPred;
       //        if(SphP[i].imat==1) SphP[i].SoundSpeed=sqrt(1.4 * SphP[i].Pressure / Particle_density_for_energy_i(i));




  /*#ifdef GLASS
  press = (1+beta)*GAMMA_MINUS1/(beta+GAMMA_MINUS1) * SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i);
#endif // GLASS
  */

#ifdef EOS_BATE
  if(SphP[i].Density<1.68e-7)
{
  press = 3.815e-5* Particle_density_for_energy_i(i);
  SphP[i].SoundSpeed=6.1767e-3;
}
else
{
  press = 0.0195463 * pow(SphP[i].Density, 1.4);
  SphP[i].SoundSpeed=sqrt(1.4*0.0195463* pow(SphP[i].Density, 0.4));
}
  
#endif  
    
    
#ifdef EOS_HELMHOLTZ
    /* pass the necessary quantities to D. Radice's wrappers for the Timms EOS */
    struct eos_input eos_in;
    struct eos_output eos_out;
    eos_in.rho  = SphP[i].Density;
    eos_in.eps  = SphP[i].InternalEnergyPred;
    eos_in.Ye   = SphP[i].Ye;
    eos_in.Abar = SphP[i].Abar;
    eos_in.temp = SphP[i].Temperature;
    int ierr = eos_compute(&eos_in, &eos_out);
    assert(!ierr);
    press              = eos_out.press;
    SphP[i].SoundSpeed = eos_out.csound;
    SphP[i].Temperature= eos_out.temp;
#endif

    
#ifdef EOS_ENFORCE_ADIABAT
    press = EOS_ENFORCE_ADIABAT * pow(SphP[i].Density, GAMMA);
#endif
    
    
    
    
#if defined(EOS_TRUELOVE_PRESSURE) || defined(TRUELOVE_CRITERION_PRESSURE)
    /* add an artificial pressure term to suppress fragmentation at/below the explicit resolution scale */
    double h_eff = DMAX(Get_Particle_Size(i), All.ForceSoftening[0]/2.8); /* need to include latter to account for inter-particle spacing << grav soft cases */
    /* standard finite-volume formulation of this (note there is some geometric ambiguity about whether there should be a "pi" in the equation below, but this 
        can be completely folded into the (already arbitrary) definition of NJeans, so we just use the latter parameter */
    double NJeans = 4; // set so that resolution = lambda_Jeans/NJeans -- fragmentation with Jeans/Toomre scales below this will be artificially suppressed now
    double xJeans = (NJeans * NJeans / GAMMA) * All.G * h_eff*h_eff * SphP[i].Density * SphP[i].Density * All.cf_afac1/All.cf_atime;
    if(xJeans>press) press=xJeans;
    SphP[i].SoundSpeed = sqrt(GAMMA * press / Particle_density_for_energy_i(i));
#endif
    
#ifdef EOS_TILLOTSON
    press = tillPressureSound(Mattable[SphP[i].imat],Particle_density_for_energy_i(i),SphP[i].InternalEnergy, &SphP[i].SoundSpeed);
    SphP[i].eospsi=tilldPdrho(Mattable[SphP[i].imat],Particle_density_for_energy_i(i),SphP[i].InternalEnergy);
    SphP[i].eosgamma=tilldPdu(Mattable[SphP[i].imat],Particle_density_for_energy_i(i),SphP[i].InternalEnergy);
        
    
    if(press< 1e-15) 
    {      
      press=1e-15;
      SphP[i].SoundSpeed=1e-7;
    }

#endif

#ifdef EOS_MANEOS
    if (SphP[i].InternalEnergyPred <2.0e-3)
    {
      SphP[i].InternalEnergy=2.0e-3;
      SphP[i].InternalEnergyPred=2.0e-3;
      }

    double rhotemp, utemp, Cv;
    rhotemp=Particle_density_for_energy_i(i);
    utemp=SphP[i].InternalEnergyPred;
#ifndef MOONRELAX  //to relax the initial condtion we have pressumed temperature profile
    ANEOSInterpolateRhoU(Mattable[SphP[i].imat], rhotemp, SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Temperature, &SphP[i].Entropy);

        //        if((rhotemp>21.)&&(rhotemp<26)&&(utemp>1.1)&&(utemp<5.5)) ANEOSInterpolateRhoU(Mattable[2], rhotemp, utemp, &press, &SphP[i].SoundSpeed, &SphP[i].Temperature, &SphP[i].Entropy, &SphP[i].eospsi, &SphP[i].eosgamma);
#else
    //  

    //    ANEOSInterpolateRhoU(Mattable[SphP[i].imat], rhotemp, SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Temperature, &SphP[i].Entropy);
      ANEOSInterpolateRhoT(Mattable[SphP[i].imat], rhotemp, SphP[i].Temperature, &SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Entropy); //from rho t to U P not working too coarse in T plane
#endif
    /*    if(P[i].imat==1)
    {
      SphP[i].Pressure=5.0;////0.4*SphP[i].InternalEnergyPred * Particle_density_for_energy_i(i); 
      SphP[i].SoundSpeed=sqrt(1.4 * SphP[i].Pressure / Particle_density_for_energy_i(i));
      }*/
    if(press< 1e-15) 
      {      
	press=1e-15;
	SphP[i].SoundSpeed=1e-7;
      }
#endif

#ifdef EOS_NANEOS

    double rhotemp, utemp, Cv;
    if (SphP[i].InternalEnergyPred <2.0e-3)
      {
	SphP[i].InternalEnergy=2.0e-3;
	SphP[i].InternalEnergyPred=2.0e-3;
      }
    rhotemp=Particle_density_for_energy_i(i);
    utemp=SphP[i].InternalEnergyPred;

#ifndef RHOT
    if (SphP[i].imat == 0)
      {ANEOSInterpolateRhoU(Mattable[0], rho0arr, t0arr, rhotemp, SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Temperature, &SphP[i].Entropy);}
    else
      {ANEOSInterpolateRhoU(Mattable[1], rho1arr, t1arr, rhotemp, SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Temperature, &SphP[i].Entropy);}
#else
    if (SphP[i].imat == 0)
      {ANEOSInterpolateRhoT(Mattable[0], rho0arr, t0arr, rhotemp, SphP[i].Temperature, &SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Entropy);}
    else
      {ANEOSInterpolateRhoT(Mattable[1], rho1arr, t1arr, rhotemp, SphP[i].Temperature, &SphP[i].InternalEnergyPred, &press, &SphP[i].SoundSpeed, &SphP[i].Entropy);}
#endif
      

      
    if(press< 1e-15) 
    {      
      press=1e-15;
      SphP[i].SoundSpeed=1e-7;
    }
#endif   
    
    return press;
}





/* trivial function to check if particle falls below the minimum allowed temperature */
void check_particle_for_temperature_minimum(int i)
{
    if(All.MinEgySpec)
    {
        if(SphP[i].InternalEnergy < All.MinEgySpec)
        {
            SphP[i].InternalEnergy = All.MinEgySpec;
            SphP[i].DtInternalEnergy = 0;
        }
    }
}



double INLINE_FUNC Particle_density_for_energy_i(int i)
{
#ifdef SPHEQ_DENSITY_INDEPENDENT_SPH
    return SphP[i].EgyWtDensity;
#endif
    return SphP[i].Density;
}




double INLINE_FUNC Particle_effective_soundspeed_i(int i)
{

#ifdef EOS_GENERAL
  return SphP[i].SoundSpeed;
#endif
#ifdef ISOTHERM_EQS 
    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
  
#endif     
    /* if nothing above triggers, then we resort to good old-fashioned ideal gas */
    /*
#ifdef EOS_BATE

  double rhotemp;
  rhotemp=Particle_density_for_energy_i(i);
  if(rhotemp>0.75)
  {
    return sqrt(2. * SphP[i].Pressure / Particle_density_for_energy_i(i));
  }
  else
  {
    return     sqrt(1.4 * SphP[i].Pressure / Particle_density_for_energy_i(i));
  }
  #endif*/

    return sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));


}







