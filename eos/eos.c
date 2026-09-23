#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"

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

#ifdef EOS_ANEOS

    /* binary .spheos table: per-query unit conversion, table shared per node */
    {
      int imat = SphP[i].imat;
      EosTableState st;
      static int num_out_of_range = 0; /* per-rank diagnostic counter */

      if(imat < 0 || imat >= EosTableSpxNumMats)
        {
          printf("EOS_ANEOS: task %d particle %d has imat=%d, but only %d material IDs are mapped (EosTableMatIds)\n",
                 ThisTask, i, imat, EosTableSpxNumMats);
          endrun(1);
        }

      st = eos_table_evaluate_code(EosTableSpx,
                                   Particle_density_for_energy_i(i),
                                   SphP[i].InternalEnergyPred,
                                   (uint32_t)EosTableSpxMatId[imat], 1,
                                   EosTableSpxUnits);
      press                 = st.pressure;
      SphP[i].SoundSpeed    = st.soundSpeed;
      SphP[i].Temperature   = st.temperature;
      SphP[i].Entropy       = st.entropy;

      if(st.status != EOS_TABLE_SUCCESS)
        {
          num_out_of_range++;
          if(num_out_of_range <= 10 || num_out_of_range % 100000 == 0)
            printf("EOS_ANEOS: task %d particle %d (imat=%d): rho=%g u=%g code units -> status=%s (occurrence %d)\n",
                   ThisTask, i, imat, Particle_density_for_energy_i(i), SphP[i].InternalEnergyPred,
                   eos_table_status_string(st.status), num_out_of_range);
        }
    }

    if(press< 1e-15)
    {
      press=1e-15;
      SphP[i].SoundSpeed=1e-7;
    }
#endif   
    
    return press;
}

#if defined(MOONRELAX) && defined(EOS_ANEOS)
/* Isentropic relaxation pin (port of SPH-EXA's relaxIsentropic): while the
 * relaxation window is active, hold each gas particle on its initial
 * isentrope. On the first call the reference entropy adopts the table entropy
 * of the current (rho,u) state (NaN = uninitialized; retried next step when
 * the state is not located); afterwards InternalEnergy/InternalEnergyPred are
 * reset to u(rho, s0) from the table, immediately before the pressure
 * evaluation in the density loop. Energy conservation is intentionally
 * violated -- relaxation runs only, never production. */
void moonrelax_isentropic_pin(int i)
{
    if(All.RelaxIsentropic == 0 || All.RelaxTimescale <= 0.0 || All.Time >= All.RelaxUntil)
        return;

    int imat = SphP[i].imat;
    if(imat < 0 || imat >= EosTableSpxNumMats)
        return; /* get_pressure() below aborts on an unmapped imat anyway */

    double rho = Particle_density_for_energy_i(i);

    if(isnan(SphP[i].RelaxEntropy0))
    {
        EosTableState st = eos_table_evaluate_code(EosTableSpx, rho, SphP[i].InternalEnergyPred,
                                                   (uint32_t)EosTableSpxMatId[imat], 1, EosTableSpxUnits);
        int located = (st.status == EOS_TABLE_SUCCESS || st.status == EOS_TABLE_DENSITY_BELOW_RANGE ||
                       st.status == EOS_TABLE_DENSITY_ABOVE_RANGE || st.status == EOS_TABLE_ENERGY_BELOW_RANGE ||
                       st.status == EOS_TABLE_ENERGY_ABOVE_RANGE);
        if(st.hasEntropy && located)
            SphP[i].RelaxEntropy0 = st.entropy;
        else if(!st.hasEntropy && located)
        {
            printf("MOONRELAX: RelaxIsentropic requires an EOS table with entropy "
                   "(task=%d particle=%d imat=%d)\n", ThisTask, i, imat);
            endrun(1);
        }
        return; /* u untouched on the initialization step */
    }

    double uPinned = SphP[i].InternalEnergy;
    int status = eos_table_invert_energy_code(EosTableSpx, rho, SphP[i].RelaxEntropy0,
                                              (uint32_t)EosTableSpxMatId[imat], EosTableSpxUnits, &uPinned);
    if(status == EOS_TABLE_UNKNOWN_MATERIAL || status == EOS_TABLE_ENTROPY_UNAVAILABLE ||
       status == EOS_TABLE_INVALID_TABLE_STATE)
    {
        printf("MOONRELAX: isentropic pin failed: task=%d particle=%d imat=%d rho=%g s0=%g status=%s\n",
               ThisTask, i, imat, rho, (double)SphP[i].RelaxEntropy0, eos_table_status_string(status));
        endrun(1);
    }
    if(status != EOS_TABLE_INVALID_INPUT)
    {
        SphP[i].InternalEnergy     = uPinned;
        SphP[i].InternalEnergyPred = uPinned;
    }
    /* invalid input (e.g. non-positive density) keeps its current energy */
}
#endif





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







