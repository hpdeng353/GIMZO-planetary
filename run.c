#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part (adding/removing calls, re-ordering some routines, and 
 * adding hooks to new elements such as particle splitting, as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void run(void)
{
    CPU_Step[CPU_MISC] += measure_time();

#ifdef MARKNGB
    int i;
    for(i = 0; i < NumPart; i++)
    {
      P[i].ngbto=0;
    }
#endif    
    if(RestartFlag != 1)		/* need to compute forces at initial synchronization time, unless we restarted from restart files */
    {
        output_log_messages();
        
        domain_Decomposition(0, 0, 0);
        
        set_non_standard_physics_for_current_time();
        
        compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */
        
        compute_hydro_densities_and_forces();	/* densities, gradients, & hydro-accels for synchronous particles */
        
        calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */
    }

    while(1)			/* main timestep iteration loop */
    {
        compute_statistics();	/* regular statistics outputs (like total energy) */
        
        write_cpu_log();		/* output some CPU usage log-info (accounts for everything needed up to the current sync-point) */
        
        if(All.Ti_Current >= TIMEBASE)	/* check whether we reached the final time */
        {
            if(ThisTask == 0)
                printf("\nFinal time=%g reached. Simulation ends.\n", All.TimeMax);
            
            restart(0);		/* write a restart file to allow continuation of the run for a larger value of TimeMax */
            
            if(All.Ti_lastoutput != All.Ti_Current)	/* make a snapshot at the final time in case none has produced at this time */
                savepositions(All.SnapshotFileCount++);	/* this will be overwritten if All.TimeMax is increased and the run is continued */
            
            break;
        }


#ifdef DISKIC
        //        All.G*=0.99;
        eq_relax5();
#endif



/* MOONRELAX relaxation now acts inside do_the_kick (moonrelax_modify_kick),
   scaled by each particle's own active timestep -- no per-loop call here. */

        find_timesteps();		/* find-timesteps */
        
        do_first_halfstep_kick();	/* half-step kick at beginning of timestep for synchronous particles */
        /*#ifdef MARKNGB

        for(i = 0; i < NumPart; i++)
        {
          if(P[i].ID==1) P[i].ngbto=100;
        }
        #endif    */    
        find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
                                             * If needed, this function will also write an output file
                                             * at the desired time.
                                             */
#ifdef MARKNGB
        int i;
        for(i = 0; i < NumPart; i++)
        {
          P[i].ngbto=0;
        }
#endif

#ifdef GDSPH
        int i;
        for(i = 0; i < NumPart; i++)
        {
          SphP[i].temprho=SphP[i].Density;
        }
#endif
	
        output_log_messages();	/* write some info to log-files */
        
        set_non_standard_physics_for_current_time();	/* update auxiliary physics for current time */
        
        if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)	/* check whether we have a big step */
        {
            domain_Decomposition(0, 0, 1);	/* do domain decomposition if step is big enough, and set new list of active particles  */
        }
        else
        {
            force_update_tree();	/* update tree dynamically with kicks of last step so that it can be reused */
            
            make_list_of_active_particles();	/* now we can set the new chain list of active particles */
        }
        
        compute_grav_accelerations();	/* compute gravitational accelerations for synchronous particles */

#ifdef GALSF_SUBGRID_DMDISPERSION
        // Need to figure out how frequently we calculate this; below is pretty rough //
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
        {
            disp_density(); /* compute the DM velocity dispersion around gas particles every 20 PM steps, should be sufficient */
        }
#endif

#if (defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE))
        /* flag particles which will be feedback centers, so kernel lengths can be computed for them */
        determine_where_SNe_occur();
#endif
        
        compute_hydro_densities_and_forces();	/* densities, gradients, & hydro-accels for synchronous particles */
        
        do_second_halfstep_kick();	/* this does the half-step kick at the end of the timestep */
        
        calculate_non_standard_physics();	/* source terms are here treated in a strang-split fashion */

        /* Check whether we need to interrupt the run */
        int stopflag = 0;
#ifdef IO_REDUCED_MODE
        if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
        if(ThisTask == 0)
        {
            FILE *fd;
            char stopfname[1000];
            sprintf(stopfname, "%sstop", All.OutputDir);
            if((fd = fopen(stopfname, "r")))	/* Is the stop-file present? If yes, interrupt the run. */
            {
                fclose(fd);
                stopflag = 1;
                unlink(stopfname);
            }
            
            if(CPUThisRun > 0.85 * All.TimeLimitCPU)	/* are we running out of CPU-time ? If yes, interrupt run. */
            {
                printf("reaching time-limit. stopping.\n");
                stopflag = 2;
            }
        }
        
        MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if(stopflag)
        {
            restart(0);		/* write restart file */
            MPI_Barrier(MPI_COMM_WORLD);
            
            if(stopflag == 2 && ThisTask == 0)
            {
                FILE *fd;
                char contfname[1000];
                sprintf(contfname, "%scont", All.OutputDir);
                if((fd = fopen(contfname, "w")))
                    fclose(fd);
                
                if(All.ResubmitOn)
                    execute_resubmit_command();
            }
            return;
        }
        
        if(ThisTask == 0)
        {
            /* is it time to write one of the regularly space restart-files? */
            if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
            {
                All.TimeLastRestartFile = CPUThisRun;
                stopflag = 3;
            }
            else
                stopflag = 0;
        }
        
        MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if(stopflag == 3)
        {
            restart(0);		/* write an occasional restart file */
            stopflag = 0;
            All.TimeLastRestartFile += report_time();
        }
        
        set_random_numbers();	/* draw a new list of random numbers */
        
        report_memory_usage(&HighMark_run, "RUN");
        
    }
    
}



void set_non_standard_physics_for_current_time(void)
{
#ifdef SET_alpha
  int i;
  double radius=0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
        {
          radius=sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1]);
          if(All.Time<0.1){
	    //	    SphP[i].Eta_ShearViscosity=0.1*0.02*0.02*pow(radius,0.5)*SphP[i].Density;
	    SphP[i].Eta_ShearViscosity=0.1*0.02*0.02*pow(radius,-1.5)*SphP[i].Density;
	    //            SphP[i].Eta_ShearViscosity=0.06*0.05*0.05*pow(radius,0.4)*SphP[i].Density;
	    //            SphP[i].Zeta_BulkViscosity=1*0.02*0.02*pow(radius,-1.5)*SphP[i].Density;
          }
        }
    }
#endif
    
}



void calculate_non_standard_physics(void)
{
#ifdef PARTICLE_EXCISION
    apply_excision();
#endif
        
    
    
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
#ifdef BH_WIND_SPAWN
    if(GlobNumForceUpdate > All.TreeDomainUpdateFrequency * All.TotNumPart)
    {
        spawn_bh_wind_feedback();
        rearrange_particle_sequence();
        force_treebuild(NumPart, NULL);
    }
#endif
#endif // ifdef BLACK_HOLES or GALSF_SUBGRID_VARIABLEVELOCITY


#ifdef INFALL
     if((All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) && (ThisTask ==0))
    {
      InfallID =rand() % (240000 + 1 - 80000) + 80000;      
      //InfallID =rand() % (32768 + 1 - 1) + 1;      
    }
     MPI_Bcast(&InfallID, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif



#if (defined(BETACOOL) || defined(SET_Ohmic) )
    NextParticle = FirstActiveParticle;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        while(1)
        {
            int i, exitFlag = 0;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
            {
                if(NextParticle<0) {exitFlag = 1;} else {i=NextParticle; NextParticle=NextActiveParticle[NextParticle];}
            }
            if(exitFlag) {break;}

            /* here apply any conditional statements about whether we should or should not enter the cooling loop */
            if(P[i].Type != 0) {continue;} /* only gas cools */
            if(P[i].Mass <= 0) {continue;} /* only non-zero mass particles cool */
            //            if((SphP[i].Density >1.6e-3)&&(P[i].r2star2>400)){continue;}

#ifdef SET_Ohmic
            if(All.Time<0.5)
            {
              //              SphP[i].Eta_MHD_OhmicResistivity_Coeff=1.11*2e-4/20*pow(P[i].r2star2,0.75);
              //SphP[i].Eta_MHD_OhmicResistivity_Coeff=0.0001/100.*pow(P[i].r2star2,0.5);
              //SphP[i].Eta_MHD_OhmicResistivity_Coeff=0.0001/20.*pow(P[i].r2star2,0.5); 
              //              SphP[i].Eta_MHD_OhmicResistivity_Coeff=0.0004/20.;
              // SphP[i].Eta_MHD_OhmicResistivity_Coeff=1.11*0.032/100.*pow(P[i].r2star2,0.125);
            }
#endif


            //            if(SphP[i].Density >1.6e-3){continue;}
            //            if(P[i].Mass >5.e-9){continue;}
            //            if(P[i].Mass >5.e-8){continue;}
            //            if(P[i].Mass >2.e-9){continue;}
            double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
            double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
            double newu=0, z_offset=0,deltau=0;
            newu = SphP[i].InternalEnergy;
#ifdef ANALYTIC_STAR
            P[i].r2star2=P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2];
#endif

            /*#ifdef BINARY
            P[i].r2star2=P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2];
            #endif*/
#ifdef BETACOOL

#ifndef SHEARING_BOX
            newu = SphP[i].InternalEnergy - SphP[i].InternalEnergy /All.betacool * dtime * pow(P[i].r2star2,-0.75);

#ifdef IRR
            newu = SphP[i].InternalEnergy - (SphP[i].InternalEnergy-All.InitGasU*pow(P[i].r2star2, -0.25))/All.betacool * dtime * pow(P[i].r2star2,-0.75);
#endif

#else
#if (SHEARING_BOX!=4)
            //            newu = SphP[i].InternalEnergy - SphP[i].InternalEnergy/All.betacool * dtime;
            newu = SphP[i].InternalEnergy  - (SphP[i].InternalEnergy-(1/(GAMMA -1.)))/6.28 * dtime;
#else 


            z_offset=fabs(P[i].Pos[2]-boxHalf_Z);
#ifndef SBETACOOL
            if(z_offset<20.) newu = SphP[i].InternalEnergy - (SphP[i].InternalEnergy-(1/(GAMMA -1.)))/All.betacool * dtime;
            //            newu = SphP[i].InternalEnergy - SphP[i].InternalEnergy/All.betacool * dtime;
#else
            if(z_offset<20.) newu = SphP[i].InternalEnergy - SphP[i].DeltaU/All.betacool * dtime; //apply to stratified shearing box body
            
#endif
#endif
#endif            

              if(newu <0.9*SphP[i].InternalEnergy)
              {
                SphP[i].InternalEnergy *= 0.9;
                printf("warning cool floor");
              }
              else {SphP[i].InternalEnergy = newu;}
            
              if(SphP[i].InternalEnergy < 1e-30){SphP[i].InternalEnergy=1e-30;}
              SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
              SphP[i].Pressure = get_pressure(i);
              
#endif              


        } /* while bracket */
    } /* omp bracket */

#endif
} 






void compute_statistics(void)
{
    if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
    {
#if !defined(EVALPOTENTIAL)          // DAA: compute_potential is not defined if EVALPOTENTIAL is on... check!
#ifdef COMPUTE_POTENTIAL_ENERGY
        compute_potential();
#endif
#endif
#ifndef IO_REDUCED_MODE
        energy_statistics();	/* compute and output energy statistics */
#endif
        
        All.TimeLastStatistics += All.TimeBetStatistics;
    }
}



void execute_resubmit_command(void)
{
    char buf[1000];
    sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
    system(buf);
#endif
    if(ThisTask ==0) {printf("[%s] is the resubmit commd\n",buf);}
}



/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system drifts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, i, prev;
  integertime dt_bin, ti_next_for_bin, ti_next_kick, ti_next_kick_global;
  int highest_active_bin, highest_occupied_bin;
  double timeold;

  timeold = All.Time;

  All.NumCurrentTiStep++;	/* we are now moving to the next sync point */

  /* find the next kick time */
  for(n = 0, ti_next_kick = TIMEBASE, highest_occupied_bin = 0; n < TIMEBINS; n++)
    {
      if(TimeBinCount[n])
	{
	  if(n > 0)
	    {
	      highest_occupied_bin = n;
	      dt_bin = (((integertime) 1) << n);
	      ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
	    }
	  else
	    {
	      dt_bin = 0;
	      ti_next_for_bin = All.Ti_Current;
	    }

	  if(ti_next_for_bin < ti_next_kick)
	    ti_next_kick = ti_next_for_bin;
	}
    }

  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  while(ti_next_kick_global >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      All.Ti_Current = All.Ti_nextoutput;

      if(All.ComovingIntegrationOn)
          All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      else
          All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

      set_cosmo_factors_for_current_time();


      move_particles(All.Ti_nextoutput);

      CPU_Step[CPU_DRIFT] += measure_time();

#ifdef OUTPUTPOTENTIAL
#if !defined(EVALPOTENTIAL) || (defined(EVALPOTENTIAL) && defined(RECOMPUTE_POTENTIAL_ON_OUTPUT))
      domain_Decomposition(0, 0, 0);

      compute_potential();
#endif
#endif


#ifndef IO_REDUCED_MODE
      mpi_printf("\n\n\nI found the last snapshot call...\n\n\n");
#endif
        savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }


  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  set_cosmo_factors_for_current_time();
#ifdef SHEARING_BOX
    calc_shearing_box_pos_offset();
#endif


  All.TimeStep = All.Time - timeold;

  /* mark the bins that will be active */
  for(n = 1, TimeBinActive[0] = 1, NumForceUpdate = TimeBinCount[0], highest_active_bin = 0; n < TIMEBINS;
      n++)
    {
      dt_bin = (((integertime) 1) << n);
      if((ti_next_kick_global % dt_bin) == 0)
	{
	  TimeBinActive[n] = 1;
	  NumForceUpdate += TimeBinCount[n];
	  if(TimeBinCount[n])
	    highest_active_bin = n;
	}
      else
	TimeBinActive[n] = 0;
    }

  sumup_large_ints(1, &NumForceUpdate, &GlobNumForceUpdate);
  MPI_Allreduce(&highest_active_bin, &All.HighestActiveTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&highest_occupied_bin, &All.HighestOccupiedTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(GlobNumForceUpdate == All.TotNumPart)
    {
      Flag_FullStep = 1;
      if(All.HighestActiveTimeBin != All.HighestOccupiedTimeBin)
	terminate("Something is wrong with the time bins.\n");
    }
  else
    Flag_FullStep = 0;




  /* move the new set of active/synchronized particles */
  /* Note: We do not yet call make_list_of_active_particles(), since we
   * may still need to old list in the dynamic tree update
   */
  for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
      if(TimeBinActive[n])
	{
	  for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
	    {
	      drift_particle(i, All.Ti_Current);
	    }
	}
    }

}


void make_list_of_active_particles(void)
{
    int i, n, prev;
    /* make a link list with the particles in the active time bins */
    FirstActiveParticle = -1;
    
    for(n = 0, prev = -1; n < TIMEBINS; n++)
    {
        if(TimeBinActive[n])
        {
            for(i = FirstInTimeBin[n]; i >= 0; i = NextInTimeBin[i])
            {
                if(P[i].Mass <= 0)
                    continue;
                
                if(prev == -1)
                    FirstActiveParticle = i;
                
                if(prev >= 0)
                    NextActiveParticle[prev] = i;
                
                prev = i;
            }
        }
    }
    
    if(prev >= 0)
        NextActiveParticle[prev] = -1;
}





/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
integertime find_next_outputtime(integertime ti_curr)
{
  int i, iter = 0;
  integertime ti, ti_next;
  double next, time;

  DumpFlag = 1;
  ti_next = -1;


  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  if(time >= All.TimeBegin && time <= All.TimeMax)
	    {
	      if(All.ComovingIntegrationOn)
		ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	      else
		ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

	      if(ti >= ti_curr)
		{
		  if(ti_next == -1)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }

		  if(ti_next > ti)
		    {
		      ti_next = ti;
		      DumpFlag = All.OutputListFlag[i];
		      if(i > All.SnapshotFileCount)
			All.SnapshotFileCount = i;
		    }
		}
	    }
	}
    }
  else
    {
      if(All.ComovingIntegrationOn)
	{
	  if(All.TimeBetSnapshot <= 1.0)
	    {
	      printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      else
	{
	  if(All.TimeBetSnapshot <= 0.0)
	    {
	      printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
	      endrun(13123);
	    }
	}
      time = All.TimeOfFirstSnapshot;

      iter = 0;

      while(time < All.TimeBegin)
	{
	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(110);
	    }
	}
      while(time <= All.TimeMax)
	{
	  if(All.ComovingIntegrationOn)
	    ti = (integertime) (log(time / All.TimeBegin) / All.Timebase_interval);
	  else
	    ti = (integertime) ((time - All.TimeBegin) / All.Timebase_interval);

	  if(ti >= ti_curr)
	    {
	      ti_next = ti;
	      break;
	    }

	  if(All.ComovingIntegrationOn)
	    time *= All.TimeBetSnapshot;
	  else
	    time += All.TimeBetSnapshot;

	  iter++;

	  if(iter > 1000000)
	    {
	      printf("Can't determine next output time.\n");
	      endrun(111);
	    }
	}
    }


  if(ti_next == -1)
    {
      ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

      if(ThisTask == 0)
	printf("\nThere is no valid time for a further snapshot file.\n");
    }
  else
    {
      if(All.ComovingIntegrationOn)
	next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
      else
	next = All.TimeBegin + ti_next * All.Timebase_interval;

      if(ThisTask == 0)
	printf("\nSetting next time for snapshot file to Time_next= %g  (DumpFlag=%d)\n\n", next, DumpFlag);

    }

  return ti_next;
}




/*! This routine writes for every synchronisation point in the timeline information to two log-files:
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdTimebins we inform about the distribution of particles over the timebins, and which timebins are active on this step.
 * code is stored.
 */
void output_log_messages(void)
{
  double z;
  int i, j;
  long long tot, tot_sph;
  long long tot_count[TIMEBINS];
  long long tot_count_sph[TIMEBINS];
  long long tot_cumulative[TIMEBINS];
  int weight, corr_weight;
  double sum, avg_CPU_TimeBin[TIMEBINS], frac_CPU_TimeBin[TIMEBINS];

  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinCountSph, tot_count_sph);

    if(ThisTask == 0)
    {
        if(All.ComovingIntegrationOn)
        {
            z = 1.0 / (All.Time) - 1;
#ifndef IO_REDUCED_MODE
            fprintf(FdInfo, "\nSync-Point %d, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
                    All.NumCurrentTiStep, All.Time, z,
                    (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
                    All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
            fflush(FdInfo);
            fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
                    All.NumCurrentTiStep, All.Time, z, All.TimeStep,
                    log(All.Time) - log(All.Time - All.TimeStep));
#endif
            printf("\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
                   All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
        }
        else
        {
#ifndef IO_REDUCED_MODE
            fprintf(FdInfo, "\nSync-Point %d, Time: %g, Nf = %d%09d, Systemstep: %g\n", All.NumCurrentTiStep,
                    All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
                    All.TimeStep);
            fflush(FdInfo);
            fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
                    All.TimeStep);
#endif
            printf("\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
        }

        for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++)
            tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];


      for(i = 0; i < TIMEBINS; i++)
	{
	  for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
	    sum += All.CPU_TimeBinMeasurements[i][j];
	  if(All.CPU_TimeBinCountMeasurements[i])
	    avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
	  else
	    avg_CPU_TimeBin[i] = 0;
	}

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count[i] > 0; i--, weight *= 2)
	{
	  if(weight > 1)
	    corr_weight = weight / 2;
	  else
	    corr_weight = weight;

	  frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
	  sum += frac_CPU_TimeBin[i];
	}

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count[i] > 0; i--)
	{
	  if(sum)
	    frac_CPU_TimeBin[i] /= sum;
	}


        printf("Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin,"Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
#endif
        for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
            if(tot_count_sph[i] > 0 || tot_count[i] > 0)
            {
                printf(" %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
                       TimeBinActive[i] ? 'X' : ' ',
                       i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
                       i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
                       (i == All.HighestActiveTimeBin) ? '<' : ' ',
                       (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
                       avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);
#ifndef IO_REDUCED_MODE
                fprintf(FdTimebin,
                        " %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
                        TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
                        i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
                        (i == All.HighestActiveTimeBin) ? '<' : ' ',
                        (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
                        avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);
#endif
                if(TimeBinActive[i])
                {
                    tot += tot_count[i];
                    tot_sph += tot_count_sph[i];
                }
            }
        printf("               ------------------------\n");
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin, "               ------------------------\n");
#endif
        {
            printf("Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
#ifndef IO_REDUCED_MODE
            fprintf(FdTimebin, "Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
#endif
        }
#ifndef IO_REDUCED_MODE
        fprintf(FdTimebin, "\n");
        fflush(FdTimebin);
#endif
    }
    
  output_extra_log_messages();
}




void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;
  int i;

  CPU_Step[CPU_MISC] += measure_time();

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS; i++)
    CPU_Step[0] += CPU_Step[i];

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	avg_CPU_Step[i] /= NTask;

      put_symbol(0.0, 1.0, '#');

      for(i = 1, tsum = 0.0; i < CPU_PARTS; i++)
	{
	  if(max_CPU_Step[i] > 0)
	    {
	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
	      tsum += t1 - t0;

	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
	      tsum += t1 - t0;
	    }
	}

      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');

#ifndef IO_REDUCED_MODE
      fprintf(FdBalance, "Step=%7d  sec=%10.3f  Nf=%2d%09d  %s\n", All.NumCurrentTiStep, max_CPU_Step[0],
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), CPU_String);
      fflush(FdBalance);
#endif
        
      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
	{
	  All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
	  memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0],
		  &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
		  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
	}

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements
							    [All.HighestActiveTimeBin]++] = max_CPU_Step[0];
    }

  CPUThisRun += CPU_Step[0];

  for(i = 0; i < CPU_PARTS; i++)
    CPU_Step[i] = 0;

#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	All.CPU_Sum[i] += avg_CPU_Step[i];

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
	      "total         %10.2f  %5.1f%%\n"
	      "treegrav      %10.2f  %5.1f%%\n"
	      "   treebuild  %10.2f  %5.1f%%\n"
	      "   treeupdate %10.2f  %5.1f%%\n"
	      "   treewalk   %10.2f  %5.1f%%\n"
	      "   treecomm   %10.2f  %5.1f%%\n"
	      "   treeimbal  %10.2f  %5.1f%%\n"
#ifdef ADAPTIVE_GRAVSOFT_FORALL
	      "adaptgrav     %10.2f  %5.1f%%\n"
	      "   agsdensity %10.2f  %5.1f%%\n"
	      "   agscomm    %10.2f  %5.1f%%\n"
	      "   agsimbal   %10.2f  %5.1f%%\n"
#endif
	      "pmgrav        %10.2f  %5.1f%%\n"
	      "hydro         %10.2f  %5.1f%%\n"
	      "   density    %10.2f  %5.1f%%\n"
	      "   denscomm   %10.2f  %5.1f%%\n"
	      "   densimbal  %10.2f  %5.1f%%\n"
	      "   hydrofrc   %10.2f  %5.1f%%\n"
	      "   hydcomm    %10.2f  %5.1f%%\n"
	      "   hydmisc    %10.2f  %5.1f%%\n"
	      "   hydnetwork %10.2f  %5.1f%%\n"
	      "   hydimbal   %10.2f  %5.1f%%\n"
	      "   hmaxupdate %10.2f  %5.1f%%\n"
	      "domain        %10.2f  %5.1f%%\n"
	      "potential     %10.2f  %5.1f%%\n"
	      "predict       %10.2f  %5.1f%%\n"
	      "kicks         %10.2f  %5.1f%%\n"
	      "i/o           %10.2f  %5.1f%%\n"
	      "peano         %10.2f  %5.1f%%\n"
	      "sfrcool       %10.2f  %5.1f%%\n"
	      "blackholes    %10.2f  %5.1f%%\n"
	      "fof/subfind   %10.2f  %5.1f%%\n"
#ifdef GRAIN_FLUID
          "grains        %10.2f  %5.1f%%\n"
#endif
          "gas_return    %10.2f  %5.1f%%\n"
          "snII_fb_loop  %10.2f  %5.1f%%\n"
          "hII_fb_loop   %10.2f  %5.1f%%\n"
          "localwindkik  %10.2f  %5.1f%%\n"
          "misc          %10.2f  %5.1f%%\n",
              
    All.CPU_Sum[CPU_ALL], 100.0,
    All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2] + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
              + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2] + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE]
              + All.CPU_Sum[CPU_TREEMISC],
    (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2] + All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]
              + All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2] + All.CPU_Sum[CPU_TREEBUILD] + All.CPU_Sum[CPU_TREEUPDATE] + All.CPU_Sum[CPU_TREEMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEBUILD], (All.CPU_Sum[CPU_TREEBUILD]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEUPDATE], (All.CPU_Sum[CPU_TREEUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2], (All.CPU_Sum[CPU_TREEWALK1] + All.CPU_Sum[CPU_TREEWALK2]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV], (All.CPU_Sum[CPU_TREESEND] + All.CPU_Sum[CPU_TREERECV]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2], (All.CPU_Sum[CPU_TREEWAIT1] + All.CPU_Sum[CPU_TREEWAIT2]) / All.CPU_Sum[CPU_ALL] * 100,
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    All.CPU_Sum[CPU_AGSDENSCOMPUTE] + All.CPU_Sum[CPU_AGSDENSWAIT] + All.CPU_Sum[CPU_AGSDENSCOMM] + All.CPU_Sum[CPU_AGSDENSMISC],
              (All.CPU_Sum[CPU_AGSDENSCOMPUTE] + All.CPU_Sum[CPU_AGSDENSWAIT] + All.CPU_Sum[CPU_AGSDENSCOMM] + All.CPU_Sum[CPU_AGSDENSMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSCOMPUTE], (All.CPU_Sum[CPU_AGSDENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSCOMM], (All.CPU_Sum[CPU_AGSDENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_AGSDENSWAIT], (All.CPU_Sum[CPU_AGSDENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
    All.CPU_Sum[CPU_MESH], (All.CPU_Sum[CPU_MESH]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT] + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
              + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
              + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_HYDNETWORK],
    (All.CPU_Sum[CPU_DENSCOMPUTE] + All.CPU_Sum[CPU_DENSWAIT] + All.CPU_Sum[CPU_DENSCOMM] + All.CPU_Sum[CPU_DENSMISC]
              + All.CPU_Sum[CPU_HYDCOMPUTE] + All.CPU_Sum[CPU_HYDWAIT] + All.CPU_Sum[CPU_TREEHMAXUPDATE]
              + All.CPU_Sum[CPU_HYDCOMM] + All.CPU_Sum[CPU_HYDMISC] + All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSCOMPUTE], (All.CPU_Sum[CPU_DENSCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSCOMM], (All.CPU_Sum[CPU_DENSCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DENSWAIT], (All.CPU_Sum[CPU_DENSWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDCOMPUTE], (All.CPU_Sum[CPU_HYDCOMPUTE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDCOMM], (All.CPU_Sum[CPU_HYDCOMM]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDMISC], (All.CPU_Sum[CPU_HYDMISC]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDNETWORK], (All.CPU_Sum[CPU_HYDNETWORK]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HYDWAIT], (All.CPU_Sum[CPU_HYDWAIT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TREEHMAXUPDATE], (All.CPU_Sum[CPU_TREEHMAXUPDATE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DOMAIN], (All.CPU_Sum[CPU_DOMAIN]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_POTENTIAL], (All.CPU_Sum[CPU_POTENTIAL]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_DRIFT], (All.CPU_Sum[CPU_DRIFT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_TIMELINE], (All.CPU_Sum[CPU_TIMELINE]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_SNAPSHOT], (All.CPU_Sum[CPU_SNAPSHOT]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_PEANO], (All.CPU_Sum[CPU_PEANO]) / All.CPU_Sum[CPU_ALL] * 100,
    0.,0.,
    0.,0.,
    0.,0.,
#ifdef GRAIN_FLUID
    All.CPU_Sum[CPU_DRAGFORCE], (All.CPU_Sum[CPU_DRAGFORCE]) / All.CPU_Sum[CPU_ALL] * 100,
#endif
    All.CPU_Sum[CPU_GASRETURN], (All.CPU_Sum[CPU_GASRETURN]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_SNIIHEATING], (All.CPU_Sum[CPU_SNIIHEATING]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_HIIHEATING], (All.CPU_Sum[CPU_HIIHEATING]) / All.CPU_Sum[CPU_ALL] * 100,
    All.CPU_Sum[CPU_LOCALWIND], (All.CPU_Sum[CPU_LOCALWIND]) / All.CPU_Sum[CPU_ALL] * 100,

    All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
        
    fprintf(FdCPU, "\n");
    fflush(FdCPU);
    }
}



void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j < 0)
    j = 0;
  if(i >= CPU_STRING_LEN)
    i = CPU_STRING_LEN;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}



#ifndef IO_REDUCED_MODE
/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy,
	      "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	      All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
	      SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
	      SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
	      SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
	      SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
	      SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
	      SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
	      SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
	      SysState.MassComp[5]);

      fprintf(FdEnergy," \n");
      fflush(FdEnergy);
    }
}
#endif



void output_extra_log_messages(void)
{
    
    
    
    if(ThisTask == 0)
    {
    }
}




void check_particles_info(const char *func, const char *file, int linenr)
{
  int i,k,vok=0,pok=0,vsph=0;
  double vv;

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Checking particle data (function %s in file %s at line %d) ...\n",func,file,linenr);

  for(i = 0; i < NumPart; i++)
    {
      
      for(k = 0; k < 3; k++)
	{
	  if( P[i].Vel[k] > -1e8 && P[i].Vel[k] < 1e8)
	    {
	      vv = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
	      if(vv > 15000)
		{
		  printf("task=%d: WARNING: Large velocity for particle %d ID %llu v[%d]=%g, renormalizing it !!\n", ThisTask, i,
			 (unsigned long long) P[i].ID,k,vv);
		  fflush(stdout);
                  P[i].Vel[0] = P[i].Vel[0] / vv * 10000;
                  P[i].Vel[1] = P[i].Vel[1] / vv * 10000;
                  P[i].Vel[2] = P[i].Vel[2] / vv * 10000;
		}
	      vok++;
	    }
	  else
	    {
	      printf("task=%d:  strange value in velocity in for particle %d ID %llu , type=%d, mass=%g, v[%d]=%g\n", ThisTask, i,
		     (unsigned long long) P[i].ID,P[i].Type,P[i].Mass,k,P[i].Vel[k]);
		     fflush(stdout);
		     endrun(712401);
	    }

	  if( P[i].Pos[k] > -10000 && P[i].Pos[k] < All.BoxSize + 10000)
	    pok++;
	  else
	    {
	      printf("task=%d:  strange value in position in for particle %d ID %llu x[%d]=%g\n", ThisTask, i,
		     (unsigned long long) P[i].ID,k,P[i].Pos[k]);
		     fflush(stdout);
		     endrun(712402);
	    }
	}

      if(P[i].Type == 0)
        {
          if( (SphP[i].InternalEnergyPred > -1e20 && SphP[i].InternalEnergyPred < 1e20))
            vsph++;
          else
            printf("task=%d: Particle id=%llu, m=%e strange InternalEnergy value in hydro: %e\n",
		   ThisTask,(unsigned long long) P[i].ID,P[i].Mass,SphP[i].InternalEnergyPred);

          if( (SphP[i].DtInternalEnergy > -1e20 && SphP[i].DtInternalEnergy < 1e20))
            vsph++;
          else
            printf("task=%d: Particle id=%llu, m=%e strange DtInternalEnergy value in hydro: %e\n",
		   ThisTask,(unsigned long long) P[i].ID,P[i].Mass,SphP[i].DtInternalEnergy);

          if( (SphP[i].Pressure  > -1e20 && SphP[i].Pressure <1e20))
            vsph++;
          else
            printf("task=%d: Particle id=%llu,m=%e strange Pressure value in hydro: %e\n",
		   ThisTask,(unsigned long long) P[i].ID,P[i].Mass,SphP[i].Pressure);
        }

      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(712411);
	}

    }

  if(ThisTask == 0)
    printf("Positions and Velocities fine for (%d,%d,%d) of %d cases on task 0...\n",pok,vok,vsph,NumPart*3);

  MPI_Barrier(MPI_COMM_WORLD);
}


#ifdef GLASS /* This module is used to relax MHD IC to analytic equilibrium stats */

double taper1(double x, double x0, double sharp)
{
  double f;
  f = 1.0 - 1.0/(exp((x-x0)/sharp)+1.);
  return f;
}

double dr_taper1(double x, double x0, double sharp)
{
  double dfdr;
  dfdr = 1.0/(exp((x-x0)/sharp)+1.)/(exp((x-x0)/sharp)+1.)/sharp*exp((x-x0)/sharp);
  return dfdr;
}

double taper2(double x, double rin, double rout)
{
  double f;
  f = sqrt(rin/x) + sqrt(rout/x) - sqrt(rin*rout)/x/x -1;
  return f;
}

double dr_taper2(double x, double rin, double rout)
{
  double dfdr;
  dfdr = -0.5 * (sqrt(rin) + sqrt(rout)) * pow(x,-1.5) + sqrt(rin*rout)/x/x;
  return dfdr;
}

void eq_relax1(void)
{
  int i,k,count=0,count_tot=0; /*counting properly sampled particles */
  double X0,Sharp,V_k,C_s,H0,Aspect0, beta=20.0,Rin=0.2, Rout=10.0, Rho_Power=-3.0, Temp_Power=-1.0,mu=1.0,T0=300.0,in_over_out=0.1;
  double StarMass,rho0, Radius, vt=0, bt=0, vr=0, vtheta=0, sin_v=0, cos_v=0,Z=0,dampfactor=0.9;
  /*  double gizmo2gauss = sqrt(4.*M_PI*All.UnitPressure_in_cgs) / All.UnitMagneticField_in_gauss;*/
  double locrhoratio = 0, glbrhoratio = 0, locvrmax=0,locvzmax=0,glbvrmax=0,glbvzmax=0;
  StarMass =10.0;
  rho0 = 6.14258149044;/* 4.43476867871; */  /*must be calculated numerically using build_ic.py*/
  X0=6.0*Rin;
  V_k = sqrt(10.0/0.2)*sqrt(0.2/Rin*StarMass/10.0);
  C_s =0.41762*sqrt(T0/300.0/mu);
  H0 = C_s / V_k * sqrt((1.0+beta)/beta) * Rin;
  Aspect0 =H0/Rin;
  Sharp = (Rin-X0)/log(1./(1.0-in_over_out*pow((Rout/Rin),Rho_Power))-1.);
  
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  Radius = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1]);
	  sin_v = P[i].Pos[1]/Radius;
	  cos_v = P[i].Pos[0]/Radius;
	  Z = P[i].Pos[2];
	  vr = P[i].Vel[0] * cos_v + P[i].Vel[1] * sin_v;
	  vtheta = - P[i].Vel[0] * sin_v + P[i].Vel[1] * cos_v;
	  
	  if (P[i].Vel[2] < locvzmax) locvzmax = P[i].Vel[2];
	  if (vr < locvrmax) locvrmax = vr;
	  
	  
	  /*tangential vel from analytical model */
	  vt  = V_k * V_k * Rin / Radius  + 2.0 / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) + \
	    (1.0 + beta) / beta * C_s * C_s * Temp_Power * pow((Radius /Rin),Temp_Power) \
	    + (1.0 + beta) / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) * Rho_Power \
	    + V_k * V_k /(Radius /Rin) * (-Radius /sqrt(Radius *Radius +Z *Z )+1.0)*Temp_Power \
	    + (1.0 + beta) / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) * \
	    dr_taper1(Radius,X0,Sharp)*Radius /taper1(Radius,X0,Sharp);

	  vt = sqrt(vt);

	  vr =  vr * dampfactor;                                                             /*relax velocities  */
	  vtheta = vt + (vtheta - vt) * dampfactor;

	  /*	  if(Radius < X0)
	    {
	      vr = 0;
	      vtheta = vt;
	      } */
	  
	  P[i].Vel[2] = P[i].Vel[2] * dampfactor;
	  P[i].Vel[0] = vr * cos_v - vtheta * sin_v;
	  P[i].Vel[1] = vr * sin_v + vtheta * cos_v;
	  /*	  P[i].dp[0] = 0;
	  P[i].dp[1] = 0;
	  P[i].dp[2] = 0;*/
	  
	  bt = rho0 * exp(1./Aspect0/Aspect0 * pow((Radius /Rin),(-Temp_Power-1.))*(1.0/sqrt(1+Z*Z/Radius/Radius)-1)) * \
	    taper1(Radius,X0,Sharp) * pow((Radius /Rin),Rho_Power);	  
	  SphP[i].InternalEnergy = (1.0/(GAMMA-1)+1.0/beta) *C_s * C_s * pow((Radius/Rin), Temp_Power);
	  SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;


	  H0 *= pow((Radius/Rin),(Temp_Power+3.0)/2.0);
	  vr = SphP[i].Density / bt;
	  if  (abs(Z) < 1.5 * H0)
	    {
	      locrhoratio += vr;
	      count += 1;		    
	    }
	  /*
	  bt = sqrt(8.0 * 3.1415926 *34.55277 /beta *bt )*C_s * pow(Radius /Rin,Temp_Power/2.0)/gizmo2gauss;
	  
	  SphP[i].B[2] = 0;       
	  SphP[i].B[0] = -bt * sin_v * P[i].Mass / SphP[i].Density;
	  SphP[i].B[1] = bt * cos_v * P[i].Mass/ SphP[i].Density;	     
	  SphP[i].BPred[2] = 0;
	  SphP[i].BPred[0] = -bt * sin_v * P[i].Mass / SphP[i].Density;
	  SphP[i].BPred[1] = bt * cos_v * P[i].Mass / SphP[i].Density;
	  */
	}
    }
  
  MPI_Allreduce(&locrhoratio, &glbrhoratio, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count, &count_tot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&locvzmax, &glbvzmax, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&locvrmax, &glbvrmax, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  glbrhoratio /= count_tot;
  if (ThisTask == 0)
    {
      printf("vrmax= %f\n",glbvrmax);
      printf("vzmax= %f\n",glbvzmax);
      printf("rhoratio code/eqsolution %f\n", glbrhoratio);
    }
}


void eq_relax2(void)
{
  int i,k,count=0,count_tot=0; /*counting properly sampled particles */
  double X0,Sharp,V_k,C_s,H0,Aspect0, beta=1.0,Rin=0.2, Rout=5.0, Rho_Power=-1.5, Temp_Power=-1.0,mu=1.0,T0=300.0,in_over_out=0.1;
  double StarMass,rho0, Radius, vt=0, bt=0, vr=0, vtheta=0, sin_v=0, cos_v=0,Z=0,dampfactor=0.99;
  /*  double gizmo2gauss = sqrt(4.*M_PI*All.UnitPressure_in_cgs) / All.UnitMagneticField_in_gauss;*/
  double locrhoratio = 0, glbrhoratio = 0, locvrmax=0,locvzmax=0,glbvrmax=0,glbvzmax=0;
  StarMass =10.0;
  rho0 =  0.441918443832;   /*must be calculated numerically using build_ic.py*/
  X0=3.0*Rin;
  V_k = sqrt(10.0/0.2)*sqrt(0.2/Rin*StarMass/10.0);
  C_s =0.41762*sqrt(T0/300.0/mu);
  H0 = C_s / V_k * sqrt((1.0+beta)/beta) * Rin;
  Aspect0 =H0/Rin;
  Sharp = (Rin-X0)/log(1./(1.0-in_over_out*pow((Rout/Rin),Rho_Power))-1.);
  
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  Radius = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1]);
	  sin_v = P[i].Pos[1]/Radius;
	  cos_v = P[i].Pos[0]/Radius;
	  Z = P[i].Pos[2];
	  vr = P[i].Vel[0] * cos_v + P[i].Vel[1] * sin_v;
	  vtheta = - P[i].Vel[0] * sin_v + P[i].Vel[1] * cos_v;
	  
	  if (P[i].Vel[2] < locvzmax) locvzmax = P[i].Vel[2];
	  if (vr < locvrmax) locvrmax = vr;
	  
	  
	  H0 *= pow((Radius/Rin),(Temp_Power+3.0)/2.0);
	  /*tangential vel from analytical model */
	  if (( Radius > 1.001* Rin) && ( Radius < 0.95 * Rout)){
	    vt  = V_k * V_k * Rin / Radius  + 2.0 / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) + \
	      (1.0 + beta) / beta * C_s * C_s * Temp_Power * pow((Radius /Rin),Temp_Power) \
	      + (1.0 + beta) / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) * Rho_Power \
	      + V_k * V_k /(Radius /Rin) * (-Radius /sqrt(Radius *Radius +Z *Z )+1.0)*Temp_Power \
	      + (1.0 + beta) / beta * C_s * C_s * pow((Radius /Rin),Temp_Power) * \
	      dr_taper2(Radius,Rin,Rout)*Radius /taper2(Radius,Rin,Rout);

	    vt = sqrt(vt);
	  }
	  
	  vr =  vr * dampfactor;                                                             /*relax velocities  */
	  vtheta = vt + (vtheta - vt) * dampfactor;

	  P[i].Vel[2] = P[i].Vel[2] * dampfactor;
	  P[i].Vel[0] = vr * cos_v - vtheta * sin_v;
	  P[i].Vel[1] = vr * sin_v + vtheta * cos_v;
	  /*	  P[i].dp[0] = 0;
	  P[i].dp[1] = 0;
	  P[i].dp[2] = 0;*/

	  
	  bt = rho0 * exp(1./Aspect0/Aspect0 * pow((Radius /Rin),(-Temp_Power-1.))*(1.0/sqrt(1+Z*Z/Radius/Radius)-1)) * \
	    taper2(Radius,Rin,Rout) * pow((Radius /Rin),Rho_Power);	  
	  SphP[i].InternalEnergy = (1.0/(GAMMA-1)+1.0/beta) *C_s * C_s * pow((Radius/Rin), Temp_Power);
	  SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;



	  
	  vr = SphP[i].Density / bt;
	  if  (abs(Z) < 1.5 * H0)
	    {
	      locrhoratio += vr;
	      count += 1;		    
	    }
	  /*
	  bt = sqrt(8.0 * 3.1415926 *34.55277 /beta *bt )*C_s * pow(Radius /Rin,Temp_Power/2.0)/gizmo2gauss;
	  
	  SphP[i].B[2] = 0;       
	  SphP[i].B[0] = -bt * sin_v * P[i].Mass / SphP[i].Density;
	  SphP[i].B[1] = bt * cos_v * P[i].Mass/ SphP[i].Density;	     
	  SphP[i].BPred[2] = 0;
	  SphP[i].BPred[0] = -bt * sin_v * P[i].Mass / SphP[i].Density;
	  SphP[i].BPred[1] = bt * cos_v * P[i].Mass / SphP[i].Density;
	  */
	}
    }
  
  MPI_Allreduce(&locrhoratio, &glbrhoratio, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count, &count_tot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&locvzmax, &glbvzmax, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&locvrmax, &glbvrmax, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  glbrhoratio /= count_tot;
  if (ThisTask == 0)
    {
      printf("vrmax= %f\n",glbvrmax);
      printf("vzmax= %f\n",glbvzmax);
      printf("rhoratio code/eqsolution %f\n", glbrhoratio);
    }
}
  

#endif

#ifdef MOONRELAX
/* Artificial relaxation for settling WoMa planet realizations, ported from
   SPH-EXA's relaxation scheme (sphexa docs/relaxation.md) and adapted to the
   KDK kick. Called from do_the_kick for active gas particles only.

   Linear drag dv/dt = -v/RelaxTimescale is applied while Time < RelaxUntil.
   During the spherical stage (Time < SphericalRelaxUntil) the particle
   velocity and the applied momentum kick are projected onto the radius vector
   about the coordinate origin, i.e. the particle is kinematically constrained
   to its radial ray. Over [SphericalRelaxUntil, SphericalRelaxUntil +
   SphericalRelaxReleaseDuration] the tangential part of the kick is restored
   with the cubic smoothstep w = 3s^2 - 2s^3 (zero slope at both ends);
   velocity is NOT reprojected during the release window, and the velocity is
   never multiplied by w -- that would add a second, unintended damping. After
   the release window, ordinary 3D damping continues until RelaxUntil.

   The drag removes kinetic energy deliberately (it is not thermalized). The
   projection origin is fixed at (0,0,0): single, non-rotating, origin-centered
   planets only -- never enable the spherical stage for impact runs.

   dp is the momentum kick mass*acceleration*dt about to be applied; dt is the
   particle's own (hydro) kick interval. */
void moonrelax_modify_kick(int i, double dp[3], double mass, double dt)
{
  int j;
  double w = 1.0; /* tangential acceleration weight: 1 = full 3D dynamics */

  if(All.RelaxTimescale <= 0.0)
    return;

  if(All.Time < All.SphericalRelaxUntil)
    w = 0.0;
  else if(All.SphericalRelaxReleaseDuration > 0.0 &&
          All.Time < All.SphericalRelaxUntil + All.SphericalRelaxReleaseDuration)
    {
      double s = (All.Time - All.SphericalRelaxUntil) / All.SphericalRelaxReleaseDuration;
      w = s * s * (3.0 - 2.0 * s);
    }

  if(w < 1.0)
    {
      double r2 = P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2];
      if(r2 > 0.0)
        {
          if(w == 0.0)
            {
              /* spherical stage: keep only the radial velocity component */
              double vrad = (P[i].Vel[0]*P[i].Pos[0] + P[i].Vel[1]*P[i].Pos[1] + P[i].Vel[2]*P[i].Pos[2]) / r2;
              for(j = 0; j < 3; j++)
                {
                  P[i].Vel[j] = vrad * P[i].Pos[j];
                  SphP[i].VelPred[j] = P[i].Vel[j];
                }
            }
          /* project the applied kick: keep its radial part, scale the tangential part by w */
          double dp_rad = (dp[0]*P[i].Pos[0] + dp[1]*P[i].Pos[1] + dp[2]*P[i].Pos[2]) / r2;
          for(j = 0; j < 3; j++)
            dp[j] = dp_rad * P[i].Pos[j] + w * (dp[j] - dp_rad * P[i].Pos[j]);
        }
    }

  if(All.Time < All.RelaxUntil)
    for(j = 0; j < 3; j++)
      dp[j] -= mass * P[i].Vel[j] / All.RelaxTimescale * dt;
}

#endif // MOONRELAX




#ifdef DISKIC
void eq_relax5(void)
{
  int i,j;
  double dampfactor=0.999,dampfactor1=0.999;
  double radius=0, vr=0,zzz=0,vphi=0,sigma_in=2.65e-4;
  double nx,ny,nz;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
        {
          radius=sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1]);
          zzz=P[i].Pos[2];
          if(radius>0.5)
          {
	    // SphP[i].InternalEnergy=0.6*pow(radius,-0.59)*0.7*0.7;
	    //	    SphP[i].InternalEnergy=2.5*pow(radius,-1);
	    //SphP[i].InternalEnergy=3.75e-3*pow(radius,-1)*(1+1/radius/radius);
	    	    SphP[i].InternalEnergy=0.4*pow(radius,-1)*(1+1./radius/radius);
	    //SphP[i].InternalEnergy=6e-4*pow(radius,-1)*(1+1./radius/radius);
          }
          /*if(P[i].imat==1)
          {
            P[i].Mass =(1. - exp(-0.5*pow((All.Time/10),0.5)))*4*0.000217;
          }
          else
          {
            P[i].Mass =(1. - exp(-0.5*pow((All.Time/10),0.5)))*4*0.000107;
            }*/
          


          nx=P[i].Pos[0]/radius;
          ny=P[i].Pos[1]/radius;
          vphi=-P[i].Vel[0]*ny + P[i].Vel[1]*nx;          
          vr=P[i].Vel[0]*nx + P[i].Vel[1]*ny; //
          //          SphP[i].InternalEnergy=0.02*0.02*pow(radius,-1.5)/1.001/0.001;
          
          //          vr=1.;
          if(radius>0)
	    {
	      P[i].Vel[0] =0.- vphi*ny + vr  * nx *dampfactor1; 
	      P[i].Vel[1] =vphi*nx + vr  * ny *dampfactor1; 
	      P[i].Vel[2] *=dampfactor;
            }
	  
        }
    }
}
#endif
#ifdef NEWTEST
void eq_relax3(void)
{
  int i,count;
  double vt=0,dampfactor=0.99,locrhoratio=0,glbrhoratio=0;
  count = 0;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  vt = 0. - SHEARING_BOX_Q * SHEARING_BOX_OMEGA_BOX_CENTER *(P[i].Pos[0] - boxHalf_X);
	  P[i].Vel[2] = P[i].Vel[2] * dampfactor;
	  P[i].Vel[0] = P[i].Vel[0] * dampfactor;
	  P[i].Vel[1] = vt + (P[i].Vel[1] - vt) * dampfactor;

	  SphP[i].InternalEnergy = 1.0/(GAMMA-1)/GAMMA;
	  SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;          
	  if(fabs(P[i].Pos[2]-boxHalf_Z)<0.3)
	    {
	      locrhoratio += SphP[i].Density/exp(- 0.5*(P[i].Pos[2]-boxHalf_Z)*(P[i].Pos[2]-boxHalf_Z));
	      count += 1;
	    }
	}
    }
  locrhoratio /= count;
  MPI_Allreduce(&locrhoratio, &glbrhoratio, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (ThisTask == 0)
    {
      printf("rhoratio code/eqsolution %f\n", glbrhoratio);
    }
}

#endif
