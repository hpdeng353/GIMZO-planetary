# coding: utf-8
import woma
import matplotlib.pyplot as plt
import h5py as h5py
import numpy as np

R_earth = 6.371e6   # m
M_earth = 5.9724e24  # kg

def plot_spherical_profiles(planet):    
    fig, ax = plt.subplots(2, 2, figsize=(8,8))
    
    ax[0, 0].plot(planet.A1_r / R_earth, planet.A1_rho)
    ax[0, 0].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[0, 0].set_ylabel(r"Density, $\rho$ [kg m$^{-3}$]")
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlim(0, None)
    
    ax[1, 0].plot(planet.A1_r / R_earth, planet.A1_m_enc / M_earth)
    ax[1, 0].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[1, 0].set_ylabel(r"Enclosed Mass, $M_{<r}$ $[M_\oplus]$")
    ax[1, 0].set_xlim(0, None)
    ax[1, 0].set_ylim(0, None)
    
    ax[0, 1].plot(planet.A1_r / R_earth, planet.A1_P)
    ax[0, 1].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[0, 1].set_ylabel(r"Pressure, $P$ [Pa]")
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_xlim(0, None)
    
    ax[1, 1].plot(planet.A1_r / R_earth, planet.A1_T)
    ax[1, 1].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[1, 1].set_ylabel(r"Temperature, $T$ [K]")
    ax[1, 1].set_xlim(0, None)
    ax[1, 1].set_ylim(0, None)
    
    plt.tight_layout()
    plt.show()
    

planet = woma.Planet(
    name            = "Earth",
    A1_mat_layer    = ["ANEOS_Fe85Si15", "ANEOS_forsterite"],
    A1_T_rho_type   = ["adiabatic", "adiabatic"],
    A1_M_layer      = [0.3 * M_earth, 0.7 * M_earth],
    P_s             = 1e5,
    T_s             = 1000,
)


planet.gen_prof_L2_find_R_R1_given_M1_M2(R_min=0.9 * R_earth, R_max=1.05 * R_earth)

#plot_spherical_profiles(planet)

coredge=planet.A1_idx_layer[0]
cmbratio=1.0*planet.A1_rho[coredge]/planet.A1_rho[coredge+1]

# here we merge models of two different mass resolution so that the core particle mass is cmbratio times larger than the mantle
#and the particles have no abrupt change in the smoothing length.
lresn=1e6
hresn=int(cmbratio*lresn)
hres = woma.ParticlePlanet(planet, hresn, verbosity=0)
lres = woma.ParticlePlanet(planet, lresn, verbosity=0)

mc=lres.A1_m[np.where(lres.A1_mat_id==402)]
mm=hres.A1_m[np.where(hres.A1_mat_id==400)]

m=np.concatenate((mc,mm))
m/=9.56072e22

uc=lres.A1_u[np.where(lres.A1_mat_id==402)]
um=hres.A1_u[np.where(hres.A1_mat_id==400)]
u=np.concatenate((uc,um))
u/=1e6


idc=np.arange(0,len(mc))
idm=np.arange(0,len(mm))+int(1e8)
ids=np.concatenate((idc,idm))

temp=np.zeros(len(ids))
entr=np.zeros(len(ids))

imat=np.zeros(len(ids))
imat[:]=0
imat[np.where(ids<1e8)]=1

posc=lres.A2_pos[np.where(lres.A1_mat_id==402)]
posm=hres.A2_pos[np.where(hres.A1_mat_id==400)]
pos=np.concatenate((posc,posm))
pos/= 6.37869e+06

velc=lres.A2_vel[np.where(lres.A1_mat_id==402)]
velm=hres.A2_vel[np.where(hres.A1_mat_id==400)]
vel=np.concatenate((velc,velm))

print (vel.shape,pos.shape)

fname='Earth.hdf5'

Ngas=len(ids)
file = h5py.File(fname,'w') 
npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others
        
h = file.create_group("Header");
h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
h.attrs['NumPart_Total'] = npart; # npart set as above
h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
h.attrs['Time'] = 0.0;  # initial time
h.attrs['Redshift'] = 0.0; # initial redshift
h.attrs['BoxSize'] = 1.0; # box size
h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
h.attrs['Omega0'] = 0.0; # z=0 Omega_matter
h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
h.attrs['HubbleParam'] = 0.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    ## ok, that ends the block of 'useless' parameters
    
    # Now, the actual data!


p = file.create_group("PartType0")

p.create_dataset("Coordinates",data=pos)

p.create_dataset("Velocities",data=vel)

p.create_dataset("ParticleIDs",data=ids)

p.create_dataset("Materials",data=imat)

p.create_dataset("Masses",data=m)

p.create_dataset("InternalEnergy",data=u)

p.create_dataset("Temperature",data=temp)

p.create_dataset("Entropy",data=entr)


file.close()
    # no P
