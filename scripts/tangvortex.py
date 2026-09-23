################################################################################
###### This is an example script to generate HDF5-format ICs for GIZMO
######  The specific example below is obviously arbitrary, but could be generalized
######  to whatever IC you need. 
################################################################################
################################################################################

## load libraries we will use 
import numpy as np
import h5py as h5py
import random
# the main routine. this specific example builds an N-dimensional box of gas plus 
#   a collisionless particle species, with a specified mass ratio. the initial 
#   gas particles are distributed in a uniform lattice; the initial collisionless 
#   particles laid down randomly according to a uniform probability distribution 
#   with a specified random velocity dispersion
#

def makeIC_box():
    DIMS=3; # number of dimensions 
    N_1D=128; # 1D particle number (so total particle number is N_1D^DIMS)
    fname='box2d.hdf5'; # output filename 

    Lbox = 1.0 # box side length
    Lx=1.
    Ly=1.
    Lz=0.125
    rx=128
    ry=128
    rz=16
    Omega=1.
    Cs=1.
    rho_desired = 1.0 # box average initial gas density
    P_desired = rho_desired * Cs * Cs # initial gas pressure
    gamma_eos = 5./3  # polytropic index of ideal equation of state the run will assume
    Kz=np.pi*2.
    theta=np.pi/4.
    sin=np.sin(theta)
    cos=np.cos(theta)
    VaKz=np.sqrt(3.*sin*sin*(1.-0.75*cos*cos))
    print VaKz
    v0=1.5/Kz*sin*sin
    ample=0.0
    ample1=0.05
    # first we set up the gas properties (particle type 0)
    
    # make a regular 1D grid for particle locations (with N_1D elements and unit length)
    x0=np.arange(0.,1.,1./rx); x0+=0.5*(1.-x0[-1]); # shift it so the lattice is distributed in the center of box instead of occupying oneside.
    y0=np.arange(0.,1.,1./ry); y0+=0.5*(1.-y0[-1]);
    z0=np.arange(0.,1.,1./rz); z0+=0.5*(1.-z0[-1]);
    # now extend that to a full lattice in DIMS dimensions
    if(DIMS==3):
        xv_g, yv_g, zv_g = np.meshgrid(x0,y0,z0, sparse=False, indexing='xy')
    if(DIMS==2):
        xv_g, yv_g = np.meshgrid(x0,x0, sparse=False, indexing='xy'); zv_g = 0.0*xv_g
    if(DIMS==1):
        xv_g=x0; yv_g = 0.0*xv_g; zv_g = 0.0*xv_g; 
    # the gas particle number is the lattice size: this should be the gas particle number
    Ngas = xv_g.size
    # flatten the vectors (since our ICs should be in vector, not matrix format): just want a
    #  simple list of the x,y,z positions here. Here we multiply the desired box size in
    xv_g=xv_g.flatten()*Lx

    yv_g=yv_g.flatten()*Ly
    zv_g=zv_g.flatten()*Lz
    # set the initial velocity in x/y/z directions (here zero)
    vx_g=0.*xv_g; vy_g=-1.5*Omega*(xv_g-0.5*Lx); vz_g=0.*xv_g;

    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*xv_g; by_g=0.*xv_g; bz_g=0.*xv_g;
    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    #   since their space-density is uniform this gives a uniform density, at the desired value
    mv_g=rho_desired/((1.*Ngas)/(Lx*Ly*Lz)) + 0.*xv_g

    vx_g=-np.sin(2*np.pi*yv_g); vy_g=np.sin(2*np.pi*xv_g); vz_g=0.*xv_g #-1.5*(xv_g-0.5); # r z phi
    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=-np.sin(2*np.pi*yv_g); by_g=np.sin(4*np.pi*xv_g); bz_g=0.*xv_g;
    bx_g/=np.sqrt(4*np.pi)
    by_g/=np.sqrt(4*np.pi)
    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    uv_g=P_desired/((gamma_eos-1.)*rho_desired) + 0.*xv_g
    # set the gas IDs: here a simple integer list
    print uv_g[0]
    id_g=np.arange(1,Ngas+1)

    B0=np.sqrt(4*np.pi)*VaKz/Kz
    B0=np.sqrt(8.*np.pi/400)*Cs
    beta=B0*B0/8./np.pi
    beta=1./beta

    Qz=np.pi*2*VaKz/Kz*rz
    print "beta=",beta,"Qz=",Qz,"v0=",v0
    
    bz_g=B0+0*xv_g
#    bz_g[np.where((xv_g<0.75)&(xv_g>0.25)&(yv_g<0.75)&(yv_g>0.25))]=B0
    by_g=0. - ample*B0*np.cos(Kz*zv_g)*cos # there should be a minus sign due to sin(theta-pi/2) but due to cord rotation by is z 
    bx_g= ample*B0*np.cos(Kz*zv_g)*sin
    vy_g += ample*v0*np.sin(Kz*zv_g)*sin
    vx_g = ample*v0*np.sin(Kz*zv_g)*cos


    bx_g += ample1*np.random.uniform(-1.,1.,Ngas)*B0
    by_g += ample1*np.random.uniform(-1.,1.,Ngas)*B0
    bz_g += ample1*np.random.uniform(-1.,1.,Ngas)*B0
    vx_g += ample1*np.random.uniform(-1.,1.,Ngas)*v0
    vy_g += ample1*np.random.uniform(-1.,1.,Ngas)*v0
    vz_g += ample1*np.random.uniform(-1.,1.,Ngas)*v0

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");
    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it 
                                        # allows for physics which can change particle masses 
    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 2.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    ## ok, that ends the block of 'useless' parameters
    
    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs
    
    # start with particle type zero. first (assuming we have any gas particles) create the group 
    p = file.create_group("PartType0")
    # now combine the xyz positions into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=xv_g; q[:,1]=yv_g; q[:,2]=zv_g;
    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)
    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=vx_g; q[:,1]=vy_g; q[:,2]=vz_g;
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q)
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)
    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=uv_g)
    # combine the xyz magnetic fields into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=bx_g; q[:,1]=by_g; q[:,2]=bz_g;
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)


    file.close()
    # all done!
makeIC_box()
