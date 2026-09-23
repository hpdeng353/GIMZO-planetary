import h5py
import numpy as np


beta=25
gamma_eos=5.0/3.
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta



s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvmhd2/init/grvmhd2-ics.hdf5','r+')

#mag = s["PartType0"]["Velocities"][:]
pos = s["PartType0"]["Coordinates"][:]
u = s["PartType0"]["InternalEnergy"][:]
del s["PartType0"]["InternalEnergy"]
rho = s["PartType0"]["Density"][:]
#id1= s["PartType1"]["ParticleIDs"][:]
#id1[0]=50000000
mass=s["PartType0"]["Masses"][:]
#del s["PartType0"]["Masses"]

mag=s["PartType0"]["MagneticField"]
#del s["PartType0"]["InternalEnergy"]
#del s["PartType1"]["ParticleIDs"]

me=(mag[:,0]*mag[:,0]+mag[:,1]*mag[:,1]+mag[:,2]*mag[:,2])/8/3.14

p_fiducial=rho*u*(gamma_eos-1)*5.27*1e6#np.average(rho1*cs21)*5.27*1e6
beta=p_fiducial/me
u=u*(1+1./beta)
#0.002 pure vertical

#u*=66.
#s.create_dataset("PartType0/MagneticField",data=mag)
#s.create_dataset("PartType0/Masses",data=mass)
s.create_dataset("PartType0/InternalEnergy",data=u)
#s.create_dataset("PartType1/ParticleIDs",data=id1)
s.flush()
s.close()
            

