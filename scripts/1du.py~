import h5py
import numpy as np
import matplotlib.pyplot as plt
Lx=1.414
Lz=24.
P0=4*np.pi
R0=Lx/4.
betay=25
betap=1600
gamma_eos=5.0/3.
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta
ample=0.


s = h5py.File('/scratch/snx3000/hpdeng/gizmout/moon1/init/moon1-ics.hdf5','r+')

mat = s["PartType0"]["Materials"][:]
#del s["PartType0"]["Materials"]

u=s["PartType0"]["InternalEnergy"][:]
del s["PartType0"]["InternalEnergy"]

pos= s["PartType0"]["Coordinates"][:]

#mass=s["PartType0"]["Masses"][:]
#del s["PartType0"]["Masses"]
#totalmass=np.sum(mass)
#vel[:,1]=-1.5*(pos[:,0]-0.5)
#vel[:,0]=0
#vel[:,2]=0


radius=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])


rbins=np.arange(0,0.85,0.02)
uprofile=np.zeros(len(rbins)-1)
for i in range(0,len(rbins)-1):
    uprofile[i]=np.average(u[np.where((rbins[i]<radius)&(rbins[i+1]>radius))])

    
for i in range(0,len(rbins)-2):
    u[np.where((rbins[i]<radius)&(rbins[i+1]>radius))]=uprofile[i]


#u[np.where(u<1)]=2.4
#mass[np.where(radius<2)]=0
#mass[np.where(radius>6)]=0
#s.create_dataset("PartType0/Materials",data=mat)
s.create_dataset("PartType0/InternalEnergy",data=u)
#s.create_dataset("PartType0/Velocities",data=vel)
#s.create_dataset("PartType0/Masses",data=mass)

s.flush()
s.close()
          
