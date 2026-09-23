import h5py
import pynbody
import numpy as np
import matplotlib.pyplot as plt
#import diskpy
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta
ample=0.

#s = h5py.File('clump.hdf5','r+')
s = h5py.File('/scratch/e1000/hpdeng/gizmout/sgwarp1/init/sgwarp1-ics.hdf5','r+')

#temp = s["PartType0"]["Temperature"][:]
#rho = s["PartType0"]["Density"][:]
#del s["PartType0"]["Temperature"]
mass= s["PartType0"]["Masses"][:]
del s["PartType0"]["Masses"]
mass*=0.5
#mass1=0.
#mag= s["PartType0"]["MagneticField"][:]
#del s["PartType0"]["MagneticField"]
#mat = s["PartType0"]["Materials"][:]
#del s["PartType0"]["Materials"]
#vel=s["PartType0"]["Velocities"][:]
#del s["PartType0"]["Velocities"]
#vel1=s["PartType2"]["Velocities"][:]
#del s["PartType2"]["Velocities"]
#ids=np.arange(1,len(rho)+1)
#ids=s["PartType2"]["ParticleIDs"][:]
#del s["PartType2"]["ParticleIDs"]
#mag=vel
#ids[0]=50000000

u=s["PartType0"]["InternalEnergy"][:]
#del s["PartType0"]["InternalEnergy"]

#u[np.where(rho>0.1)]*=20.
#entr=s["PartType0"]["Entropy"][:]
#mass = s["PartType2"]["Masses"][:]
#del s["PartType2"]["Masses"]
#mass[0]=0
#massmin=np.min(mass)

#vfrac=s["PartType0"]["Vfrac"][:]

pos=s["PartType0"]["Coordinates"][:]
#del s["PartType0"]["Coordinates"]


#xx=pos[:,0]
#yy=pos[:,1]
radius=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])

#mag[:,:]*=600
#u*=666.6667
#=6e-4/radius*(1+3./radius/radius)





#s.create_dataset("PartType0/Coordinates",data=pos)
#s.create_dataset("PartType2/Coordinates",data=pos1)
#
#s.create_dataset("PartType0/Materials",data=mat)
#s.create_dataset("PartType0/InternalEnergy",data=u)
#s.create_dataset("PartType0/Temperature",data=temp)


s.create_dataset("PartType0/Masses",data=mass)
#s.create_dataset("PartType0/MagneticField",data=mag)
#s.create_dataset("PartType0/Velocities",data=vel)
#s.create_dataset("PartType2/Velocities",data=vel1)
#s.create_dataset("PartType2/ParticleIDs",data=ids)


s.flush()
s.close()
"""            
xx=radius[np.where((zz<5)&(zz>-5)&(mat==1))]

yy=entr[np.where((zz<5)&(zz>-5)&(mat==1))]
yy*=1e6
plt.scatter(xx,yy, s=1,c="#fb0a2a")
plt.show()
"""
