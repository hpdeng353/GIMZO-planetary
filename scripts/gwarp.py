import h5py
import pynbody
import numpy as np
import matplotlib.pyplot as plt
import diskpy
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta


#s = h5py.File('clump.hdf5','r+')
s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk1/init/grvdisk1-ics.hdf5','r+')

#rho = s["PartType0"]["Density"][:]

width=300
center=15.
aspect=0.05
amp=3
vel=s["PartType0"]["Velocities"][:]
del s["PartType0"]["Velocities"]

pos=s["PartType0"]["Coordinates"]
del s["PartType0"]["Coordinates"]

#warning  the midplane rad rad1
zz=np.copy(pos[:,2])
xx=np.copy(pos[:,0])
yy=np.copy(pos[:,1])
vx=np.copy(vel[:,0])
vz=np.copy(vel[:,2])
vy=np.copy(vel[:,1])



rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])
rr=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
tempx=4./3./aspect*np.power(rad, 1.5)
center=4./3./aspect*np.power(center, 1.5)
sint=rad/rr
cost=zz/rr
sinp=yy/rad
cosp=xx/rad
omega=np.power(rr,-1.5)
vr=-rr*omega*sint*cost*cosp*np.exp(-(tempx-center)**2/2./width/width)*amp
vp=0.5*rr*omega*sint*cost*sinp*np.exp(-(tempx-center)**2/2./width/width)*amp
ut=-aspect*rr*cosp*np.exp(-(tempx-center)**2/2./width/width)/width/width*(center-tempx)*amp


vx+=vr*sint*cosp + ut*cost*cosp  - vp*sinp
vy+=vr*sint*sinp + ut*cost*sinp  + vp*cosp
vz+=vr*cost      - ut*sint       + 0.


lx=aspect*np.exp(-(tempx-center)**2/2./width/width)*amp
lz=np.sqrt(1-lx*lx)



pos[:,0]=lz*xx+lx*zz
pos[:,2]=-lx*xx+lz*zz
vel[:,0]=lz*vx+lx*vz
vel[:,2]=-lx*vx+lz*vz
vel[:,1]=vy





s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType0/Velocities",data=vel)




s.flush()
s.close()
