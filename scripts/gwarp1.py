import h5py
#import pynbody
import numpy as np
#import matplotlib.pyplot as plt
import pynbody
#import diskpy
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta


#s = h5py.File('clump.hdf5','r+')


s = h5py.File('/scratch/e1000/hpdeng/gizmout/grvdisk4/init/grvdisk4-ics.hdf5','r+')

vel=s["PartType0"]["Velocities"][:]
del s["PartType0"]["Velocities"]

pos=s["PartType0"]["Coordinates"]
del s["PartType0"]["Coordinates"]


#mass= s["PartType2"]["Masses"][:]
#del s["PartType2"]["Masses"]


#s1=pynbody.load('/scratch/e1000/hpdeng/gizmout/grvdisk4/init/grvdisk4-ics1.hdf5')
#pynbody.analysis.halo.center(s1,'hyb')



#mass[0]=0
#warning  the midplane rad rad1
#pos=s1.g['pos']
#vel=s1.g['vel']

xx=pos[:,0]
yy=pos[:,1]
zz=pos[:,2]
vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]


rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])
rr=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])

sint=rad/rr
cost=zz/rr
sinp=yy/rad
cosp=xx/rad
omega=np.power(rr,-1.5)

'''
amp=0.1
lx=amp*(1+np.tanh((rad-12)/5))/2
lz=np.sqrt(1-lx*lx)

'''
amp=0.25
lx=amp*np.exp(-1/0.02/rad/rad)
lz=np.sqrt(1-lx*lx)


vphi=-sinp*vx+cosp*vy
vr=cosp*vx+sinp*vy
vr+=np.sqrt(1./rad)*zz/(0.02*rad)*lx*sinp
vphi+=0.5*np.sqrt(1./rad)*zz/(0.02*rad)*lx*cosp
vx=-vphi*sinp+vr*cosp
vy=vphi*cosp+vr*sinp


pos[:,0]=lz*xx+lx*zz
pos[:,2]=-lx*xx+lz*zz
vel[:,0]=lz*vx+lx*vz
vel[:,2]=-lx*vx+lz*vz
vel[:,1]=vy



#vel1=s1.s['vel']
#pos1=s1.s['pos']


s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType0/Velocities",data=vel)
#s.create_dataset("PartType2/Coordinates",data=pos1)
#s.create_dataset("PartType2/Velocities",data=vel1)

#s.create_dataset("PartType2/Masses",data=mass)                

s.flush()
s.close()
