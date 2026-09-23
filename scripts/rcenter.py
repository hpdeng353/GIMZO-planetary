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
s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk1/init/grvdisk1-ics.hdf5','r+')

npart=s["Header"].attrs["NumPart_ThisFile"]
Npart=s["Header"].attrs["NumPart_Total"]
Ngas=npart[0]

npart = npart + [0,1,0,0,0,0]
Npart = Npart + [0,1,0,0,0,0]

s["Header"].attrs.modify("NumPart_ThisFile",npart)
s["Header"].attrs.modify("NumPart_Total",Npart)


vel=s["PartType0"]["Velocities"][:]
del s["PartType0"]["Velocities"]
pos=s["PartType0"]["Coordinates"][:]
del s["PartType0"]["Coordinates"]
"""
vel1=s["PartType1"]["Velocities"][:]
del s["PartType1"]["Velocities"]
pos1=s["PartType1"]["Coordinates"][:]
del s["PartType1"]["Coordinates"]
"""

vel2=s["PartType2"]["Velocities"][:]
del s["PartType2"]["Velocities"]
pos2=s["PartType2"]["Coordinates"][:]
del s["PartType2"]["Coordinates"]

mass=s["PartType2"]["Masses"][:]
del s["PartType2"]["Masses"]
mass[:]=0

pos-=pos2[0]
vel-=vel2[0]
"""
pos1[0]-=pos2[0]
vel1[0]-=vel2[0]
"""
pos2[0]-=pos2[0]
vel2[0]-=vel2[0]


ids=[50000001,]
mass1=[0.33,]

s.create_dataset("PartType1/Masses",data=mass1)
s.create_dataset("PartType1/ParticleIDs",data=ids)



s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType0/Velocities",data=vel)


s.create_dataset("PartType1/Coordinates",data=pos2)
s.create_dataset("PartType1/Velocities",data=vel2)

s.create_dataset("PartType2/Coordinates",data=pos2)
s.create_dataset("PartType2/Velocities",data=vel2)
s.create_dataset("PartType2/Masses",data=mass)

s.flush()
s.close()