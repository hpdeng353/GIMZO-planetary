import h5py
import numpy as np


#oringinal 1.12 Lem

s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk1/init/grvdisk1-ics.hdf5','r+')

rho=s["PartType0"]["Density"][:]
pos =s["PartType0"]["Coordinates"][:]
pos1=s["PartType1"]["Coordinates"][:]
vel=s["PartType0"]["Velocities"][:]
#vel*=0.71/0.82*1.07
#del s["PartType0"]["Velocities"]
del s["PartType0"]["Coordinates"]
del s["PartType1"]["Coordinates"]


pos[:,0]-=pos1[0,0]
pos[:,1]-=pos1[0,1]
pos[:,2]-=pos1[0,2]
pos1[:,:]=0

#s.create_dataset("PartType0/Velocities",data=vel)
s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType1/Coordinates",data=pos1)
s.flush()
s.close()
            

