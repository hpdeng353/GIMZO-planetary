# coding: utf-8
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pynbody




s=pynbody.load("test.hdf5")
pynbody.analysis.halo.center(s,mode="hyb")
p=pynbody.analysis.profile.Profile(s.g,min=4,max=38,ndim=3,nbins=200)
bins=p["rbins"]
vphi=p["vphi"][:]
pos=np.copy(s.g["pos"])
rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
delta=(p["rbins"][1]-p["rbins"][0])/2.
posx=pos[:,0]
posy=pos[:,1]
vx=np.copy(s.g["vel"][:,0])
vy=np.copy(s.g["vel"][:,1])

for i in range(0,len(p["rbins"])):
    vy[np.where((rad>bins[i]-delta)&(rad<bins[i]+delta))]-=posx[np.where((rad>bins[i]-delta)&(rad<bins[i]+delta))]/bins[i]*vphi[i]
    vx[np.where((rad>bins[i]-delta)&(rad<bins[i]+delta))]+=posy[np.where((rad>bins[i]-delta)&(rad<bins[i]+delta))]/bins[i]*vphi[i]
    


s1 = h5py.File('test1.hdf5','r+')
vel=s1["PartType0"]["Velocities"][:]
del s1["PartType0"]["Velocities"]

vel[:,1]=vy
vel[:,0]=vx

s1.create_dataset("PartType0/Velocities",data=vel)    
s1.flush()
s1.close()