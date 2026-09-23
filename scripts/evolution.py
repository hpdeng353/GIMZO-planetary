# coding: utf-8
import numpy as np
import h5py as h5py
import math as m
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pyplot as plt


alpha=[]
beta=[]
time=[]
startnum=0
step=1
timestep=0.1
numbsnaps=50
fname_base="/scratch/snx3000/hpdeng/gizmout/TSPHR0.2T2NX64/snapshot_"

fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s = h5py.File(fname,"r+")
        density = s["PartType0"]["Density"][:]
        pos = s["PartType0"]["Coordinates"][:]
        mass = s["PartType0"]["Masses"][:]
    
        pos[:,0]-=0.5
        pos[:,1]-=0.5
        Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1]
        Radius = np.sqrt(Radius)
#        vortm= np.sqrt(
        vort =vort[:,0]*vort[:,0]+vort[:,1]*vort[:,1]+vort[:,2]*vort[:,2])
#        vort1=vortm[np.where(Radius>8.)]
#        rho=density[np.where(Radius<0.2)]
        masses=mass[np.where(Radius<0.2)]
        time.append(i*timestep)
        alpha.append(np.sum(masses))
 #       beta.append(np.max(vort1))
        print i
