# coding: utf-8
import numpy as np
import h5py as h5py
import math as m
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pyplot as plt

me=[]
rhov2=[]
rhou=[]
mstress=[]
rstress=[]
startnum=0
time=[]
step=1
timestep=0.5/6.28
numbsnaps=170
fname_base="snapshot_"
#fname_base ="/home/ics/hpdeng/snapshot_"                                                                                                                                                                           
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s = h5py.File(fname,"r+")
        density = s["PartType0"]["Density"][:]
        pos = s["PartType0"]["Coordinates"][:]
        mag = s["PartType0"]["MagneticField"][:]
        vel = s["PartType0"]["Velocities"][:]
        Mass = s["PartType0"]["Masses"][:]
        u = s["PartType0"]["InternalEnergy"][:]
        Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1]
        Radius = np.sqrt(Radius)
        velx=vel[:,0]
        vely=vel[:,1]
        velz=vel[:,2]
  
        Bx=mag[:,0]
        By=mag[:,1]
        
        x=pos[:,0]-0.5
        
        vely=vely-1.5*x

        ME = mag[:,0]*mag[:,0]+mag[:,1]*mag[:,1]+mag[:,2]*mag[:,2]
        ME/=8.
        ME/=np.pi
        me.append(np.average(ME))

        Rhov2= density*(velx*velx+vely*vely+velz*velz)*0.5
        rhov2.append(np.average(Rhov2))
        
        Rhou=density*u
        rhou.append(np.average(Rhou))
        
        Mstress=-Bx*By/4./np.pi
        mstress.append(np.average(Mstress))
        
        Rstress=density*velx*vely
        rstress.append(np.average(Rstress))
        
        print i
        time.append(i*timestep)


plt.figure()
plt.plot(time,me)
plt.title("me")
plt.figure()
plt.plot(time,rhov2)
plt.title("ke")
plt.figure()
plt.plot(time,rhou)
plt.title("iu")
plt.figure()
plt.plot(time,mstress)
plt.title("MS")
plt.figure()
plt.plot(time,rstress)
plt.title("RS")
