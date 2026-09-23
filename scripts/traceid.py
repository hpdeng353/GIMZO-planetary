import h5py
import numpy as np
import math as m
s = h5py.File('snapshot_000.hdf5','r+')

id0=455176#1677364#419339 # ## #
ids = s["PartType0"]["ParticleIDs"][:]  
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
rho = s["PartType0"]["Density"][:]
mat = s["PartType0"]["Materials"][:]

e = "/PartType0/newpgrp" in s
if e :
    del s["PartType0"]["newpgrp"]

Ngas=len(mass)
pgrp=np.zeros(Ngas)
pgrp[:]=1

#mass[np.where(mass>0.00022)]=0.0001

xx=pos[:,0]
xx1=xx[np.where(ids<id0)]
x0= np.average(xx1)
xx1-= x0
xx-=x0

yy=pos[:,1]
yy1=yy[np.where(ids<id0)]
y0=np.average(yy1)
yy1-=y0
yy-=y0

zz=pos[:,2]
zz1=zz[np.where(ids<id0)]

rr1 = np.sqrt(xx1*xx1+yy1*yy1+zz1*zz1)
rr = np.sqrt(xx*xx+yy*yy+zz*zz)

radiustarget=max(rr1)
radiuscore=max(rr[np.where((mat==1)&(ids<id0))])

rmid=(radiuscore+radiustarget)/2

pgrp[np.where(rr>radiuscore)]=2
pgrp[np.where(rr>rmid)]=2
pgrp[np.where((ids>id0)&(mat==1))]=3
pgrp[np.where((ids>id0)&(mat==0))]=4








s.create_dataset("PartType0/newpgrp",data=pgrp)
s.flush()
s.close()
            

startnum=800
step=1
numbsnaps=100000000
fname_base="snapshot_"
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s = h5py.File(fname,'r+')
        pos = s["PartType0"]["Coordinates"][:]
        vel = s["PartType0"]["Velocities"][:]
        rho = s["PartType0"]["Density"][:]
        del  s["PartType0"]["Coordinates"]
        
        xx=pos[:,0]
        yy=pos[:,1]
        zz=pos[:,2]
        pos[:,0]-=np.average(xx[np.where(rho>31)])
        pos[:,1]-=np.average(yy[np.where(rho>31)])
        pos[:,2]-=np.average(zz[np.where(rho>31)])
        ids1=s["PartType0"]["ParticleIDs"][:]
        pgrp1=pgrp
        
        temp=np.zeros((2,len(ids)))
        temp1=np.zeros((2,len(ids)))
        temp[0,:]=ids 
        temp[1,:]=pgrp

        temp1[0,:]=ids1
        temp1[1,:]=pgrp1
        
        ii=np.argsort(ids)
        jj=np.argsort(ids1)

        temp1[1,jj]=temp[1,ii]

        pgrp1=temp1[1,:]
        e = "/PartType0/newpgrp" in s
        if e :
            del s["PartType0"]["newpgrp"]

        s.create_dataset("PartType0/newpgrp",data=pgrp1)
        s.create_dataset("PartType0/Coordinates",data=pos)

        s.flush()
        s.close()
        print i