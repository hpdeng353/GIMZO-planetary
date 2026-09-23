# -*- coding: utf-8 -*-
import h5py
import numpy as np
import math as m
s = h5py.File('snapshot_000.hdf5','r+')

rho0=8.
rhoe=15.20
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
rho = s["PartType0"]["Density"][:]
Ngas=len(mass)
pgrp=np.zeros(Ngas)
#1 escape; 0 planet; 2 disk
pgrp[:]=2.
pgrp[np.where(rho>rho0)]=0.
xx=pos[:,0]
yy=pos[:,1]
zz=pos[:,2]
vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]


#subgrp_planet=pgrp[(pgrp>0.5)&(pgrp<1.5)]
#Nplanet=
print pgrp
mm=mass[np.where(pgrp<0.5)]
print mm,len(mm)
Mp=np.sum(mm)
print Mp 
Mp0=0.

while (np.fabs(Mp-Mp0)>0.001):
    Mp0=Mp
    rp=np.cbrt(3*Mp/(4.*np.pi*rhoe))
    xx_p=xx[np.where(pgrp<0.5)]
    yy_p=yy[np.where(pgrp<0.5)]
    zz_p=zz[np.where(pgrp<0.5)]
    x0=np.sum(xx_p*mm)/Mp
    y0=np.sum(yy_p*mm)/Mp
    z0=np.sum(zz_p*mm)/Mp
    vx_p=vx[np.where(pgrp<0.5)]
    vy_p=vy[np.where(pgrp<0.5)]
    vz_p=vz[np.where(pgrp<0.5)]
    vx0=np.sum(vx_p*mm)/Mp
    vy0=np.sum(vy_p*mm)/Mp
    vz0=np.sum(vz_p*mm)/Mp
    for i in range(0,Ngas):
        if pgrp[i] != 0:
            xxr=xx[i]-x0
            yyr=yy[i]-y0
            zzr=zz[i]-z0
            vxr=vx[i]-vx0
            vyr=vy[i]-vy0
            vzr=vz[i]-vz0
            Radius=np.sqrt(xxr*xxr + yyr*yyr + zzr*zzr)
            if Radius < rp:
                pgrp[i]=0.
            else:
                E=0.5*(vxr*vxr+vyr*vyr+vzr*vzr)-Mp/Radius
                if E>=0:
                    pgrp[i]=1.
                else:
                    pgrp[i]=0.
                    majora= -Mp/(2.0*E)
                    j2=(xxr*vyr-yyr*vxr)*(xxr*vyr-yyr*vxr)+(yyr*vzr-zzr*vyr)*(yyr*vzr-zzr*vyr)+(xxr*vzr-zzr*vxr)*(xxr*vzr-zzr*vxr)
                    ecent=np.sqrt(1.-j2/(Mp*majora))
                    if majora*(1-ecent)>rp:
                        pgrp[i]=2.0
    
    mm=mass[np.where(pgrp<0.5)]
    Mp=np.sum(mm)
    mdisk=mass[np.where(pgrp>1.5)]
    Md=np.sum(mdisk)
    print "pos=",x0,y0,z0
    print "vel=", vx0,vy0,vz0
    print "Mp=", Mp,"Rp",rp
    print "Md=",Md,"Md/Mp=",Md/Mp

s.create_dataset("PartType0/pgrp",data=pgrp)
"""
s.create_dataset("PartType0/E",data=E)
s.create_dataset("PartType0/major",data=majora)
s.create_dataset("PartType0/j2",data=j2)
s.create_dataset("PartType0/ecent",data=ecent)
"""
s.flush()
s.close()
            

