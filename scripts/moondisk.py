# -*- coding: utf-8 -*-
import h5py
import numpy as np
import math as m
s = h5py.File('snapshot_000.hdf5','r+')

rho0=15
rhoe=15.20
ids = s["PartType0"]["ParticleIDs"][:]  
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
rho = s["PartType0"]["Density"][:]
Ngas=len(mass)
pgrp=np.zeros(Ngas)
pgrp[:]=1
#1 escape; 2 planet; 0 disk
pgrp[np.where(rho>rho0)]=2.
xx=pos[:,0]
yy=pos[:,1]
zz=pos[:,2]
vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]


#subgrp_planet=pgrp[(pgrp>0.5)&(pgrp<1.5)]
#Nplanet=
print pgrp
mm=mass[np.where(pgrp==2)]
print mm,len(mm)
Mp=np.sum(mm)
print Mp 
Mp0=0.

while (np.fabs(Mp-Mp0)>0.001):
    Mp0=Mp
    rp=np.cbrt(3*Mp/(4.*np.pi*rhoe))
    xx_p=xx[np.where(pgrp>1)]
    yy_p=yy[np.where(pgrp>1)]
    zz_p=zz[np.where(pgrp>1)]

    x0=np.sum(xx_p*mm)/Mp
    y0=np.sum(yy_p*mm)/Mp
    z0=np.sum(zz_p*mm)/Mp


    vx_p=vx[np.where(pgrp>1)]
    vy_p=vy[np.where(pgrp>1)]
    vz_p=vz[np.where(pgrp>1)]

    vx0=np.sum(vx_p*mm)/Mp
    vy0=np.sum(vy_p*mm)/Mp
    vz0=np.sum(vz_p*mm)/Mp

    xxr=xx-x0
    yyr=yy-y0
    zzr=zz-z0
    vxr=vx-vx0
    vyr=vy-vy0
    vzr=vz-vz0


    Radius=np.sqrt(xxr*xxr + yyr*yyr + zzr*zzr)

    E=0.5*(vxr*vxr+vyr*vyr+vzr*vzr)-Mp/Radius

    pgrp[:]=1.

    pgrp[np.where(E<0)]=2.
# only for particles with E<0

    majora= -Mp/(2.0*E)
    j2=(xxr*vyr-yyr*vxr)*(xxr*vyr-yyr*vxr)+(yyr*vzr-zzr*vyr)*(yyr*vzr-zzr*vyr)+(xxr*vzr-zzr*vxr)*(xxr*vzr-zzr*vxr)
    ecent=np.sqrt(1.-j2/(Mp*majora))
#if E>0 we set artificially majora>0 to bypass majora*(1-ecent) argument
    print majora,
    print ecent
    pgrp[np.where((majora*(1-ecent)>rp)&(E<0)&(rho<50.3))]=0.

    mm=mass[np.where(pgrp>1)]
    Mp=np.sum(mm)

    mdisk=mass[np.where(pgrp<1)]
    Md=np.sum(mdisk)
    print "pos=",x0,y0,z0
    print "vel=", vx0,vy0,vz0
    print "Mp=", Mp,"Rp",rp
    print "Md=",Md,"Md/Mp=",Md/Mp, "Md/Ml=", Md/0.769

mmd1=mass[np.where((pgrp==0)&(ids<443124))]
#mmd2=mass[np.where((pgrp==0)&(ids<443124))]
mme1=mass[np.where((pgrp==2)&(ids<443124))]

ft= np.sum(mmd1)/Md
deltaft=ft/(np.sum(mme1)/Mp)-1
print "ft",ft, "deltaft",deltaft
s.create_dataset("PartType0/pgrp",data=pgrp)
s.flush()
s.close()
            

