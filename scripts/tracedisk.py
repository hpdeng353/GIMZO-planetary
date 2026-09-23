# -*- coding: utf-8 -*-
import h5py
import numpy as np
import math as m
s = h5py.File('test.hdf5','r+')

rho0=5
rhoe=15.20
ids = s["PartType0"]["ParticleIDs"][:]  
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
rho = s["PartType0"]["Density"][:]
mat = s["PartType0"]["Materials"][:]

del s["PartType0"]["newpgrp"]
del s["PartType0"]["Velocities"]
Ngas=len(mass)
pgrp=np.zeros(Ngas)
pgrp[:]=1

#1 escape; 2 planet; 0 disk
pgrp[np.where(rho>rho0)]=2.
xx=pos[:,0]
rho[np.where(xx<-50)]=0.01
mass[np.where(xx<-50)]=0.0
pgrp[np.where(rho>rho0)]=2.
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
    jdisk=j2[np.where(pgrp==0)]
    jdisk=np.sqrt(jdisk)
    Ld=np.sum(mdisk*jdisk)/57.4
    print "LD=", Ld

mdstot=np.sum(mass[np.where((pgrp==0)&(mat==0))])
mdstarget=np.sum(mass[np.where((pgrp==0)&(ids<419339)&(mat==0))])

mpstot=np.sum(mass[np.where((pgrp==2)&(mat==0))])
mpstarget=np.sum(mass[np.where((pgrp==2)&(ids<419339)&(mat==0))])

ftd= mdstarget/mdstot
fte= mpstarget/mpstot
deltaft=ftd/fte -1
print "ftd",ftd, "fte",fte, "deltaft",deltaft
iron=mass[np.where((pgrp==0)&(mat==1))]
#iron=mass[np.where((pgrp==2)&(mat==1))]
print "iron fraction=", np.sum(iron)/Md



vel[:,0] -=vx0
s.create_dataset("PartType0/newpgrp",data=pgrp)
s.create_dataset("PartType0/Velocities",data=vel)
s.flush()
s.close()
            
id0=ids[np.where(pgrp==0)]
id1=ids[np.where(pgrp==1)]
id2=ids[np.where(pgrp==2)]
"""
startnum=0
step=1
numbsnaps=801
fname_base="snapshot_"
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s = h5py.File(fname,'r+')

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
#        del s["PartType0"]["newpgrp"]                      
        s.create_dataset("PartType0/newpgrp",data=pgrp1)
        s.flush()
        s.close()
        print i
"""