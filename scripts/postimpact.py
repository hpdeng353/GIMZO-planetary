# -*- coding: utf-8 -*-
import h5py
import numpy as np
import math as m
import matplotlib.pyplot as plt
s = h5py.File('test.hdf5','r+')
R0=1.
id0=1677364#419339#### ##455176
rho0=5
rhoe=15.20
entropy=s["PartType0"]["Entropy"][:]  
ids = s["PartType0"]["ParticleIDs"][:]  
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
rho = s["PartType0"]["Density"][:]
u = s["PartType0"]["InternalEnergy"][:]
mat = s["PartType0"]["Materials"][:]
e = "/PartType0/newpgrp" in s
if e :
    del s["PartType0"]["newpgrp"]
#del s["PartType0"]["Velocities"]
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
mdstarget=np.sum(mass[np.where((pgrp==0)&(ids<id0)&(mat==0))])

#mpstot=np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R0))])
#mpstarget=np.sum(mass[np.where((pgrp==2)&(ids<id0)&(mat==0)&(Radius>R0))])
mpstot=np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R0))])
mpstarget=np.sum(mass[np.where((pgrp==2)&(ids<id0)&(mat==0)&(Radius>R0))])
ftd= mdstarget/mdstot
fte= mpstarget/mpstot
deltaft=ftd/fte -1
print "ftd",ftd, "fte",fte, "deltaft",deltaft
iron=mass[np.where((pgrp==0)&(mat==1))]
#iron=mass[np.where((pgrp==2)&(mat==1))]
print "iron fraction=", np.sum(iron)/Md
print "upper mantle=",np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R0))])


vel[:,0] -=vx0
s.create_dataset("PartType0/newpgrp",data=pgrp)
#s.create_dataset("PartType0/Velocities",data=vel)
s.flush()
s.close()

bins=np.arange(0,1.4,0.01)
binentropy=np.arange(0,1.4,0.01)
binmixing=np.arange(0,1.4,0.01)
binu=np.arange(0,1.4,0.01)
binrho=np.arange(0,1.4,0.01)
for i in range(0,len(bins)-1):
    binentropy[i]=np.average(entropy[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binu[i]=np.average(u[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binrho[i]=np.average(rho[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binmixing[i]=np.sum(mass[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1])&(ids<id0))])/np.sum(mass[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
binmixing[-1]=binmixing[-2]
binentropy[-1]=binentropy[-2]
binu[-1]=binu[-2]
binrho[-1]=binrho[-1]

#plt.plot(bins,binu)
plt.plot(bins,binrho)    
#plt.plot(bins,binentropy)
#plt.plot(bins,binmixing)

plt.show()