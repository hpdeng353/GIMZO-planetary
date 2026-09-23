# -*- coding: utf-8 -*-
import h5py
import numpy as np
import math as m
import matplotlib.pyplot as plt
s = h5py.File('test.hdf5','r+')
R1=1.01
R0=0.83
id0=455176#419339######1677364####### ##
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
#pgrp[np.where(rho>rho0)]=2.
xx=pos[:,0]
#rho[np.where(xx<-50)]=0.01
#mass[np.where(xx<-50)]=0.0
pgrp[np.where(rho>rho0)]=2.
yy=pos[:,1]
zz=pos[:,2]
vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]


#subgrp_planet=pgrp[(pgrp>0.5)&(pgrp<1.5)]
#Nplanet=
print pgrp
mp=mass[np.where(pgrp==2)]
print mp,len(mp)
Mp=np.sum(mp)
print Mp 
Mp0=0.
f=0.1
while (np.fabs(Mp-Mp0)>0.001):
    Mp0=Mp
    rp=np.cbrt(3*Mp/(4.*np.pi*rhoe*(1-f)))

    xx_p=xx[np.where(pgrp>1)]
    yy_p=yy[np.where(pgrp>1)]
    zz_p=zz[np.where(pgrp>1)]

    x0=np.sum(xx_p*mp)/np.sum(mp)
    y0=np.sum(yy_p*mp)/np.sum(mp)
    z0=np.sum(zz_p*mp)/np.sum(mp)


    vx_p=vx[np.where(pgrp>1)]
    vy_p=vy[np.where(pgrp>1)]
    vz_p=vz[np.where(pgrp>1)]

    vx0=np.sum(vx_p*mp)/np.sum(mp)
    vy0=np.sum(vy_p*mp)/np.sum(mp)
    vz0=np.sum(vz_p*mp)/np.sum(mp)

    xxr=xx-x0
    yyr=yy-y0
    zzr=zz-z0
    vxr=vx-vx0
    vyr=vy-vy0
    vzr=vz-vz0


    Radius=np.sqrt(xxr*xxr + yyr*yyr + zzr*zzr)
    Radius1=np.sqrt(xxr*xxr + yyr*yyr)
    E=0.5*(vxr*vxr+vyr*vyr+vzr*vzr)-Mp/Radius

    pgrp[:]=1.

    pgrp[np.where(E<0)]=2.
# only for particles with E<0

    majora= -Mp/(2.0*E)
    j2=(xxr*vyr-yyr*vxr)*(xxr*vyr-yyr*vxr)+(yyr*vzr-zzr*vyr)*(yyr*vzr-zzr*vyr)+(xxr*vzr-zzr*vxr)*(xxr*vzr-zzr*vxr)
    jz=(xxr*vyr-yyr*vxr)*(xxr*vyr-yyr*vxr)
    jz=np.sqrt(jz)

    ecent=np.sqrt(1.-j2/(Mp*majora))
#if E>0 we set artificially majora>0 to bypass majora*(1-ecent) argument
    print majora,
    print ecent
    pgrp[np.where((majora*(1-ecent)>rp)&(E<0)&(rho<50.3))]=0.

    mp=mass[np.where(pgrp>1)]

    Mp=np.sum(mp)

    mdisk=mass[np.where(pgrp<1)]
    Md=np.sum(mdisk)

    jdisk=jz[np.where(pgrp==0)]


    Jpz1=jz[np.where((pgrp!=1)&(ids>id0)&(mat==0))]
    Jpz2=jz[np.where((pgrp!=1)&(ids<id0)&(mat==0))]
    
    Jp=jz[np.where(pgrp==2)]

    Mp1=mass[np.where((pgrp!=1)&(ids>id0)&(mat==0))]           
    Mp2=mass[np.where((pgrp!=1)&(ids<id0)&(mat==0))] 

    
    Lp1=np.sum(Mp1*Jpz1)/57.4
    Lp2=np.sum(Mp2*Jpz2)/57.4
    Lp=np.sum(mp*Jp)/57.4
    Ld=np.sum(mdisk*jdisk)/57.4

    f=2.5*(Lp/3.2)*(Lp/3.2)/(1+(2.5-3.75*0.35)**2)

    print "pos=",x0,y0,z0
    print "vel=", vx0,vy0,vz0
    print "Mp=", Mp,"Rp",rp,"f",f
    print "Md=",Md,"Md/Mp=",Md/Mp, "Md/Ml=", Md/0.769
    print "LD=", Ld
    print "L_red=",Lp1, "L_green=",Lp2,"Lp=",Lp



mdstot=np.sum(mass[np.where((pgrp==0)&(mat==0))])
mdstarget=np.sum(mass[np.where((pgrp==0)&(ids<id0)&(mat==0))])

#mpstot=np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R0))])
#mpstarget=np.sum(mass[np.where((pgrp==2)&(ids<id0)&(mat==0)&(Radius>R0))])
mpstot=np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R1))])
mpstarget=np.sum(mass[np.where((pgrp==2)&(ids<id0)&(mat==0)&(Radius>R1))])
ftd= mdstarget/mdstot
fte= mpstarget/mpstot
deltaft=ftd/fte -1
print "ftd",ftd, "fte",fte, "deltaft",deltaft
iron=mass[np.where((pgrp==0)&(mat==1))]
#iron=mass[np.where((pgrp==2)&(mat==1))]
print "iron fraction=", np.sum(iron)/Md
print "upper mantle=",np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R1))])/np.sum(mass[np.where((pgrp==2)&(mat==0))])
print "upper layer=",np.sum(mass[np.where((pgrp==2)&(mat==0)&(Radius>R0))])/np.sum(mass[np.where((pgrp==2)&(mat==0))])
print "mantle mass=",np.sum(mass[np.where((pgrp==2)&(mat==0))])

pdpth=np.average(Radius[np.where((pgrp==2)&(ids>id0)&(mat==0)&(Radius<1.3))])

mass1000=np.sum(Radius[np.where((pgrp==2)&(ids>id0)&(mat==0)&(Radius>R1))])
massred=np.sum(Radius[np.where((pgrp==2)&(ids>id0)&(mat==0))])

print "fraction 1000=", mass1000/massred

print "penetration=", pdpth
corer=np.max(Radius[np.where((pgrp==2)&(Radius<1.3))])
print "core radius=", corer
vel[:,0] -=vx0
s.create_dataset("PartType0/newpgrp",data=pgrp)
#s.create_dataset("PartType0/Velocities",data=vel)
s.flush()
s.close()

bins=np.arange(0,1.3,0.01)
binentropy=np.arange(0,1.3,0.01)
binmixing=np.arange(0,1.3,0.01)
binu=np.arange(0,1.3,0.01)
binrho=np.arange(0,1.3,0.01)
binjz=np.arange(0,1.3,0.01)
binimpf=np.arange(0,1.3,0.01)
binimpfc=np.arange(0,1.3,0.01)
bincore=np.arange(0,1.3,0.01)
for i in range(0,len(bins)-1):
    binentropy[i]=np.average(entropy[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binu[i]=np.average(u[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binjz[i]=np.average(jz[np.where((pgrp==2)&(Radius1>bins[i])&(Radius1<bins[i+1]))])
    binrho[i]=np.average(rho[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binmixing[i]=np.sum(mass[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1])&(ids<id0))])/np.sum(mass[np.where((pgrp==2)&(Radius>bins[i])&(Radius<bins[i+1]))])
    binimpf[i]=np.sum(mass[np.where((pgrp==2)&(Radius<bins[i+1])&(Radius>bins[i])&(ids>id0)&(mat==0))])/62.46/0.2/0.7
    binimpfc[i]=np.sum(mass[np.where((pgrp==2)&(Radius<bins[i+1])&(ids>id0)&(mat==0))])/62.46/0.2/0.7
    bincore[i]=np.sum(mass[np.where((pgrp==2)&(Radius<bins[i+1])&(mat==1))])/62.46/0.3

binmixing[-1]=binmixing[-2]
binentropy[-1]=binentropy[-2]
binu[-1]=binu[-2]
binrho[-1]=binrho[-2]
binjz[-1]=binjz[-2]
binimpf[-1]=binimpf[-2]
bincore[-1]=bincore[-2]
#plt.plot(bins,binu)
#plt.plot(bins,binrho)    
#plt.plot(bins,binjz)    
#plt.plot(bins,3*bins*bins)
#plt.plot(bins,binentropy)
#plt.plot(bins,binmixing)

#plt.plot(bins,binimpf)
#plt.show()
