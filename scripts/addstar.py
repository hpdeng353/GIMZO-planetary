import h5py
import numpy as np

s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk6/init/grvdisk6-ics.hdf5','r+')


npart=s["Header"].attrs["NumPart_ThisFile"]
Npart=s["Header"].attrs["NumPart_Total"]
Ngas=npart[0]

npart = npart + [0,1,0,0,0,0]
Npart = Npart + [0,1,0,0,0,0]

s["Header"].attrs.modify("NumPart_ThisFile",npart)
s["Header"].attrs.modify("NumPart_Total",Npart)

rho = s["PartType0"]["Density"][:]
mass = s["PartType0"]["Masses"][:]
pos = s["PartType0"]["Coordinates"][:]

xx=pos[:,0]
yy=pos[:,1]
zz=pos[:,2]
xx1=xx[np.where(rho>1)]
yy1=yy[np.where(rho>1)]
zz1=zz[np.where(rho>1)]
x0=np.average(xx1)
y0=np.average(yy1)
z0=np.average(zz1)

pos1=[[x0,y0,z0],]

pos = s["PartType0"]["Velocities"][:]

xx=pos[:,0]
yy=pos[:,1]
zz=pos[:,2]
xx1=xx[np.where(rho>1)]
yy1=yy[np.where(rho>1)]
zz1=zz[np.where(rho>1)]
x0=np.average(xx1)
y0=np.average(yy1)
z0=np.average(zz1)

vel1=[[x0,y0,z0],]


starmass=np.sum(mass[np.where(rho>1)])
mass[np.where(rho>1)]=0
ids = np.asarray([50000001,])
mass1 = np.asarray([starmass,])

print mass1[0]


del s["PartType0"]["Masses"]



s.create_dataset("PartType1/Masses",data=mass1)
s.create_dataset("PartType1/ParticleIDs",data=ids)
s.create_dataset("PartType1/Coordinates",data=pos1)
s.create_dataset("PartType1/Velocities",data=vel1)

s.create_dataset("PartType0/Masses",data=mass)

print "one dark particle added"
    

s.flush()
s.close()
            

