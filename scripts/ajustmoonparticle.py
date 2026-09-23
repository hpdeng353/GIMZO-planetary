import h5py
import numpy as np

s = h5py.File('/scratch/snx3000/hpdeng/gizmout/moon4/init/moon4-ics.hdf5','r+')


npart=s["Header"].attrs["NumPart_ThisFile"]
Npart=s["Header"].attrs["NumPart_Total"]
Ngas=npart[0]


pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]
u = s["PartType0"]["InternalEnergy"][:]
ids= s["PartType0"]["ParticleIDs"][:]

Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1]
Radius = np.sqrt(Radius)

gmass=s["PartType0"]["Masses"][:]







x=pos[:,0]
y=pos[:,1]
z=pos[:,2]
x1=x[np.where(ids>28059)]
y1=y[np.where(ids>28059)]
z1=z[np.where(ids>28059)]



vx=vel[:,0]
vy=vel[:,1]
vz=vel[:,2]
vx1=vx[np.where(ids>28059)]
vy1=vy[np.where(ids>28059)]
vz1=vz[np.where(ids>28059)]


u1=u[np.where(ids>28059)]
u=np.concatenate((u,u1))



newpart=len(anulus2)
ids1=np.arange(2000000,newpart+2000001)
ids2=np.arange(1000000,newpart+1000001)
print "ngas=, newpart=", len(gmass),newpart


ids=np.concatenate((ids,ids1))
ids=np.concatenate((ids,ids2))


gmass=np.concatenate((gmass,anulus1))
gmass=np.concatenate((gmass,anulus2))

vx=np.concatenate((vx,vx1))
vx=np.concatenate((vx,vx2))

vy=np.concatenate((vy,vy1))
vy=np.concatenate((vy,vy2))

vz=np.concatenate((vz,vz1))
vz=np.concatenate((vz,vz2))




x=np.concatenate((x,x1))
x=np.concatenate((x,x2))

y=np.concatenate((y,y1))
y=np.concatenate((y,y2))


z=np.concatenate((z,z1))
z=np.concatenate((z,z2))

npart = npart + [2*newpart,0,0,0,0,0]
Npart = Npart + [2*newpart,0,0,0,0,0]

s["Header"].attrs.modify("NumPart_ThisFile",npart)
s["Header"].attrs.modify("NumPart_Total",Npart)


pos=np.zeros((len(gmass),3))
pos[:,0]=x
pos[:,1]=y
pos[:,2]=z

vel=np.zeros((len(gmass),3))
vel[:,0]=vx
vel[:,1]=vy
vel[:,2]=vz


del s["PartType0"]["Masses"]
del s["PartType0"]["Coordinates"]
del s["PartType0"]["ParticleIDs"]
del s["PartType0"]["InternalEnergy"]
del s["PartType0"]["Velocities"]

del s["PartType1"]["Coordinates"]

s.create_dataset("PartType0/Masses",data=gmass)
s.create_dataset("PartType0/ParticleIDs",data=ids)
s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType0/InternalEnergy",data=u)
s.create_dataset("PartType0/Velocities",data=vel)

s.create_dataset("PartType1/Coordinates",data=pos1)
print "infall added"
    

s.flush()
s.close()
            

