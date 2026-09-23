import h5py
import numpy as np

b=0.574
#deltav=12.
s = h5py.File('sgaia.hdf5','r+')
s1= h5py.File('stheia.hdf5','r+')


npart=s["Header"].attrs["NumPart_ThisFile"]
Npart=s["Header"].attrs["NumPart_Total"]

npart1=s1["Header"].attrs["NumPart_ThisFile"]
Npart1=s1["Header"].attrs["NumPart_Total"]

Ngas=npart[0]
Ngas1=npart1[0]

npart = npart + npart1
Npart = Npart + Npart1


print "nags=",Ngas,"ngas1=",Ngas1
s["Header"].attrs.modify("NumPart_ThisFile",npart)
s["Header"].attrs.modify("NumPart_Total",Npart)

id=s["PartType0"]["ParticleIDs"][:]
id1=s1["PartType0"]["ParticleIDs"][:]
id1 += Ngas #  !!!!! two new sphere need this.0-N1 0-N2

ids=np.concatenate((id,id1))

chi=s["PartType0"]["ParticleChildIDsNumber"][:]
chi1=s1["PartType0"]["ParticleChildIDsNumber"][:]
chi=np.concatenate((chi,chi1))

gen=s["PartType0"]["ParticleIDGenerationNumber"][:]
gen1=s1["PartType0"]["ParticleIDGenerationNumber"][:]
gen=np.concatenate((gen,gen1))

mat=s["PartType0"]["Materials"][:]
mat1=s1["PartType0"]["Materials"][:]

mat=np.concatenate((mat,mat1))

rho=s["PartType0"]["Density"][:]
rho1=s1["PartType0"]["Density"][:]
rho=np.concatenate((rho,rho1))


vel=s["PartType0"]["Velocities"][:]
vel1=s1["PartType0"]["Velocities"][:]

pos = s["PartType0"]["Coordinates"][:]
pos1 = s1["PartType0"]["Coordinates"][:]

x0=np.average(pos[:,0])
y0=np.average(pos[:,1])
z0=np.average(pos[:,2])
x1=np.average(pos1[:,0])
y1=np.average(pos1[:,1])
z1=np.average(pos1[:,2])

pos[:,0]-=x0
pos[:,1]-=y0
pos[:,2]-=z0
pos1[:,0]-=x1
pos1[:,1]-=y1
pos1[:,2]-=z1

Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1] + pos[:,2]*pos[:,2]
Radius = np.sqrt(Radius)
rr=np.max(Radius)

Radius1 = pos1[:,0] * pos1[:,0] + pos1[:,1] * pos1[:,1] + pos1[:,2]*pos1[:,2]
Radius1 = np.sqrt(Radius1)
rr1=np.max(Radius1)

print "Re=",rr, "Rt=",rr1


u=s["PartType0"]["InternalEnergy"][:]
u1=s1["PartType0"]["InternalEnergy"][:]
uv_g= np.concatenate((u,u1))


mass=s["PartType0"]["Masses"][:]
mass1=s1["PartType0"]["Masses"][:]
#mass1[np.where(mat1==0)]=np.min(mass)
#mass1[:]=mass[0]
Masses= np.concatenate((mass,mass1))
mass=np.sum(mass)
mass1=np.sum(mass1)



gamma=mass1/(mass+mass1)

print "mass=", mass1/62.46, mass/62.46
print "gamma=", gamma

deltax=(rr+rr1)*1.02*np.sqrt(1.-b*b)
deltay=(rr+rr1)*b

xx= -deltax*gamma
xx1= deltax*(1-gamma)
yy= -deltay*gamma
yy1=deltay*(1-gamma)

"""
xx=-0.179
yy=-0.155
xx1=1.2
yy1=1.04
"""

pos[:,0]+=xx 
pos[:,1]+= yy
pos1[:,0] += xx1
pos1[:,1] += yy1

pos = np.concatenate((pos,pos1),axis=0)

deltav=np.sqrt(1.1/(rr+rr1))*11.18*1.2
#deltav=np.sqrt(1.1/1.62)*11.18*1.2
#deltav=9.2*1.15
vv= deltav*gamma
vv1= -deltav*(1-gamma)

#vv=1.085
#vv1=-7.26

vx=np.average(vel[:,0])
vx1=np.average(vel1[:,0])

vel[:,0]-=vx
vel1[:,0]-=vx1

vel[:,0] += vv
vel1[:,0] += vv1

print vv,vv1
vel = np.concatenate((vel,vel1),axis=0)

Limp=-np.sum(vel[:,0]*pos[:,1]*Masses[:])/57.4
print "Limp=", Limp
del s["PartType0"]["ParticleIDs"]
del s["PartType0"]["Materials"]
del s["PartType0"]["Masses"]
del s["PartType0"]["Coordinates"]
del s["PartType0"]["Velocities"]
del s["PartType0"]["InternalEnergy"]
del s["PartType0"]["Density"]
del s["PartType0"]["ParticleChildIDsNumber"]
del s["PartType0"]["ParticleIDGenerationNumber"]

s.create_dataset("PartType0/ParticleIDs",data=ids)
s.create_dataset("PartType0/Materials",data=mat)
s.create_dataset("PartType0/Masses",data=Masses)
s.create_dataset("PartType0/Coordinates",data=pos)
s.create_dataset("PartType0/Velocities",data=vel)
s.create_dataset("PartType0/InternalEnergy",data=uv_g)
s.create_dataset("PartType0/Density",data=rho)
s.create_dataset("PartType0/ParticleChildIDsNumber",data=chi)
s.create_dataset("PartType0/ParticleIDGenerationNumber",data=gen)
    

s.flush()
s.close()
            

