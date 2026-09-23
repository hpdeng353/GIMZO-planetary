import h5py
import numpy as np


beta=25
gamma_eos=5.0/3.
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta



s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvmhd5/init/grvmhd5-ics.hdf5','r+')

mag = s["PartType0"]["MagneticField"][:]
pos = s["PartType0"]["Coordinates"][:]
u = s["PartType0"]["InternalEnergy"][:]
del s["PartType0"]["MagneticField"]
rho = s["PartType0"]["Density"][:]
vel=s["PartType0"]["Velocities"][:]
"""
pos1= s["PartType1"]["Coordinates"][:]
vel1=s["PartType1"]["Velocities"][:]
ids1= s["PartType1"]["Coordinates"][:]

ids1[:]=50000000
#print pos1[0][0],pos1[0][1],pos1[0][2]
#print pos[:,0]

pos[:,0] -= pos1[0][0]
pos[:,1] -= pos1[0][1]
pos[:,2] -= pos1[0][2]
#print pos[:,0]
pos1[:,:]=0



vel[:,0] -= vel1[0][0]
vel[:,1] -= vel1[0][1]
vel[:,2] -= vel1[0][2]

"""
mag=vel[:,:]

#mag=s["PartType0"]["MagneticField"][:]
#del s["PartType0"]["MagneticField"]
rr=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
array = np.asarray(rr)
order = array.argsort()
ranks = order.argsort()

sin = pos[:,1]/rr
cos = pos[:,0]/rr





#cs2 = u*(gamma_eos-1)*gamma_eos
#rho1=rho[np.where((rr>8)&(rr<10)&(pos[:,2]<0.01)&(pos[:,2]>-0.01))]
#cs21=cs2[np.where((rr>8)&(rr<10)&(pos[:,2]<0.01)&(pos[:,2]>-0.01))]
p_fiducial=rho*u*(gamma_eos-1)*5.27*1e6#np.average(rho1*cs21)*5.27*1e6
bample=np.sqrt(p_fiducial/beta*8*np.pi)

bx= - bample * sin 
by= bample * cos
#mag0[np.where(rr<5)]=bx[np.where(rr<5)]
#mag0[np.where((pos[:,2]<0.1)&(pos[:,2]>-0.1)&(rr<5))]=bx[np.where((pos[:,2]<0.1)&(pos[:,2]>-0.1)&(rr<5))]
#mag1[np.where(rr<5)]=by[np.where(rr<5)]
#mag1[np.where((pos[:,2]<0.1)&(pos[:,2]>-0.1)&(rr<5))]=by[np.where((pos[:,2]<0.1)&(pos[:,2]>-0.1)&(rr<5))]
#mag2[np.where(rr<5)]=0
bx[np.where(rho<1e-6)]==0
by[np.where(rho<1e-6)]==0

#bx=mag[:,0]
#by=mag[:,1]
bz=0.04*np.power(rr/0.05,-1.5)
#bx[np.where(rr<10)]=0
#by[np.where(rr<10)]=0
#bz[np.where(rr<10)]=0

#del s["PartType1"]["Velocities"]
#del s["PartType0"]["Velocities"]
#del s["PartType0"]["ParticleIDs"]
#del s["PartType1"]["ParticleIDs"]
#del s["PartType0"]["Coordinates"]
#mag[:,0]=bx
#mag[:,1]=by
mag[:,2]+=bz

u*=88.
s.create_dataset("PartType0/MagneticField",data=mag)
#s.create_dataset("PartType0/InternalEnergy",data=u)

#s.create_dataset("PartType0/Coordinates",data=pos)
#s.create_dataset("PartType0/ParticleIDs",data=ranks)
#s.create_dataset("PartType1/ParticleIDs",data=ids1)
#s.create_dataset("PartType1/Coordinates",data=pos1)
#s.create_dataset("PartType1/Velocities",data=vel1)
#s.create_dataset("PartType0/Velocities",data=vel)
s.flush()
s.close()
            

