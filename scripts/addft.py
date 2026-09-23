import h5py
import numpy as np

Lx=1.
Lz=1.
P0=4*np.pi
R0=Lx/4.
betay=100
betap=1600
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta



s = h5py.File('/scratch/snx2000tds/hpdeng/gizmout/boxtest3/snapshot_000.hdf5','r+')
pos = s["PartType0"]["Coordinates"][:]
mag = s["PartType0"]["Velocities"][:]

del s["PartType0"]["MagneticField"]
del s["PartType0"]["Coordinates"]

pos[:,0] -= Lx/2.
pos[:,2] -= Lz/2.

rr=np.sqrt(pos[:,0]*pos[:,0] + pos[:,2]*pos[:,2])
print min(rr)

mag[:,0] = - np.sqrt(2*P0/betap)*np.sin(np.pi*rr/R0)*pos[:,2]/rr
mag[:,2] =  np.sqrt(2*P0/betap)*np.sin(np.pi*rr/R0)*pos[:,0]/rr
mag[:,1] = np.sqrt(2*P0/betay - 2*P0/betap*np.sin(np.pi*rr/R0)*np.sin(np.pi*rr/R0))

for i in range(0,len(rr)):
    if (i % 500000) == 0:
        print i
    if rr[i] > R0 :
        mag[i,0] = 0.
        mag[i,1] = 0.
        mag[i,2] = 0.


                                                         

pos[:,2] += Lz/2.
pos[:,0] += Lx/2.

s.create_dataset("PartType0/MagneticField",data=mag)
s.create_dataset("PartType0/Coordinates",data=pos)

s.flush()
s.close()
            

