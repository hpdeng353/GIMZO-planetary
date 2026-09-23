import h5py
import numpy as np
import matplotlib.pyplot as plt

#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta
ample=0.

#s = h5py.File('snapshot_1400.hdf5','r+')
s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk1/snapshot_100.hdf5','r+')




del s["PartType0"]["ParticleIDs"]

pos=s["PartType0"]["Coordinates"]
Radius=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
array = np.asarray(Radius)
order = array.argsort()
ranks = order.argsort()



print max(ranks)
s.create_dataset("PartType0/ParticleIDs",data=ranks)


s.flush()
s.close()
            
