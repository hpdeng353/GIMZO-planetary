import h5py
import numpy as np



s = h5py.File('/scratch/snx3000/hpdeng/gizmout/grvdisk1/init/grvdisk1-ics.hdf5','r+')


vel = s["PartType0"]["Velocities"][:]
pos = s["PartType0"]["Coordinates"][:]
del s["PartType0"]["Velocities"]


rr=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
sin = pos[:,1]/rr
cos = pos[:,0]/rr


vel[:,0] += 1/np.sqrt(rr)*sin
vel[:,1] -= 1/np.sqrt(rr)*cos



s.create_dataset("PartType0/Velocities",data=vel)
s.flush()
s.close()
            

