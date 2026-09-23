import h5py
import numpy as np

Lx=1.414
Lz=24.
P0=4*np.pi
R0=Lx/4.
betay=25
betap=1600
gamma_eos=5.0/3.
#B0=np.sqrt(8.*np.pi/400)
#beta=B0*B0/8./np.pi
#beta=1./beta
ample1=0.1


s = h5py.File('hrisph.hdf5','r+')


vort = s["PartType0"]["Vorticity"][:]
u = s["PartType0"]["InternalEnergy"][:]
#hsml=s["PartType0"]["SmoothingLength"][:]
u = np.abs(vort[:,2])
#u =np.sqrt(vort[:,0]*vort[:,0]+vort[:,1]*vort[:,1]+vort[:,2]*vort[:,2])

del s["PartType0"]["InternalEnergy"]
#del s["PartType0"]["SmoothingLength"]

s.create_dataset("PartType0/InternalEnergy",data=u)
#s.create_dataset("PartType0/SmoothingLength",data=hsml)

s.flush()
s.close()
            

