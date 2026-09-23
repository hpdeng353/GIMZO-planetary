import h5py
import numpy as np

import numpy as np
import h5py as h5py
import math as m
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pyplot as plt

def taper1(x,x0,sharp): #default setting center of taper x0=3*rin
    f = 1.0 - 1.0/(m.exp((x-x0)/sharp)+1.)
    return f
def dr_taper1(x,x0,sharp):
    dfdr = 1.0/(m.exp((x-x0)/sharp)+1.)/(m.exp((x-x0)/sharp)+1.)/sharp*m.exp((x-x0)/sharp)
    return dfdr

def Polar2Cartisian(radius,theta):
    xx=np.zeros(len(radius))
    yy = np.zeros(len(radius))
    for i in range(0,len(radius)):
        xx[i] = radius[i] * m.cos(theta[i])
        yy[i] = radius[i] * m.sin(theta[i])
    return xx, yy

def Veltransform(vr,vt,theta):
    vx = np.zeros(len(theta))
    vy = np.zeros(len(theta))
    for i in range(len(theta)):
        #vx[i] = vr[i] * m.cos(theta[i]) - vt[i] * m.sin(theta[i])
        #vy[i] = vr[i] * m.sin(theta[i]) + vt[i] * m.cos(theta[i])
        vx[i] =  - vt[i] * m.sin(theta[i])
        vy[i] =  vt[i] * m.cos(theta[i])
    return vx, vy

    
Ngas= 5000000
T0=300
mu=1.0
gamma=5.0/3.0
selfgravity = 0 # selfgravity = 0 no dark particle in ic
StarMass = 10.0
Rho_Power = -3.0 # diff from sigma power rho_per
Temp_Power = -1.0
Rin = 0.2 #in code unit corespond to 1.25AU temp0=300
Rout = 10.0 #in code unit corespond to 62.5AU
X0=6*Rin
#    Sharp=0.05*X0
in_over_out = 0.1
Sharp = (Rin-X0)/np.log(1./(1.0-in_over_out*m.pow((Rout/Rin),Rho_Power))-1.)
diskmass = 0.1 #mass unit = 0.1 mass of sun
#    soft =0.02
V_k = m.sqrt(10.0/0.2)*m.sqrt(0.2/Rin*StarMass/10.0) # rin =0.2code unit =1.25AU, 7.07 in code unit
C_s =0.41762*m.sqrt(T0/300.0/mu) # corespond to 300 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!isothemal sound speed
beta =1.
betaphi = beta*100./99
betaz=100*beta
H0 = C_s / V_k * m.sqrt((1.0+beta)/beta) * Rin    #caution h0**2 = k*T0/mu*m_proton /omega_k**2 *(1+beta)/beta 0.037
Aspect0 =H0/Rin
rho0 = 4.43476867871


s = h5py.File("/zbox/data/hpdeng/gizmout/ngrvmhd2/snapshot_001.hdf5","r+")
density = s["PartType0"]["Density"][:]
pos = s["PartType0"]["Coordinates"][:]
mag = s["PartType0"]["Coordinates"][:]

for i in range(0,len(density)):
    Radius = m.sqrt(pos[i,0] * pos[i,0] + pos[i,1] * pos[i,1])
    sin = pos[i,1]/Radius
    cos = pos[i,0]/Radius
    press = density[i]* C_s * C_s * m.pow(Radius/Rin, Temp_Power)
    magt = m.sqrt(8.*3.141592/betaphi*34.55277*press)
    magz = m.sqrt(8.*3.141592/betaz*34.55277*press)
    mag[i,2] = magz
    mag[i,0] = -magt * sin
    mag[i,1] = magt * cos
del s["PartType0"]["MagneticField"]
s.create_dataset("PartType0/MagneticField",data=mag)
s.flush()
s.close()
            

