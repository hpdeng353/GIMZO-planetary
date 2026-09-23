import h5py
import numpy as np
import math as m
s = h5py.File('/users/hpdeng/mysoft/pynbody/pynbody/scripts/disk.hdf5','r+')

case ="pertrho"
mass = s["PartType0"]["Masses"][:]  
pos = s["PartType0"]["Coordinates"][:]
vel = s["PartType0"]["Velocities"][:]

T0=300
mu=1.0
StarMass = 10.0
Rho_Power = -3.0
Temp_Power = -0.9
Rin = 7 #in code unit corespond to 1.25AU temp0=300
Rout = 30.0 #in code unit corespond to 62.5AU
V_k = m.sqrt(10.0/0.2)*m.sqrt(0.2/Rin*StarMass/10.0) # rin =0.2code unit =1.25AU, 7.07 in code unit
C_s =m.sqrt(0.0015)
beta=20
betaz =200.

H0 = C_s / V_k * m.sqrt((1.0+beta)/beta) * Rin    #caution h0**2 = k*T0/mu*m_proton /omega_k**2 *(1+beta)/beta 0.037
C_s =m.sqrt(0.0015)
kz=0. #2*np.pi/0.1

kr=2*np.pi/25.*5
mphi=3.
ample= 0.05
r1=0.3


def Cartisian2Polar(vx,vy,cos_v,sin_v):
    vr=np.zeros(len(vx))
    vphi=np.zeros(len(vx))
    #vr=vx*cos_v + vy*sin_v
    vphi = -vx*sin_v + vy*cos_v
    return vr,vphi
    
def Polar2Cartisian(vr,vphi,cos_v,sin_v):
    vx=np.zeros(len(vr))
    vy=np.zeros(len(vr))
    vx=vr*cos_v - vphi*sin_v
    vy=vr*sin_v + vphi*cos_v
    return vx,vy


#del s["PartType0"]["InternalEnergy"]
Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1]
Radius = np.sqrt(Radius)
H=H0*np.power((Radius/Rin),(Temp_Power+3.)/2)/m.sqrt(betaz/2.)
sin = pos[:,1]/Radius
cos = pos[:,0]/Radius
theta= np.arcsin(sin)
for i in range(0,len(Radius)):
    if (cos[i]<0):
        theta[i] = np.pi - theta[i]
    if (cos[i]>0) and (sin[i]<0):
        theta[i] = np.pi * 2 + theta[i]
if case =="pertvel":
    del s["PartType0"]["Velocities"]
    Vr, Vphi= Cartisian2Polar(vel[:,0],vel[:,1],cos,sin)
    Vr=ample* C_s* np.power((Radius/Rin),Temp_Power/2.) *np.cos(kr*Radius+mphi*theta+kz*pos[:,2])*kz/kr
    vel[:,2] = ample* C_s* np.power((Radius/Rin),Temp_Power/2.) *np.cos(kr*Radius+mphi*theta+kz*pos[:,2]) #np.random.uniform(-1.,1.,len(mass)) 
    vel[:,0], vel[:,1] = Polar2Cartisian(Vr,Vphi,cos,sin)
    s.create_dataset("PartType0/Velocities",data=vel)

if case=="pertrho":
    del s["PartType0"]["Coordinates"]
    pos[:,0] -= ample *sin[:]*np.cos(kr*Radius+mphi*theta+kz*pos[:,2]) *Radius
    pos[:,1] += ample *cos[:]*np.cos(kr*Radius+mphi*theta+kz*pos[:,2]) *Radius
    s.create_dataset("PartType0/Coordinates",data=pos)

s.flush()
s.close()
            

