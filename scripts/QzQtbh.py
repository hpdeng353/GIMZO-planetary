import numpy as np
import h5py as h5py
import math as m
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pyplot as plt

def taper1(x,x0,sharp): #default setting center of taper x0=3*rin
    f = 1.0 - 1.0/(m.exp((x-x0)/sharp)+1.)
    return f
def taper2(x,rin,rout):
    f = m.sqrt(rin/x) + m.sqrt(rout/x)-m.sqrt(rin*rout/x/x)-1.
    return f
    
    
def dr_taper1(x,x0,sharp):
    dfdr = 1.0/(m.exp((x-x0)/sharp)+1.)/(m.exp((x-x0)/sharp)+1.)/sharp*m.exp((x-x0)/sharp)
    return dfdr
def dr_taper2(x,rin,rout):
    dfdr = -0.5*(m.sqrt(rin)+m.sqrt(rout))*m.pow(x,-1.5)+m.sqrt(rin*rout)/x/x
    return dfdr

Ngas= 3000000
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
beta =1
betaphi = beta*10./9
betaz=10*beta
H0 = C_s / V_k * m.sqrt((1.0+beta)/beta) * Rin    #caution h0**2 = k*T0/mu*m_proton /omega_k**2 *(1+beta)/beta 0.037
Aspect0 =H0/Rin
step = 1.0/1000 * (Rout - Rin)
rr = np.arange(Rin,Rout,step)
H = np.arange(Rin,Rout,step)
normal= 0
Hmax = H0 * m.pow((Rout/Rin),(Temp_Power+3.0)/2.0) * 4
stepz = 1.0/400*Hmax
zz = np.arange(-Hmax,Hmax,stepz)

for i in range(0,1000):
    Vt = V_k * V_k * Rin / rr[i] + 2.0 / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) + \
         (1.0 + beta) / beta * C_s * C_s * Temp_Power * m.pow((rr[i]/Rin),Temp_Power) \
         + (1.0 + beta) / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) * Rho_Power \
         +(1.0 + beta) / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) * \
         dr_taper1(rr[i],X0,Sharp)/taper1(rr[i],X0,Sharp)*rr[i]
    Vt1= V_k * V_k * Rin / rr[i] + 2.0 / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) + \
         (1.0 + beta) / beta * C_s * C_s * Temp_Power * m.pow((rr[i]/Rin),Temp_Power) \
         + (1.0 + beta) / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) * Rho_Power 
    H[i] = H0 * m.pow((rr[i]/Rin),(Temp_Power+3.0)/2.0) 
rho0 = 4.43476867871

s1=pynbody.load("/zbox/data/hpdeng/gizmout/grvmhd1/snapshot_000.hdf5")


pynbody.analysis.angmom.sideon(s1.gas)
pynbody.analysis.halo.center(s1.gas,mode="com")
sph.image(s1.g,qty="rho",width=25,cmap="Greys",log=True,resolution=3000,vmin=1e-6,vmax=8e-2)
RR,ZZ=np.meshgrid(rr,zz)
density = RR*ZZ
Qz = RR*ZZ
Qphi = RR*ZZ

for i in range(0,len(rr)):
    for j in range(0,len(zz)):
        density[j][i]=rho0*m.exp(1./Aspect0/Aspect0 * m.pow((rr[i]/Rin),(-Temp_Power-1.))*(1.0/m.sqrt(1+zz[j]*zz[j]/rr[i]/rr[i])-1)) * taper1(rr[i],X0,Sharp)*m.pow((rr[i]/Rin),Rho_Power)
        Qz[j][i] = 2*1.414*3.1415/m.sqrt(betaz) * m.pow((diskmass/Ngas/density[j][i]),-1.0/3.0) * H0 * m.pow((rr[i]/Rin),(Temp_Power+3.0)/2.0)
        Qphi[j][i] = 2*1.414* 3.1415 / m.sqrt(betaphi) * m.pow((diskmass/Ngas/density[j][i]),-1.0/3.0) * H0 * m.pow((rr[i]/Rin),(Temp_Power+3.0)/2.0)

plt.plot(rr,H,'r--')
plt.plot(rr,-H,'r--')
CS1=plt.contour(RR,ZZ,Qphi,4)
plt.clabel(CS1,inline=1,fontsize=10)
plt.title("Qt")
plt.show()
print Sharp

