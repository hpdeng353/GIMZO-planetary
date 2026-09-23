import numpy as np
import h5py as h5py
import math as m
import pynbody
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
Sharp = (Rin-X0)/np.log(1./(1.0-m.pow((Rout/Rin),Rho_Power))-1.)
diskmass = 0.1 #mass unit = 0.1 mass of sun
#    soft =0.02
V_k = m.sqrt(10.0/0.2)*m.sqrt(0.2/Rin*StarMass/10.0) # rin =0.2code unit =1.25AU, 7.07 in code unit
C_s =0.41762*m.sqrt(T0/300.0/mu) # corespond to 300 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!isothemal sound speed
beta =1.
H0 = C_s / V_k * m.sqrt((1.0+beta)/beta) * Rin    #caution h0**2 = k*T0/mu*m_proton /omega_k**2 *(1+beta)/beta 0.037
Aspect0 =H0/Rin


step = 1.0/1000 * (Rout - Rin)
sigma = np.zeros(1000)
jprofile = np.zeros(1000)
jprofile1 = np.zeros(1000)

rr = np.arange(Rin,Rout,step)
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
    jprofile[i] = rr[i] * m.sqrt(Vt)
    jprofile1[i] = rr[i] * m.sqrt(Vt1)
    
    
    for j in range(0,800):
        normal += m.exp(1./Aspect0/Aspect0 * m.pow((rr[i]/Rin),(-Temp_Power-1.))*(1.0/m.sqrt(1+zz[j]*zz[j]/rr[i]/rr[i])-1)) * taper1(rr[i],X0,Sharp) * m.pow((rr[i]/Rin),Rho_Power)*step*stepz * 2. * np.pi * rr[i]
        sigma[i] += m.exp(1./Aspect0/Aspect0 * m.pow((rr[i]/Rin),(-Temp_Power-1.))*(1.0/m.sqrt(1+zz[j]*zz[j]/rr[i]/rr[i])-1)) * taper1(rr[i],X0,Sharp) * m.pow((rr[i]/Rin),Rho_Power)*stepz 
rho0 = diskmass /normal
print 'normal=', normal
for i in range(0,len(sigma)):
    sigma[i] *= rho0

print 'rho0',rho0   
vphi = np.zeros(1000)
vphi1 = np.zeros(1000)

for i in range(0,1000): #rhoprofile is suface density testsigma is surfdensity times geometric factor
    
    vphi[i] = V_k * V_k * Rin / rr[i] + 2.0 / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) + \
              (1.0 + beta) / beta * C_s * C_s * Temp_Power * m.pow((rr[i]/Rin),Temp_Power) \
              + (1.0 + beta) / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) * Rho_Power 
    vphi[i] = m.sqrt(vphi[i])
    #if (rr[i] < Rout*0.99) and (rr[i]>Rin*1.01) : #at outer edge /rho cause sigularity
    vphi1[i]= (1.0 + beta) / beta * C_s * C_s * m.pow((rr[i]/Rin),Temp_Power) * \
              dr_taper1(rr[i],X0,Sharp)/taper1(rr[i],X0,Sharp)*rr[i]
    vphi1[i] = m.sqrt(vphi1[i])
        
plt.figure()
plt.plot(rr,vphi)
plt.plot(rr,vphi1,'r')

s1=pynbody.load("/zbox/data/hpdeng/gizmout/ngrvmhd2/snapshot_001.hdf5")
p=pynbody.analysis.profile.Profile(s1.g, nbins=500, ndims=3,max=10)
plt.figure()
plt.plot(p["rbins"], p["density"],label='snapshot')
plt.plot(rr,sigma,'-r',label='analytical')
plt.title("surface density")
plt.legend()

plt.figure()
jtot_fixed=np.ma.fix_invalid(p["jtot"],fill_value=0)
plt.plot(p["rbins"],jtot_fixed.data)

plt.plot(rr,jprofile,'-r',label="analytical")
plt.plot(rr,jprofile1,label="non_taper_component")
plt.title("Jtot profile")
plt.legend()


plt.show()
print Sharp