import numpy as np
import matplotlib.pyplot as plt
import math as m
Rin=0.01
Rout=5
Rho_Power = -3
X0=6.*Rin
#Sharp=0.05*X0
in_over_out = 0.1
#Sharp = (Rin-X0)/np.log(1./(1.0-in_over_out*m.pow((Rout/Rin),Rho_Power))-1.)
A=10
P=8
Soft=0.1
Ms=10
#print Sharp
def taper1(x,x0,sharp): #default setting center of taper x0=3*rin
    f = 1.0 - 1.0/(m.exp((x-x0)/sharp)+1.)
    return f

def taper2(x,rin,rout):
    f = m.sqrt(rin/x) + m.sqrt(rout/x)-m.sqrt(rin*rout/x/x)-1.
    return f
    
    
def dr_taper1(x,x0,sharp=0.1):
    dfdr = 1.0/(m.exp((x-x0)/sharp)+1.)/(m.exp((x-x0)/sharp)+1.)/sharp*m.exp((x-x0)/sharp)
    return dfdr

def dr_taper2(x,rin,rout):
    dfdr = -0.5*(m.sqrt(rin)+m.sqrt(rout))*m.pow(x,-1.5)+m.sqrt(rin*rout)/x/x
    return dfdr

def soft1(x,soft,ms,a,p):
    f = ms/x/x*m.exp(-a*m.pow(soft/x,p))
    return f
def soft(x,soft,ms):
    f= ms/(x*x+soft*soft)
    return f
    
rr=np.arange(Rin,Rout,0.0001)
shape=np.zeros(len(rr))
dshape=np.zeros(len(rr))
rhoprofile=np.zeros(len(rr))
force1=np.zeros(len(rr))
force=np.zeros(len(rr))

for i in range(0,len(rr)):
 #   shape[i] = taper1(rr[i],X0,Sharp)
  #  shape[i] = taper2(rr[i], Rin, Rout)
   # dshape[i] = dr_taper1(rr[i],X0,Sharp)*rr[i]/shape[i]
#    dshape[i] = dr_taper2(rr[i], Rin, Rout)*rr[i]/shape[i]
    #rhoprofile[i]=m.pow((rr[i]/Rin),Rho_Power)*shape[i]
    force[i] = soft(rr[i],Soft,Ms)
    force1[i] = soft1(rr[i],Soft,Ms,A,P)
#rhofloor=rhoprofile[0]*0.99
#for i in range(0,len(rr)):
#    rhoprofile[i] -= rhofloor

plt.figure()
plt.plot(rr,force)
plt.plot(rr,force1,'-r')
plt.title("soft shape")

"""
plt.figure()
plt.plot(rr,dshape)
plt.title("r*dfdr/f")

plt.figure()
plt.plot(rr,rhoprofile)
plt.title("middle plane rho profile")
"""
plt.show()