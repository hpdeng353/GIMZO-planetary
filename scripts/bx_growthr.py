import numpy as np
import h5py as h5py
import math as m
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


alpha=[]
time=[]
startnum=0
step=1
timestep=0.5
numbsnaps=14
fname_base="/scratch/snx3000/hpdeng/gizmout/boxtest2/channel48/snapshot_"
#fname_base ="/home/ics/hpdeng/snapshot_"
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s = h5py.File(fname,"r+")
        density = s["PartType0"]["Density"][:]
        pos = s["PartType0"]["Coordinates"][:]
        mag = s["PartType0"]["MagneticField"][:]
        vel = s["PartType0"]["Velocities"][:]
        Radius = pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1]
        Radius = np.sqrt(Radius)
        alpha1 = np.log(mag[:,0]*mag[:,0])
        alpha2 = np.average(alpha1)
        time.append(i*timestep)
        alpha.append(alpha2)
        print i


fiducial = np.zeros(len(time))+0.75
alpha3=np.zeros(len(time))+0.75
for i in range(1,len(time)):
    alpha3[i]=0.5*(alpha[i]-alpha[i-1])/timestep

maxde1=0
maxde2=0

for i in range(0,11):
    if abs(alpha3[i]-0.75)/0.75 > maxde1 :
        maxde1 = abs(alpha3[i]-0.75)/0.75
for i in range(0,numbsnaps):
    if abs(alpha3[i]-0.75)/0.75 > maxde2 :
        maxde2 = abs(alpha3[i]-0.75)/0.75

print "max deviation of growth rate in linear regime:",maxde1
print "max deviation of growth rate in nonlinear regime:",maxde2
plt.figure()
plt.plot(time,alpha3)
plt.plot(time,fiducial,'-r')
plt.title("Nz32")
plt.show()


#plt.savefig('alpha')

    
            

