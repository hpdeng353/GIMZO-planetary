import numpy as np
import matplotlib.pyplot as plt
import math as m
Num_Sample=30000
Rho_Power=-1 # negetive
Temp_Power=-1
soft =2
Rin=0.5
Rout=10.0
# calculate the normalized accumulative distribution function and use findloc solveline to get radius
def SIGMA_profile(num_sample,rho_power, temp_power,rin,rout):
    rr=np.arange(0,1.,1./num_sample) * (rout - rin) + rin
    rho0=np.zeros(num_sample)
#    powerfile=np.zeros(num_sample)
    temp=np.zeros(num_sample)
 #   powerfile1=np.zeros(num_sample)
    step = 1./num_sample * (rout -rin)
    normal=0
#    normal1=0
    probability= np.zeros(num_sample)

    for i in range(0,len(rr)):
        temp[i] = m.pow(rr[i]/rin )
        rho0[i] = (m.sqrt(rout/rr[i]) + m.sqrt(rin/rr[i]) - m.sqrt(rout * rin / rr[i] / rr[i]) -1) * m.pow((rr[i]/rin),power) * rr[i]/rin #note the final term corespond to the geometry factor 2*pi*r*rho
        powerfile[i] = m.pow((rr[i]/rin),power)* rr[i]/rin
#            powerfile1[i] = m.pow((rr[i]/rin),power)
#            if powerfile1[i] > 25 :
#                powerfile1[i] = 0

        normal += step * SIGMA[i]
        normal1 += step * powerfile[i]
        if i < len(rr) - 1 :
            probability[i+1] = probability[i] + step * SIGMA[i]
    probability = probability / normal # accumulative probability distribution
    SIGMA /= normal
    powerfile /= normal1
    return rr, SIGMA, powerfile, probability

#use montel 
    
def Cartisian2polar(pos):
    pos[0] = m.sqrt(pos[0]*pos[0] + pos[1] * pos[1] + pos[2] * pos[2])
    pos[1] = m.atan(pos[1], pos[0]) # between -pi to pi

def Polar2Cartisian(radius,theta):
    xx=np.zeros(len(radius))
    yy = np.zeros(len(radius))
    for i in range(0,len(radius)):
        xx[i] = radius[i] * m.cos(theta[i])
        yy[i] = radius[i] * m.sin(theta[i])
    return xx, yy

def Testpowerfile(num_sample,power,rin,rout):
    radius = np.random.uniform(0,1.0,num_sample)
    theta = np.random.uniform(-1.0,1.0,num_sample) * m.pi 
    for i in range(0, len(radius)):
        radius[i] = m.pow(rout/rin, radius[i]) * rin
    return radius, theta

def Testrand(num_sample, power, rin, rout):
    rr = np.arange(0,1.0,1.0/num_sample) * ( rout - rin) + rin
    step = 1.0/num_sample
    rho=np.zeros(num_sample)
    normal=0
    prob = np.zeros(num_sample)
    for i in range(0,len(rr)):
            rho[i] = m.pow((rr[i]/rin),power)#(m.sqrt(rout/rr[i]) + m.sqrt(rin/rr[i]) - m.sqrt(rout * rin / rr[i] / rr[i]) -1) * m.pow((rr[i]/rin),power)
            normal += step * rho[i]
            
    prob = rho / normal
    radius = np.random.choice(rr,prob) #this doesn't work because array is too large
    theta = np.random.uniform(-1.0,1.0,num_sample) * m.pi 
    return radius, theta

def solveline(x1,x2,y1,y2,y):
    k = (y2-y1)/(x2-x1)
    x = (y-y1)/k + x1
    return x

def findloc(arr,val):
    j=len(arr)
    i=0
    trial= j/2
    while True :
        if m.fabs(j-i) == 1: 
            break
        if arr[trial] < val :
            i = trial
            trial = (i+j)/2
        else:
            j = trial
            trial = (i+j)/2
    return i,j

def findpos(num_sample,probability,rin,rout):
    uniformarr=np.random.uniform(0,1.0,num_sample)
    rr = np.arange(0,1.0,1.0/num_sample) * ( rout - rin) + rin
    radius =np.zeros(num_sample)
    for i in range(0,len(uniformarr)):
        index1,index2=findloc(probability,uniformarr[i])
        radius[i] = solveline(rr[index1],rr[index2],probability[index1],probability[index2],uniformarr[i])
#        if i<10:
 #           print uniformarr[i],radius[i],index1,index2
    theta = np.random.uniform(-1.0,1.0,num_sample) * m.pi 
    return radius, theta

    



Radius, Rho, Powerfile, Probability = SIGMA_profile(Num_Sample, Rho_Power, Temp_Power, Rin, Rout)
#Radius,Theta = Testrand(Num_Sample, Power, Rin, Rout)
#Radius,Theta = findpos(Num_Sample,Probability,Rin,Rout)
#XX, YY = Polar2Cartisian(Radius,Theta)

plt.figure()
plt.plot(Radius,Powerfile)
plt.plot(Radius,Rho)
#plt.scatter(XX, YY, s=1)
plt.axis('equal')
plt.title('Rin= %.2f Rout= %.2f Power= %.2f' %(Rin, Rout, Power))
plt.show()

        
        
    
    
