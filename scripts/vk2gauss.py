import numpy as np
import matplotlib.pyplot as plt
import math as m
xx= np.arange(-10,10,0.001)
yy= np.arange(-10,10,0.001)
yy[:]=0
zz= np.arange(-10,10,0.001)
zz[:]=0
aspect=0.3
R=3
H=R*aspect
for i in range(0,len(xx)):
    yy[i]=m.exp(-xx[i]*xx[i]/2./H/H)
    zz[i]=m.exp(1./aspect/aspect*(1./m.sqrt(1.+xx[i]*xx[i]/R/R)-1))

plt.figure
plt.plot(xx,yy)
plt.plot(xx,zz,'r')
plt.show()