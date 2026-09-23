# coding: utf-8
# %load lx.py
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                     \
                                                                                                                                                                                                            
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pylab as pl
from pynbody.filt import *
startnum=00
step=40
numbsnaps=221

fname_base="snapshot_"
fname_ext=".hdf5"
rr=np.logspace(0,np.log10(40),num=100)
colors=plt.cm.cool(np.linspace(0,1,10))
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        pynbody.analysis.halo.center(s,'com')
        pos=s.g['pos']
        vel=s.g['vel']
        jvec=np.cross(pos,vel)
        rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])

        jtot=[]
        lx=[]
        ly=[]
        for l in range (1,len(rr)):
            jv=np.average(jvec[np.where((rad>rr[l-1])&(rad<rr[l]))],axis=0)
            jj=np.sqrt(jv[0]*jv[0]+jv[1]*jv[1]+jv[2]*jv[2])
            jtot.append(jj)
            lx.append(jv[0]/jj)
            ly.append(jv[1]/jj)
        plt.plot(rr[:-1],lx,c=colors[i//step])
        
        
        
