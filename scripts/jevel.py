# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                     \
                                                                                                                                                                                                            
import pynbody
import pynbody.plot.sph as sph
import matplotlib.pylab as pl
from pynbody.filt import *
startnum=100
step=20
numbsnaps=221
fname_base="snapshot_"
fname_ext=".hdf5"
rbins=np.arange(2,50,0.5)
colors=plt.cm.cool(np.linspace(0,1,7))
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        outname=fname.replace(".hdf5",".png")

        jtheta=[]
        jphi=[]
        jtot=[]

        pynbody.analysis.angmom.sideon(s.gas)
        pynbody.analysis.halo.center(s,mode="hyb")
        s.g["vel"]-=s.s["vel"]
        for j in range (1,len(rbins)):
            annul=s.g[Annulus(rbins[j-1],rbins[j])]
            jv=pynbody.analysis.angmom.ang_mom_vec(annul)/np.sum((annul.g["mass"]))
            jj=np.sqrt(jv[0]*jv[0]+jv[1]*jv[1]+jv[2]*jv[2])
            jtot.append(jj)
            jtheta.append(math.acos(jv[2]/jj)*180./3.1415)
            jphi.append(math.atan2(jv[1],jv[0])*180./3.1415)

        pynbody.analysis.angmom.faceon(s.gas)
        disk=s[Disc(50,20)]
        pynbody.analysis.angmom.faceon(disk.gas)
        pynbody.analysis.halo.center(disk,mode="hyb")
        p=pynbody.analysis.profile.Profile(disk.g,min=2,max=50)

        ax1=plt.subplot(4,1,1)
        ax1.plot(rbins[1:],jtot,color=colors[(i-100)/20])
        
        ax2=plt.subplot(4,1,2)
        ax2.plot(rbins[1:],jtheta,color=colors[(i-100)/20])
        
        ax3=plt.subplot(4,1,3)
        ax3.plot(rbins[1:],jphi,color=colors[(i-100)/20])
        
        ax4=plt.subplot(4,1,4)
        ax4.plot(p["rbins"],p["density"],color=colors[(i-100)/20])
        print i
plt.show()