# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')
from pynbody.filt import * 
import pynbody
import pynbody.plot.sph as sph
startnum=170
step=10
numbsnaps=201
fname_base="snapshot_"
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        outname=fname.replace(".hdf5",".png")
        pynbody.analysis.angmom.faceon(s.gas)
        pynbody.analysis.halo.center(s,mode="hyb")



#        mag=s.g["MF"]                                                                                                                                                                                              
#        s.g["divB"]=s.g["vphi"]                                                                                                                                                                                    
#        sph.image(s.g,qty="rho",width=60,cmap="bwr",log=True,resolution=1000,vmin=1e-5,vmax=0.1)
        print i
        

        disk=s[Disc(80,20)]
        pynbody.analysis.angmom.faceon(disk)
        pynbody.analysis.halo.center(disk,mode="hyb")

        p=pynbody.analysis.profile.Profile(disk.g,min=2,max=30,ndim=2)
        cs=sqrt(1.4*0.0195463* pow(p["density"], 0.4))
        kappa=np.zeros(len(vphi))
        vphi=p["vphi"]
        vphir=vphi*p["rbins"]
        kappa=np.zeros(len(vphi))
        for j in range (0,len(vphi)-1):
            kappa[j]=(vphir[j+1]-vphir[j])/(p["rbins"][j+1]-p["rbins"][j])
    
        kappa*=2
        kappa*=vphi
        kappa/=p["rbins"]
        kappa/=p["rbins"]
        kappa=np.sqrt(kappa)
        p=pynbody.analysis.profile.Profile(disk.g,min=2,max=30,ndim=2)
        sigma=p["density"]

        
        QQ=cs*kappa/3.14/sigma
        plt.plot(p["rbins"],QQ)
        