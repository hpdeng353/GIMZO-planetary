# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                          
import pynbody
import pynbody.plot.sph as sph
import pynbody.filt as filt
jx=[]
jy=[]
jz=[]
startnum=0
step=1
numbsnaps=660
fname_base="/scratch/snx3000/hpdeng/gizmout/grvdisk1/snapshot_"
fname_ext=".hdf5"
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        outname=fname.replace(".hdf5",".png")
        pynbody.analysis.halo.center(s,mode="hyb")     
        pynbody.analysis.angmom.faceon(s.gas)
        h=s[filt.Sphere(20)]
        j=pynbody.analysis.angmom.ang_mom_vec(h.g)
        jx.append(j[0])
        jy.append(j[1])
        jz.append(j[2])
        print i
        
