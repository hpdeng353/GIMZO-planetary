# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                      
import pynbody
import pynbody.plot.sph as sph
startnum=0
step=1
numbsnaps=98
fname_base="snapshot_"
fname_ext=".hdf5"
temp=[]





for i in range(startnum,numbsnaps):
    if (i % step ==0):
        if (i<10): fname=fname_base+'00'+str(i)+fname_ext     
        if (i>=10): fname=fname_base + '0'+str(i) + fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        j=i+1
        if (j<10): fname1=fname_base+'00'+str(j)+fname_ext
        if (j>=10): fname1=fname_base+'0'+str(j)+fname_ext
        if (j>=100): fname1=fname_base + str(j) + fname_ext

        s=pynbody.load(fname)
        pos=s.g["pos"]-s.s["pos"]                
        rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1])
        
        s1=pynbody.load(fname1)
        print i 



        ids1=s.g["iord"]
        ids2=s1.g["iord"]
        mask=np.isin(ids1,ids2,invert=True)

        rad1=rad[mask]
        temp.append(len(rad1[np.where((rad1>10)&(rad1<20))]))
        
