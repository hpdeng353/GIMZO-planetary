# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                      
import pynbody
import pynbody.plot.sph as sph
startnum=60
step=1
numbsnaps=144
fname_base="snapshot_"
fname_ext=".hdf5"
temp=[]
posa=[]
temp1=[]
maga=[]
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        print i 
        ids=s.g["iord"]
        posall=s.g["pos"]
        zz=np.copy(posall[:,2])
        mag1=s.g["MF"]
        mag=np.sqrt(mag1[:,0]*mag1[:,0]+mag1[:,1]*mag1[:,1]+mag1[:,2]*mag1[:,2])
        temp=[]
        temp1=[]
        for j in range(0,len(idsm)):
            temp.append(zz[np.where(ids==idsm[j])][0])
            temp1.append(mag[np.where(ids==idsm[j])][0])
            
        
        posa.append(temp)
        maga.append(temp1)
        
