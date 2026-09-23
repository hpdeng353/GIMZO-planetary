# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                                                       
import pynbody
import pynbody.plot.sph as sph
startnum=0
step=1
numbsnaps=91
fname_base="snapshot_"
fname_ext=".hdf5"


for i in range(startnum,numbsnaps):
    if (i % step ==0):
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        bins=np.arange(1,50,0.1)
        rho=s.g['rho']
        vel=s.g['vel']
        pos=s.g["pos"]
        rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
        lx=[]
        ly=[]
        for k in range(0,len(bins)-1):
            vel1=vel[np.where((rad>bins[k])&(rad<bins[k+1]))]
            pos1=pos[np.where((rad>bins[k])&(rad<bins[k+1]))]
            jj=np.cross(pos1,vel1)
            rho1=rho[np.where((rad>bins[k])&(rad<bins[k+1]))]
            j1=jj[np.where(rho1>np.max(rho1)*0.9)]
            j2=np.average(j1,axis=0)

            lx.append(j2[0]/np.sqrt(j2.dot(j2)))
            ly.append(j2[1]/np.sqrt(j2.dot(j2)))
        plt.plot(bins[:-1],lx,label=str(i))
        plt.plot(bins[:-1],ly,label=str(i),linestyle='--')
        print i
        
