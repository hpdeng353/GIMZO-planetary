# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                      
import pynbody
import pynbody.plot.sph as sph
startnum=180
step=30
numbsnaps=182
fname_base="snapshot_"
fname_ext=".hdf5"
temp1=[]
temp2=[]
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        clumps=[]
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        outname=fname.replace(".hdf5",".clump")        
        pynbody.analysis.halo.center(s,mode="com")
        fname1=outname+'.npy'

        pos=s.g['pos']
        poss=s.s['pos']        
        ids=s.g['iord']
        xx=pos[:,0]
        yy=pos[:,1]
        clump=np.load(fname1)
        clump=clump[1:]
        nclump=np.max(clump)
        match=0
        ncomm=0
        ncomm1=0
        for jj in range(1,nclump+1):        
            ids1=ids[np.where(clump==jj)]
            ncomm1=len(np.intersect1d(ids0, ids1))
            if ncomm1>ncomm:
                ncomm=ncomm1
                match=jj
        
        print i, ncomm
        xx0=np.average(xx[np.where(clump==match)])
        yy0=np.average(yy[np.where(clump==match)])
        s1=s[pynbody.filt.Cuboid(xx0-1,yy0-1,-1,xx0+1,yy0+1,1)]
        pynbody.analysis.halo.center(s1,mode="hyb")
        p = pynbody.analysis.profile.Profile(s1.g,max=1,ndim=3)
        rho1=s1.g['rho']*5.94e-7
        mass1=s1.g['mass']
        pos1=s1.g['pos']
        rad1=np.sqrt(pos1[:,0]*pos1[:,0]+pos1[:,1]*pos1[:,1]+pos1[:,2]*pos1[:,2])
        
        mass22=mass1[np.where((rad1<0.7)&(rho1>1.e-9))]
        mass12=mass1[np.where((rad1<0.7)&(rho1>1.e-10))]

        temp1.append(np.sum(mass12)/1e-3)
        temp2.append(np.sum(mass22)/1e-3)
        

        plt.plot(p['rbins'],p['density']*6e-7, label=fname)
        
