# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                      
import pynbody
import pynbody.plot.sph as sph
startnum=115
step=1
numbsnaps=300
fname_base="snapshot_"
fname_ext=".hdf5"
temp1=[]
temp2=[]
temp3=[]
for i in range(startnum,numbsnaps):
    if (i % step ==0):
        clumps=[]
        fname=fname_base+'00'+str(i)+fname_ext
        if (i>=10): fname=fname_base+'0'+str(i)+fname_ext
        if (i>=100): fname=fname_base + str(i) + fname_ext
        s=pynbody.load(fname)
        outname=fname.replace(".hdf5",".clump")        
        pynbody.analysis.angmom.sideon(s)
        pynbody.analysis.halo.center(s,mode="com")
        pos=s.g['pos']
        vel=s.g['vel']
        ll=np.cross(pos,vel)  #caution y is in the total AM direction
        vv=np.sqrt(vel[:,0]*vel[:,0]+vel[:,1]*vel[:,1]+vel[:,2]*vel[:,2])
        rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
        ll2=ll[:,0]*ll[:,0]+ll[:,1]*ll[:,1]+ll[:,2]*ll[:,2]
        ll2=ll2.tolist()
        ll2=np.asarray(ll2)
        vv=vv.tolist()
        vv=np.asarray(vv)
        rad=rad.tolist()
        rad=np.asarray(rad)
        E=0.5*vv*vv-1/rad
        majora=-0.5/E
        
        ecc=np.sqrt(1-ll2/majora)
        cosi=ll[:,1]/np.sqrt(ll2)
        AMD=np.sqrt(majora)*(1-np.sqrt(1-ecc*ecc)*cosi)
        temp=np.average(AMD)
        print i,temp
        temp1.append(temp)
        np.savetxt('AMDH3ccc.txt',temp1)
        
