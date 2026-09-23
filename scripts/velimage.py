# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
#matplotlib.use('Agg')                                                                                                                                                                                      
import pynbody
import pynbody.plot.sph as sph
startnum=49
step=1
numbsnaps=50
fname_base="snapshot_"
fname_ext=".hdf5"
temp=[]
eps=[]
        
s=pynbody.load('snapshot_080.hdf5')
pynbody.analysis.halo.center(s,mode="com")
s1=s[pynbody.filt.Cuboid(10,-1,-1,16,5,1)]
pynbody.analysis.halo.center(s1,mode="hyb")
s2=s1[pynbody.filt.Cuboid(-2,-2,-1,2,2,1)]
pynbody.analysis.halo.center(s2,mode="hyb")
cmap=plt.get_cmap('seismic',8)
sph.velocity_image(s2.g,qty="rho",width=4,cmap=cmap,log=True,resolution=1000,units='g cm^-3',mode='stream', show_cbar=True,vmin=1e-12,vmax=1.e-8,vector_color='C2',vector_resolution=100,density=2)


