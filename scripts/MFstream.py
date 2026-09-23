# coding: utf-8
import  matplotlib
import  matplotlib.pyplot as plt
import pynbody
import pynbody.plot.sph as sph
startnum=0
step=5
numbsnaps=40
fname_base="snapshot_"
fname_ext=".hdf5"

s=pynbody.load('snapshot_080.hdf5')
pynbody.analysis.halo.center(s,mode="com")
s1=s[pynbody.filt.Cuboid(10,-1,-1,16,5,1)]
pynbody.analysis.halo.center(s1,mode="hyb")
s2=s1[pynbody.filt.Cuboid(-2,-2,-1,2,2,1)]
pynbody.analysis.halo.center(s2,mode="hyb")
rho2=s2.g['rho']
u2=s2.g['u']
press=rho2*u2*(5./3-1)*5.27*1e6
MF=s2.g['MF']
me=(MF[:,0]*MF[:,0]+MF[:,1]*MF[:,1]+MF[:,2]*MF[:,2])/8/3.14
beta=press/me

s2.g['eps']=beta
s.g['vel']=s.g['MF']
sph.velocity_image(s2.g,qty="eps",width=4,cmap='viridis',log=True,resolution=1000,mode='stream', show_cbar=True,vmin=1e-1,vmax=1.e3,vector_color='C3',vector_resolution=100,density=2)
