# coding: utf-8
# %load halovis.py
#import yt
fname = 'snapshot_000.hdf5'

unit_base = {'UnitLength_in_cm'         : 3.08568e+18,
             'UnitMass_in_g'            :   1.989e+33,
             'UnitVelocity_in_cm_per_s' :      6559.14}

bbox_lim = 1e3 #kpc

bbox = [[-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim]]
 
ds = yt.load(fname,unit_base=unit_base,bounding_box=bbox)
ds.index
ad= ds.all_data()
sorted(ds.field_list)
#hsml= ad[("PartType1","AGS-Softening")]
mass= ad[("PartType1","Masses")]
#rho1=mass/hsml/hsml/hsml*47.75
pos= ad[("PartType1","Coordinates")]
x0=np.average(pos[:,0])
pos[:,0]-=x0
y0=np.average(pos[:,1])
pos[:,1]-=y0
z0=np.average(pos[:,2])
pos[:,2]-=z0
rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
#plt.scatter(rad,rho1,s=0.05)

rbins=np.linspace(0,250,200)
dr=rbins[1]-rbins[0]
rho=[]
for i in range(0,len(rbins)-1):
    temp=np.sum(mass[np.where((rad>rbins[i])&(rad<rbins[i+1]))])
    temp=temp/4/3.1415/rbins[i]/rbins[i]/dr
    rho.append(temp)
    
    
rho0=1000/4./3.1415/16.31**3/(np.log(16.)-15/16.)
#plt.plot(rbins,rho0*16.31/rbins/(1+rbins/16.31)**2)
plt.plot(rbins[1:],np.log(rho))
rad=np.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
zz=pos[:,2]
cosi=zz/rad
cf=(1-cosi)/2
rad1=rad[np.where(rad<50)]
cf1=cf[np.where(rad<50)]
