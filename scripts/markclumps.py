# coding: utf-8
import pynbody
import diskpy
s=pynbody.load("/scratch/snx3000/hpdeng/gizmout/grvdisk3/snapshot_133.hdf5")
clump=diskpy.clumps.clumpfinding.find_clumps(s)
a=clump[np.where(clump>0.5)]
len(a)/1193
len(a)/1193.
diskpy.clumps.clumpfinding.clump_im(s,clump_array=clump,width="30 au")
