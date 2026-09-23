# coding: utf-8
points=np.array([pos1[:,0],pos1[:,1],pos1[:,2]]).T.reshape(-1,1,3)
seg=np.concatenate([points[:-1],points[1:]],axis=1)
norm=plt.Normalize(mag.min(),mag.max())

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
lc=Line3DCollection(seg,cmap="copper",norm=norm)
line=ax.add_collection(lc)

lc.set_linewidth(2)
lc.set_array(mag)
fig.colorbar(line,ax=ax)

ax.set_zlim3d(-1,1)
ax.set_ylim3d(-30,30)
ax.set_xlim3d(-30,30)
