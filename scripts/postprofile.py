# coding: utf-8

get_ipython().magic(u'run postimpact2.py')
binentropy*=1e6
plt.plot(bins,binentropy,label="Hit-and-run",color="#fb0a2a")
plt.tick_params(axis='both',labelsize=16)
plt.xlabel(r"$M/M_{\oplus}$",fontsize=16)
plt.ylabel(r"Entropy (J/kg/K)",fontsize=16)

a=np.arange(0,0.3,0.01)
b=np.zeros(len(a))
b[:]=1200
plt.plot(a,b,"k--")
a=np.arange(0.3,1.,0.01)
b=np.zeros(len(a))
b[:]=2700
plt.plot(a,b,"k--",label="Initial Value")
plt.legend(loc="upper left",fontsize=16)


get_ipython().magic(u'run postimpact2.py')

plt.figure(figsize=(8,6))
plt.plot(bins,binmixing,label="Canonical",color="#02adea")
get_ipython().magic(u'run postimpact2.py')
plt.plot(bins,binmixing,label="Hit-and-run",color="#fb0a2a")
plt.tick_params(axis='both',labelsize=16)
plt.xlabel(r"$M/M_{\oplus}$",fontsize=16)
plt.ylabel(r"$F_{tar}$",fontsize=16)
plt.legend(fontsize=16)
