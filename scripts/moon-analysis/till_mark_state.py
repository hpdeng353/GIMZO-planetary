"""
This script reads a tipsy binary and writes the state (condensed, expansed: cold, intermediate or hot).
"""
from matplotlib import *
from numpy import *
from matplotlib.pyplot import *
from sys import *

# Used to read tipsy binary files
import pynbody

# For the colormaps
from brewer2mpl import *

# Set a font
rcParams['font.family'] = 'serif'
#rcParams['font.size'] = 14.0

if len(argv) != 2:
		print "Usage: till_mark_state.py tipsyfile.std"
		exit(1)

tipsyfile = argv[1]

"""
We do not convert code units to cgs here because L=1RE is a good
unit to visualize a planet. But we want to convert the density
to cgs.
"""

bin = pynbody.load(tipsyfile)

pos = bin.gas["pos"]
x = pos[:,0]
y = pos[:,1]
z = pos[:,2]
rho = bin.gas["rho"]
u = bin.gas["temp"]
iMat = bin.gas["metals"]

# The cold curve for basalt and iron
cold_granite = loadtxt('cold_basalt.txt')
cold_iron = loadtxt('cold_iron.txt')

# The region, where P<0
#np_granite = loadtxt('pressneg.basalt.fANEOS.txt')
#np_iron = loadtxt('pressneg.iron.txt')

"""
rhop_granite = np_granite[:,0]
up_granite = np_granite[:,1]
rhop_iron = np_iron[:,0]
up_iron = np_iron[:,1]
"""

# Load cold curve
rhocold1 = cold_granite[:,0]
ucold1 = cold_granite[:,1]

rhocold2 = cold_iron[:,0]
ucold2 = cold_iron[:,1]

"""
# Unit convertion factors
ErgPerGmUnit = 9998228982.69		# erg/g
GmPerCcUnit = 0.368477421278		# g/cc
SecUnit = 6378.6921465				# s
Lunit = 637812728.278				# cm
Munit = 9.5607162e+25				# g

# Convert to cgs
u *= ErgPerGmUnit
rho *= GmPerCcUnit
"""

# Values for iron
us_core = 1.42025
us2_core = 8.4515
rho0_core = 21.331

# Values for granite
us_mantle = 3.5
us2_mantle = 18.0
rho0_mantle = 7.33

"""
# Values for basalt
us_mantle = 4.72
us2_mantle = 487.0
rho0_mantle = 7.33
"""
i = where(iMat == 1)
j = where(iMat == 0)

#figure(figsize=(13,10))

"""
First plot the iron particles.
"""
"""
subplot(1,2,1)
xmax = max(rho[i])*1.1
ymax = max(u[i])*1.1

xlim(0,xmax)
ylim(0,ymax)

title(r'Iron')
xlabel('Density')
ylabel('Internal energy')

plot(rhocold2,ucold2,'-',color='red',linewidth=2,label='Cold curve (Iron)')
fill_between(rhop_iron,up_iron,color='orange', linewidth=0,alpha=0.5)
scatter(rho[i],u[i],s=16,c='blue',linewidth=0.01)

# Core
plot([0,rho0_core],[us_core,us_core],'r--')
plot([0,rho0_core],[us2_core,us2_core],'r--')
plot([rho0_core,rho0_core],[0,ymax],'r--')
"""
"""
# We need 300 dpi for the small format and 150 for the large one
savefig(tipsyfile+'.iron.png', dpi=150, bbox_inches='tight')
close()
"""

"""
Then the basalt particles.
"""
subplot(1,2,2)
xmax = max(rho[j])*1.1
ymax = max(u[j])*1.1

xlim(0,xmax)
ylim(0,ymax)

title(r'Basalt')
xlabel('Density')
#ylabel('Internal energy')

plot(rhocold1,ucold1,'-',color='red',linewidth=2,label='Cold curve (Granite)')
#fill_between(rhop_granite,up_granite,color='orange', linewidth=0,alpha=0.5)
scatter(rho[j],u[j],s=16,c='green',linewidth=0.01)

# Mantle
plot([0,rho0_mantle],[us_mantle,us_mantle],'r--')
plot([0,rho0_mantle],[us2_mantle,us2_mantle],'r--')
plot([rho0_mantle,rho0_mantle],[0,ymax],'r--')

#show()
"""
# We need 300 dpi for the small format and 150 for the large one
savefig(tipsyfile+'.basalt.png', dpi=150, bbox_inches='tight')
"""

# We need 300 dpi for the small format and 150 for the large one
savefig(tipsyfile+'.composition.png', dpi=150, bbox_inches='tight')

close()

print "Done."
