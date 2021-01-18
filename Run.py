from Generate import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
import matplotlib
matplotlib.use('pdf')

# Data for phase diagram generation
NAME = 'PD_Triangle_L_5'
numin = 0.34
numax = 0.99
Gammamin = 0.
Gammamax = 1.25
NpointsNu = 3
NpointsGamma=3
Nmax = 1500
Wmax = 30

L = 5.
PTYPE = 'Triangle'
EPS = 0.01
G = Generate(L,EPS,PTYPE)
#Gamma,nu,Color = G.MakePhaseDiagram(numin,numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,Wmax)
Gamma,nu,Color = np.loadtxt('Gamma.txt',dtype=float),np.loadtxt('nu.txt',dtype=float),np.loadtxt('Color.txt',dtype=float)
#np.savetxt('Gamma.txt',Gamma)
#np.savetxt('nu.txt',nu)
#np.savetxt('Color.txt',Color)
print(Color)

# Create the figure, and save it

BulkColor = (86./255,42./255.,132./255)
fig,ax = plt.subplots(figsize=(16,10))
cmap = cm.Reds#.reversed()
cmap.set_bad(color=BulkColor)
masked=np.ma.masked_array(Color,Color>=0)
psm = ax.pcolormesh(Gamma,nu,-masked,cmap=cm.Reds,norm=mcolors.LogNorm())
psm2=ax.pcolormesh(Gamma,nu,np.ma.masked_array(Color,Color<=0),cmap=cm.Blues)

#psm = ax.pcolormesh(J,Ka,-masked,cmap=cmap)#,norm=colors.LogNorm())
#psm2=ax.pcolormesh(J,Ka,np.ma.masked_array(Width,Width<=0),cmap=cm.Blues)

cba=plt.colorbar(psm)
cbb=plt.colorbar(psm2)
cba.set_label('Disk Radius',fontsize=20)
cbb.set_label('fiber width',fontsize=20)
plt.xlabel('$\Gamma$ : rescaled surface tension',fontsize=20)
plt.ylabel('$\\nu$ : Poisson\'s ratio',fontsize=20)
plt.title('Phase diagram of the analytic applied to finite aggregate',fontsize=20)
plt.savefig(NAME+'.pdf')
