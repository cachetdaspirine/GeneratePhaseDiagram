from Generate import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
import matplotlib
matplotlib.use('pdf')

# Data for phase diagram generation
NAME = 'H_L7_5_FULL'
numin = 0.34
numax = 0.99
Gammamin = 0.
Gammamax = 1.2
NpointsNu = 10
NpointsGamma=50
Nmax = 4000
Wmax = 40
OrderMax = 0

L = 7.5
PTYPE = 'Hexagon'
EPS = 0.01
G = Generate(L,EPS,PTYPE)
Gamma,nu,Color = G.MakePhaseDiagram(numin,numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,Wmax,OrderMax)
#Gamma,nu,Color = np.loadtxt('Gamma.txt',dtype=float),np.loadtxt('nu.txt',dtype=float),np.loadtxt('Color.txt',dtype=float)
np.save('Res/G_'+NAME,Gamma)
np.save('Res/nu_'+NAME,nu)
np.save('Res/C_'+NAME,Color)
#print(Color)

# Create the figure, and save it

#BulkColor = (86./255,42./255.,132./255)
cmaplist = [(86./255,42./255.,132./255.,1.), (156./255,18/255.,109/255.,1.),(204./255.,35./255.,129./255.),(0,0,0,1)]
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, 3)
fig,ax = plt.subplots(figsize=(16,10))
#cmap = cm.Reds#.reversed()
#cmap.set_bad(color=BulkColor)
#masked=np.ma.masked_array(Color,Color>=0)
#psm = ax.pcolormesh(Gamma,nu,-masked,cmap=cm.Reds)#,norm=mcolors.LogNorm())
#psm2=ax.pcolormesh(Gamma,nu,np.ma.masked_array(Color,Color<=0),cmap=cm.Blues)

#define the hexagon region :
HexRegion = np.ma.masked_array(Color[:,:,0],Color[:,:,0]<=0)
psmHex = ax.pcolormesh(Gamma,nu,HexRegion,cmap=cm.Reds)#,norm=mcolors.LogNorm())
cba=plt.colorbar(psmHex)
#Define the Fiber Region :
FiberRegion1 = np.ma.masked_array(Color[:,:,1],Color[:,:,1]<=0)
psmFiber1 = ax.pcolormesh(Gamma,nu,FiberRegion1,cmap=cm.Blues)
cbb1=plt.colorbar(psmFiber1)
FiberRegion2 = np.ma.masked_array(Color[:,:,1],Color[:,:,1]>=0)
psmFiber2 = ax.pcolormesh(Gamma,nu,FiberRegion2,cmap=cm.Greens)
cbb2=plt.colorbar(psmFiber2)
#Define the different bulk regions :
BulkRegion = np.ma.masked_array(Color[:,:,2],Color[:,:,2]<=0)
psmBulk = ax.pcolormesh(Gamma,nu,BulkRegion,cmap = cmap)
#cbc = plt.colorbar(psmBulk)
#psm = ax.pcolormesh(J,Ka,-masked,cmap=cmap)#,norm=colors.LogNorm())
#psm2=ax.pcolormesh(J,Ka,np.ma.masked_array(Width,Width<=0),cmap=cm.Blues)
#print(Color[:,:,0])
#print(Color[:,:,1])
#print(Color[:,:,2])

cba.set_label('Disk Radius',fontsize=20)
cbb1.set_label('fiber_1 width',fontsize=20)
cbb2.set_label('fiber_2 width',fontsize=20)
plt.xlabel('$\Gamma$ : rescaled surface tension',fontsize=20)
plt.ylabel('$\\nu$ : Poisson\'s ratio',fontsize=20)
plt.title('Phase diagram of the analytic applied to finite aggregate',fontsize=20)
plt.savefig('Res/'+NAME+'.pdf')
