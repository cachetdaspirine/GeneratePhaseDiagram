import numpy as np
import System as Sys
import Shape as Sh
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def line(x,a,b):
    return a*x+b
class BF:
    def __init__(self,WidthMax,P):
        #Make an array of system for each W
        self.Wmax = WidthMax
        self.Aspect = np.array([1./3.,1./5.,0.1])
        self.Systems = np.array([np.empty(3,dtype=object)#np.array([
        #Sys.System(
        #Sh.Fiber(w,int(w/A),ParticleType=P.ParticleType),
        #eps=P.epsilon,
        #Kcoupling=P.kc,
        #Kvol=P.kA,
        #ParticleType=P.ParticleType)
        #for A in self.Aspect])
        for w in range(0,WidthMax,1)])
    def CheckInfFiber(self,w,P):
        E=list()
        if w >= self.Wmax:
            w = self.Wmax
        if not self.Systems[w-1][0]:
            for n,a in enumerate(self.Aspect):
                self.Systems[w-1][n] = Sys.System(
                Sh.Fiber(w,int(w/a),ParticleType=P.ParticleType),
                eps=P.epsilon,
                Kcoupling=P.kc,
                Kvol=P.kA,
                ParticleType=P.ParticleType)
        for i,S in enumerate(self.Systems[w-1]):
            A = self.Aspect[i]
            #loop over several aspect ratio
            FiberArray = Sh.Fiber(w,int(w/A),ParticleType=P.ParticleType)
            SurfaceEnergy = Sh.SurfaceEnergy(FiberArray,J=P.J,ParticleType=P.ParticleType)
            E.append((S.Energy+SurfaceEnergy)/(w*int(w/A)))

        popt,pconv = curve_fit(line,self.Aspect,E)
        print('a, b = '+str(popt))
        fig,ax = plt.subplots()
        ax.plot(self.Aspect,line(self.Aspect,popt[0],popt[1]))
        ax.scatter(self.Aspect,E)
        return fig, ax  
    def Get_Einf(self,w,P):
        E=list()
        if w >= self.Wmax:
            w = self.Wmax
        if not self.Systems[w-1][0]:
            for n,a in enumerate(self.Aspect):
                self.Systems[w-1][n] = Sys.System(
                Sh.Fiber(w,int(w/a),ParticleType=P.ParticleType),
                eps=P.epsilon,
                Kcoupling=P.kc,
                Kvol=P.kA,
                ParticleType=P.ParticleType)
        for i,S in enumerate(self.Systems[w-1]):
            A = self.Aspect[i]
            #loop over several aspect ratio
            FiberArray = Sh.Fiber(w,int(w/A),ParticleType=P.ParticleType)
            SurfaceEnergy = Sh.SurfaceEnergy(FiberArray,J=P.J,ParticleType=P.ParticleType)
            E.append((S.Energy+SurfaceEnergy)/(w*int(w/A)))
        try :
            popt,pconv = curve_fit(line,self.Aspect,E)
        except ValueError:
            for s in self.Systems[w-1]:
                print(s.Energy)
            print(w)
            print(E)
            raise
        return popt[1]

    def Get_Best_Fiber(self,P):
        # given a set of parameter, return the width and the energy of the best
        # fiber for this set of parameter
        w=1
        Einf1 = self.Get_Einf(w,P)
        w+=1
        Einf2 = self.Get_Einf(w,P)
        while Einf2<Einf1 :
            w+=1
            Einf1=Einf2
            Einf2 = self.Get_Einf(w,P)
            if w>=self.Wmax-1:
                return w,Einf2
            #while we haven't found the best fiber
        return w-1,Einf1
