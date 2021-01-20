import numpy as np
import System as Sys
import Shape as Sh

class BD:
    def __init__(self,Nmax,P):
        self.Nmax = Nmax
        self.Range = P.HRange(Nmax)
        self.Systems = np.empty(self.Range.shape[0],dtype=object)
        #self.Systems = np.array([
        #Sys.System(
        #Sh.Parallel(r,ParticleType=P.ParticleType),
        #eps=P.epsilon,
        #Kcoupling=P.kc,
        #Kvol=P.kA,
        #ParticleType=P.ParticleType)
        #for r in self.Range])
    def Get_E(self,n,P):
        if n>= self.Range.shape[0]:
            n=-1
        DiskArray = Sh.Parallel(self.Range[n],ParticleType=P.ParticleType)
        if self.Systems[n]:
            return (self.Systems[n].Energy+Sh.SurfaceEnergy(DiskArray,J=P.J,ParticleType=P.ParticleType))/Sh.Np(DiskArray)
        else :
            self.Systems[n] =Sys.System(
            DiskArray,
            eps=P.epsilon,
            Kcoupling=P.kc,
            Kvol=P.kA,
            ParticleType=P.ParticleType)
            return (self.Systems[n].Energy+Sh.SurfaceEnergy(DiskArray,J=P.J,ParticleType=P.ParticleType))/Sh.Np(DiskArray)
    def Get_Best_Disk(self,P):
        n1 = 0
        E1 = self.Get_E(n1,P)
        n2 = 1
        E2 = self.Get_E(n2,P)
        while E2<E1:
            n1=n2
            n2+=1
            E1=E2
            E2 = self.Get_E(n2,P)
            if n2==self.Range.shape[0]-1:
                return Sh.Np(Sh.Parallel(self.Range[n2],ParticleType=P.ParticleType)),E2
        return Sh.Np(Sh.Parallel(self.Range[n1],ParticleType=P.ParticleType)),E1
