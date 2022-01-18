import numpy as np
import sys
sys.path.append("/home/hcleroy/Extra_Module_py")
sys.path.append("/home/hleroy/Simulation/Extra_Module_py")
import Conversion as Conv
from Numeric_Hex_Energy import *
from Numeric_Fiber_Energy import *
from Numeric_Bulk_Energy import *
import MeasurePoisson as MP

class Generate:
    def __init__(self,L,EPS,PTYPE,Expansion = False):
        self.Exp = Expansion
        self.L = L
        self.EPS = EPS
        self.PTYPE = PTYPE


    def GetBestAggregate(self,Nu,Gamma,bd,bf,bb,nuRatio=1):
        P = Conv.AnalyticToSimul(nu = Nu,Gamma=Gamma,l=self.L,epsilon=self.EPS,writting=False,ParticleType=self.PTYPE,nu2=Nu*nuRatio)
        w,Ef = bf.Get_Best_Fiber(P,type=1)
        n,Ed = bd.Get_Best_Disk(P)
        Order,Eb = bb.Get_Best_Bulk(P)
        if P.nu == P.nu2:
            EBulk = P.FB
        else:
            EBulk = MP.GetEBulk(Parameter=P)[1]
        if Ef<Eb and Ef<Ed and Ef<EBulk:
            return np.array([0,w,0])
        elif Ed<Eb and Ed < EBulk:
            return np.array([n,0,0])
        elif Eb < EBulk:
            return np.array([0,0,Order+2])
        else:
            return np.array([0,0,1])
    def MakePhaseDiagram(self,Numin,Numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,WidthMax,OrderMax,nuRatio=1):
        Nu,Gamma = np.linspace(Numin,Numax,NpointsNu,dtype=float), np.linspace(Gammamin,Gammamax,NpointsGamma,dtype=float)
        Gamma,Nu = np.meshgrid(Gamma,Nu)
        Color = np.array([np.array([np.zeros(3,dtype=float) for _ in range(Nu.shape[1])]) for _ in range(Nu.shape[0])])
        for i in range(Nu.shape[0]):
            P = Conv.AnalyticToSimul(nu = Nu[i,0],
                                    Gamma=Gamma[i,0],
                                    l=self.L,
                                    epsilon=self.EPS,
                                    writting=False,
                                    ParticleType=self.PTYPE,
                                    nu2=Nu[i,0]*nuRatio)
            bd = BD(Nmax,P,Expansion=self.Exp)
            bf = BF(WidthMax,P,Expansion=self.Exp)
            bb = BB(OrderMax,Nmax,P,Expansion=self.Exp)
            print(Nu[i,0])
            for j in range(Nu.shape[1]):
                #nu[i,j] = cte for i
                Color[i,j] = self.GetBestAggregate(Nu[i,j],Gamma[i,j],bd,bf,bb,nuRatio)
        return Gamma,Nu,Color
