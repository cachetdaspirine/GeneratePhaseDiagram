import numpy as np
import Conversion as Conv
from Numeric_Hex_Energy import *
from Numeric_Fiber_Energy import *
from Numeric_Bulk_Energy import *

class Generate:
    def __init__(self,L,EPS,PTYPE):
        self.L = L
        self.EPS = EPS
        self.PTYPE = PTYPE
    def GetBestAggregate(self,Nu,Gamma,bd,bf,bb):
        P = Conv.AnalyticToSimul(nu = Nu,Gamma=Gamma,l=self.L,epsilon=self.EPS,writting=False,ParticleType=self.PTYPE)
        w,Ef = bf.Get_Best_Fiber(P)
        n,Ed = bd.Get_Best_Disk(P)
        Order,Eb = bb.Get_Best_Bulk(P)
        # print('Nu = '+str(Nu))
        # print('Gamma = '+str(Gamma))
        # print('Best Fiber energy/width ='+str(Ef)+' '+str(w))
        # print('Best Disk energy/N = '+str(Ed)+' '+str(n))
        # print('Best Lacunar Bulk energy/order ='+str(Eb)+' '+str(Order))
        # print('Bulk free energy = '+str(P.FB))
        # print('Lacunar order 0 = '+str(P.Flacune))

        if Ef<Eb and Ef<Ed and Ef<P.FB :
            return np.array([0,w,0])
        elif Ed<Eb and Ed < P.FB:
            return np.array([n,0,0])
        elif Eb < P.FB:
            return np.array([0,0,Order+2])
        else :
            return np.array([0,0,1])
    def MakePhaseDiagram(self,Numin,Numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,WidthMax,OrderMax):
        Nu,Gamma = np.linspace(Numin,Numax,NpointsNu,dtype=float), np.linspace(Gammamin,Gammamax,NpointsGamma,dtype=float)
        Gamma,Nu = np.meshgrid(Gamma,Nu)
        Color = np.array([np.array([np.zeros(3,dtype=float) for _ in range(Nu.shape[1])]) for _ in range(Nu.shape[0])])
        for i in range(Nu.shape[0]):
            P = Conv.AnalyticToSimul(nu = Nu[i,0],
                                    Gamma=Gamma[i,0],
                                    l=self.L,
                                    epsilon=self.EPS,
                                    writting=False,
                                    ParticleType=self.PTYPE)
            bd = BD(Nmax,P)
            bf = BF(WidthMax,P)
            bb = BB(OrderMax,Nmax,P)
            print(Nu[i,0])
            for j in range(Nu.shape[1]):
                #nu[i,j] = cte for i
                Color[i,j] = self.GetBestAggregate(Nu[i,j],Gamma[i,j],bd,bf,bb)
        return Gamma,Nu,Color
