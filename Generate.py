#!/usr/bin/python3
import numpy as np
import Conversion as Conv
from Numeric_Hex_Energy import *
from Numeric_Fiber_Energy import *

class Generate:
    def __init__(self,L,EPS,PTYPE):
        self.L = L
        self.EPS = EPS
        self.PTYPE = PTYPE
    def GetBestAggregate(self,Nu,Gamma,bd,bf):
        P = Conv.AnalyticToSimul(nu = Nu,Gamma=Gamma,l=self.L,epsilon=self.EPS,writting=False,ParticleType=self.PTYPE)
        w,Ef= bf.Get_Best_Fiber(P)
        n,Ed = bd.Get_Best_Disk(P)
        if Ef<P.FB and Ef<Ed :
            return w
        elif Ed<P.FB:
            return -n
        else :
            return 0
    def MakePhaseDiagram(self,Numin,Numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,WidthMax):
        Nu,Gamma = np.linspace(Numin,Numax,NpointsNu,dtype=float), np.linspace(Gammamin,Gammamax,NpointsGamma,dtype=float)
        Gamma,Nu = np.meshgrid(Gamma,Nu)
        Color = np.array([np.zeros(Nu.shape[1],dtype=float) for _ in range(Nu.shape[0])])
        for i in range(Nu.shape[0]):
            P = Conv.AnalyticToSimul(nu = Nu[i,0],
                                    Gamma=Gamma[i,0],
                                    l=self.L,
                                    epsilon=self.EPS,
                                    writting=False,
                                    ParticleType=self.PTYPE)
            bd = BD(Nmax,P)
            bf = BF(WidthMax,P)
            print(Nu[i,0])
            for j in range(Nu.shape[1]):
                #nu[i,j] = cte for i
                Color[i,j] = self.GetBestAggregate(Nu[i,j],Gamma[i,j],bd,bf)
        return Gamma,Nu,Color
