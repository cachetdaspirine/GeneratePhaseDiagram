import numpy as np
import sys
sys.path.append('/home/hleroy/Simulation/Extra_Module_py/')
#sys.path.append('/home/hugo/Extra_Module_py/')
import Conversion as Conv
from Numeric_Hex_Energy import *
from Numeric_Fiber_Energy import *
from Numeric_Bulk_Energy import *
import multiprocessing as mp

class Generate:
    def __init__(self,L,EPS,PTYPE,Expansion = False):
        self.Exp = Expansion
        self.L = L
        self.EPS = EPS
        self.PTYPE = PTYPE
    def GetBestAggregate(self,Nu,Gamma,bd,bf,bb):
        P = Conv.AnalyticToSimul(nu = Nu,Gamma=Gamma,l=self.L,epsilon=self.EPS,writting=False,ParticleType=self.PTYPE)
        if P.ParticleType == 'Triangle':
            w1,Ef1 = bf.Get_Best_Fiber(P)
            w2,Ef2 = 0,np.inf
            n,Ed = bd.Get_Best_Disk(P)
            Order,Eb = bb.Get_Best_Bulk(P)
        else :
            w1,Ef1 = bf.Get_Best_Fiber(P,1)
            w2,Ef2 = bf.Get_Best_Fiber(P,2)
            n,Ed = bd.Get_Best_Disk(P)
            Order,Eb = bb.Get_Best_Bulk(P)
        # print('Nu = '+str(Nu))
        # print('Gamma = '+str(Gamma))
        # print('Best Fiber energy/width ='+str(Ef)+' '+str(w))
        # print('Best Disk energy/N = '+str(Ed)+' '+str(n))
        # print('Best Lacunar Bulk energy/order ='+str(Eb)+' '+str(Order))
        # print('Bulk free energy = '+str(P.FB))
        # print('Lacunar order 0 = '+str(P.Flacune))
        energies = [Ef1,Ef2,Ed,Eb,P.FB]
        NumBest = np.argmin(energies)
        if NumBest == 0 :
            return np.array([0,w1,0])
        elif NumBest == 1:
            return np.array([0,-w2,0])
        elif NumBest == 2:
            return np.array([n,0,0])
        elif NumBest == 3:
            return np.array([0,0,Order+2])
        elif NumBest == 4:
            return np.array([0,0,1])
    def GetColor(self,Nmax,WidthMax,OrderMax,NU,GAMMA):
        P = Conv.AnalyticToSimul(nu = NU[0],
                                Gamma=GAMMA[0],
                                l=self.L,
                                epsilon=self.EPS,
                                writting=False,
                                ParticleType=self.PTYPE)
        bd = BD(Nmax,P,Expansion = self.Exp)
        bf = BF(WidthMax,P,Expansion = self.Exp)
        bb = BB(OrderMax,Nmax,P,Expansion = self.Exp)
        print(NU[0])
        res = np.zeros((NU.shape[0],3))
        for j in range(NU.shape[0]):
            #nu[i,j] = cte for i
            res[j] = self.GetBestAggregate(NU[j],GAMMA[j],bd,bf,bb)
        return res
    def MakePhaseDiagram(self,Numin,Numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,WidthMax,OrderMax):
        Nu,Gamma = np.linspace(Numin,Numax,NpointsNu,dtype=float), np.linspace(Gammamin,Gammamax,NpointsGamma,dtype=float)
        Gamma,Nu = np.meshgrid(Gamma,Nu)
        Color = np.array([np.array([np.zeros(3,dtype=float) for _ in range(Nu.shape[1])]) for _ in range(Nu.shape[0])])
        #pool = mp.Pool(mp.cpu_count())
        #Color = pool.map(self.GetColoGamma,nu,Color = pool.apply(G.MakePhaseDiagram,args=(numin,numax,NpointsNu,Gammamin,Gammamax,NpointsGamma,Nmax,Wmax,OrderMax))r,(Nu,Gamma))
        #for i in range(Nu.shape[0]):
        #Color = np.array([pool.apply(self.GetColor,args=(Nmax,WidthMax,OrderMax,Nu[i,:],Gamma[i,:])) for i in range(Nu.shape[0])])
        Color = np.array([self.GetColor(Nmax,WidthMax,OrderMax,Nu[i,:],Gamma[i,:]) for i in range(Nu.shape[0])])
        #pool.close()
            # P = Conv.AnalyticToSimul(nu = Nu[i,0],
                                    # Gamma=Gamma[i,0],
                                    # l=self.L,
                                    # epsilon=self.EPS,
                                    # writting=False,
                                    # ParticleType=self.PTYPE)
            # bd = BD(Nmax,P,Expansion = self.Exp)
            # bf = BF(WidthMax,P,Expansion = False)
            # bb = BB(OrderMax,Nmax,P,Expansion = self.Exp)
            # print(Nu[i,0])
            # for j in range(Nu.shape[1]):
                #nu[i,j] = cte for i

        return Gamma,Nu,Color
