import numpy as np
import copy

TopologieDownHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieUpHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieDownTriangle = [(1,0),(-1,0),(0,1)]
TopologieUpTriangle = [(1,0),(-1,0),(0,-1)]

def Parallel0(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<=size and j<i and j+i>=size/2 and i+j<=3*size/2-1:
                Res[i,j]=1
    return Res

def Parallel3(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<size and j<i and j+i>size/2 and i+j<=3*size/2-3:
                Res[i,j]=1
    return Res

def Parallel1(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<size and j<i and j+i>size/2-1 and i+j<=3*size/2-2:
                Res[i,j]=1
    return Res

def Parallel2(size):
    Res=np.array([np.zeros(size//2,dtype=int) for _ in range(2*size)])
    for i in range(2*size):
        for j in range(size//2):
            if i-j<=size and j<i and j+i>size/2 and i+j<3*size/2-1:
                Res[i,j]=1
    return Res
def ParallelHex(size):
    Res=np.array([np.zeros(3*size,dtype=int) for _ in range(2*size)])
    Res[Res.shape[0]//2,Res.shape[1]//2] = 1
    for k in range(1,size):
        Resdt = copy.copy(Res)
        for i in range(Res.shape[0]):
            for j in range(Res.shape[1]):
                if Res[i,j]==1:
                    ij = (i,j)
                    for neigh in Get_Neighbors(Resdt,ij,Free = True,ParticleType = 'Hexagon'):
                        Resdt[neigh] = 1
        Res = Resdt

    return Res

#---------------------------------------------------------------
# I am not sure that it is necessary, but depending on the value
# of the size, there are 4 different way of  making  an  hexagon
# At the end the function Parallel ensure to return the  correct
# Hexagon
#---------------------------------------------------------------
#N = 6 * (Size/4)**2
def Parallel(size,ParticleType='Triangle'):
    if ParticleType=='Triangle':
        if size <4 :
            return np.array([[0,0,0],[0,1,0],[0,0,0]])
        if size%4==0:
            return Parallel0(size)
        elif size%4==1:
            return Parallel1(size)
        elif size%4==2:
            return Parallel2(size)
        elif size%4==3:
            return Parallel3(size)
    elif ParticleType=='Hexagon':
        return ParallelHex(size)
#---------------------------------------------------------------
# Given an array  of  0 and 1 : returns the number of 1, i.e the
# Number of particle
#---------------------------------------------------------------
def Np(Array):
    unique, counts = np.unique(Array, return_counts=True)
    return dict(zip(unique, counts))[1]

#---------------------------------------------------------------
# Make a fiber of certain width, and length.
#---------------------------------------------------------------
def Fiber(Width,Length,ParticleType='Triangle'):
    if ParticleType=='Triangle':
        Res=np.array([np.zeros(Width+2,dtype=int) for _ in range(Length+2)])
        for i in range(1,Length+1):
            for j in range(1,Width+1):
                Res[i,j]=1
        return Res
    elif ParticleType=='Hexagon':
        return Fiber2(Width,Length)
def Fiber2(W,L):
        Res=np.array([np.zeros(L) for _ in range(W+int(L//2+0.5))])
        for l in range(L):
            if l%2==0:
                Res[l//2:W+l//2,l]=1
            else :
                Res[l//2+1:W+l//2,l]=1
        return np.flip(np.transpose(Res),0)
def Fiber3(W,L):
        Res=np.array([np.zeros(L) for _ in range(W+int(L//2+0.5))])
        for l in range(L):
            if l%2==0:
                Res[l//2:W+l//2,l]=1
            else :
                Res[l//2:W+l//2,l]=1
        return np.flip(np.transpose(Res),0)
def Fiber4(W,L):
    Res=np.array([np.zeros(W) for _ in range(L)])
    Res[:]=1
    return np.transpose(Res)

def Get_Neighbors(Array,ij,Occupied=False,Free=False,Border=False,ParticleType = 'Triangle'):
        Lx, Ly = Array.shape[0], Array.shape[1]
        if ParticleType == 'Triangle':
            TopologieUp = TopologieUpTriangle
            TopologieDown = TopologieDownTriangle
        elif ParticleType=='Hexagon':
            TopologieUp = TopologieUpHex
            TopologieDown = TopologieDownHex
        if (ij[0]+ij[1])%2==0:
            Res = np.array(TopologieDown)+np.array(ij)
        else :
            Res = np.array(TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        if not Border:
            Resreg=np.delete(Res,np.argwhere((Res[:,0]>=Lx) | (Res[:,0]<0) | (Res[:,1]>=Ly) | (Res[:,1]<0)),0)
        else:
            Resreg=Res
        #Build a numpy array of tuple
        Resbis=np.empty(Resreg.__len__(),dtype=object)
        Resbis[:] = list(zip(Resreg[:,0],Resreg[:,1]))
        if Border:
            Resbis = list(Resbis)
        #check the occupancie or not
        if Border :
            for n in reversed(range(Resbis.__len__())):
                if Resbis[n][0]>=0 and Resbis[n][1]>=0 and Resbis[n][0]<Lx and Resbis[n][1]<Ly:
                    if Array[Resbis[n]]!=0:
                        del Resbis[n]
        else :
            if Occupied:
                Resbis=Resbis[np.array([Array[r]==1 for r in Resbis])]
            elif Free:
                Resbis = Resbis[np.array([Array[r]==0 for r in Resbis])]

        return set(Resbis)
def SurfaceEnergy(Array,J=1.,ParticleType='Triangle'):
    surface = 0
    for i in range(Array.shape[0]):
        for j in range(Array.shape[1]):
            ij = (i,j)
            if Array[ij] == 1:
                surface+=Get_Neighbors(Array,ij,Free=True,Border=True,ParticleType=ParticleType).__len__()
    return J*surface
