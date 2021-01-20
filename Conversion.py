import numpy as np
from scipy.optimize import newton
class SimulToAnalytic:
    def __init__(self,k,kA,kc,J,epsilon,writting=True,ParticleType='Triangle'):
        if not isinstance(k,float) or (not isinstance(kA,float)) or (not isinstance(kc,float)) or (not isinstance(epsilon,float)):
                print('k='+str(k))
                print('kA='+str(kA))
                print('kc='+str(kc))
                print('j='+str(J))
                print('eps='+str(epsilon))
                raise TypeError('bad type of variable')
        self.k,self.kA,self.kc,self.J,self.epsilon,self.ParticleType = k,kA,kc,J,epsilon,ParticleType #store the constant
        self.nu=(np.sqrt(3)*self.k+2*self.kA)/(3*np.sqrt(3)*self.k+2*self.kA) # Poisson's ratop
        if ParticleType=='Triangle':
            self.la=np.sqrt(3)/2*self.k+self.kA # lambda lame coef
            self.mu=np.sqrt(3)/2*self.k # mu lame coef
            self.kcouplingc=self.kc*4*np.sqrt(3) # coupling parameter
        elif ParticleType=='Hexagon':
            self.la = np.sqrt(3)/4*k+0.5*kA
            self.mu = np.sqrt(3)/4*k
            self.kcouplingc = np.sqrt(3)*kc/(1-epsilon**2)
        self.LMu=self.la+2*self.mu # lambda+2mu
        self.l=np.sqrt((self.la+2*self.mu)/(2*self.kcouplingc)) # depth length
        self.fb=0.5*(self.LMu)*(1+self.nu)*(2*self.epsilon/(1+self.epsilon))**2 # bulk free energy per area
        if self.ParticleType == 'Triangle':
            self.FB=self.fb*np.sqrt(3)/4 #Bulk free energy per particle
            self.Gamma=2*self.J/(self.l*self.fb*(1+self.nu)) # rescaled surface tension
        elif self.ParticleType == 'Hexagon':
            self.FB = self.fb * 3*np.sqrt(3)/2*(1./3.+epsilon**2)
            self.Gamma=2*self.J/(self.l*(1+self.nu)*self.fb*(1/3+epsilon**2)**0.5)
        if writting:
            print('nu='+str(self.nu))
            print('lambda='+str(self.la))
            print('mu='+str(self.mu))
            print('kcoupling_continuous='+str(self.kcouplingc))
            print('l='+str(self.l))
            print('volumique bulk free energy fb='+str(self.fb))
            print('free energy per particle FB='+str(self.FB))
            print('Gamma='+str(self.Gamma))
            print('ParticleType = '+str(self.ParticleType))
        def Range(self,Nmax):
            if self.ParticleType == 'Triangle':
                Size = int(4*(Nmax/6)**0.5+0.5)
                return np.arange(2,Size,2)
            elif self.ParticleType == 'Hexagon':
                Size = int(0.5* (1+np.sqrt(1+8*Nmax)))
                return np.arange(1,Size,1)
class AnalyticToSimul:
    def __init__(self,nu=1/3,Gamma=0,l=1.,epsilon=0.1,writting = True,ParticleType='Triangle'): #here we assumed k = 1
        if (not isinstance(nu,float)) or (not isinstance(Gamma,float)) or (not isinstance(l,float)) or (not isinstance(epsilon,float)):
            print('nu='+str(nu))
            print('Gamma='+str(Gamma))
            print('l='+str(l))
            print('epsilon='+str(epsilon))
            print(isinstance(epsilon,float))
            raise TypeError('bad type of variable')
        self.nu,self.Gamma,self.l,self.epsilon,self.ParticleType = nu,Gamma,l,epsilon,ParticleType
        self.k = 1
        self.kA = (3*np.sqrt(3)*self.nu-np.sqrt(3))/(2*(1-self.nu))
        if ParticleType=='Triangle':
            self.LMu = 3*np.sqrt(3)/2+self.kA
            self.kc = self.LMu/(8*np.sqrt(3)*self.l**2)
        elif ParticleType=='Hexagon':
            self.LMu = 3*np.sqrt(3)/4+0.5*self.kA
            self.kc = (1-epsilon**2)*self.LMu/(2*np.sqrt(3)*self.l**2)
        self.fb = 0.5*self.LMu*(1+self.nu)*(2*self.epsilon/(1+self.epsilon))**2
        if ParticleType=='Triangle':
            self.FB = self.fb*np.sqrt(3)/4 #Bulk free energy per particle
            self.J = self.l*self.fb*(1+self.nu)*self.Gamma/2
        elif ParticleType=='Hexagon':
            self.FB = self.fb* 3*np.sqrt(3)/2*(1./3.+epsilon**2/3.)
            self.J = self.Gamma*(self.l*(1+self.nu)*self.fb*(1/3+epsilon**2/3.))/2
        if writting :
            print('k='+str(self.k))
            print('kA='+str(self.kA))
            print('kc='+str(self.kc))
            print('J='+str(self.J))
            print('fb='+str(self.fb))
            print('ParticleType = '+str(self.ParticleType))
    def Range(self,Nmax):
        if self.ParticleType == 'Triangle':
            #Size = int(4*(Nmax/6)**0.5+0.5)
            Size = int(4*np.sqrt(Nmax/np.pi*0.433))+2
            return np.arange(4,Size,4)
        elif self.ParticleType == 'Hexagon':
            #Size = int(0.5* (1+np.sqrt(1+8*Nmax)))
            #return np.arange(1,Size,1)
            def NFunc(R):
                return (R/2-R**(1./3.))**2-Nmax
            Size = int(newton(NFunc,20))
            return np.arange(4,Size,2)
    def Size(self,N):
        if self.ParticleType == 'Triangle':
            return max(int(4*np.sqrt(N/np.pi*0.433))+2,1)
        elif self.ParticleType== 'Hexagon' :
            def NFunc(R):
                return (R/2-R**(1./3.))**2-N
            Size = int(newton(NFunc,20))
            return max(Size,1)
    def HRange(self,Nmax):
        if self.ParticleType == 'Triangle':
            Size = int(4*(Nmax/6)**0.5+0.5)
            return np.arange(2,Size,2)
        elif self.ParticleType == 'Hexagon':
            Size = int(1./6.* (3+np.sqrt(3*(4*Nmax-1))))
            return np.arange(2,Size,2)
    def write(self,All=False):
            return np.arange(1,Size,1)
    def HSize(self,N):
        if self.ParticleType == 'Triangle':
            return max(int(4*(N/6)**0.5),1)
        elif self.ParticleType== 'Hexagon' :
            #return max(int(0.5* (1+np.sqrt(1+8*N))),1)
            return max(int(1./6.* (3+np.sqrt(3*(4*N-1)))),1)
    def write(self,All=False):
            print('k='+str(self.k))
            print('kA='+str(self.kA))
            print('kc='+str(self.kc))
            print('J='+str(self.J))
            print('fb='+str(self.fb))
            print('ParticleType = '+str(self.ParticleType))
            if All :
                print('epsilon='+str(self.epsilon))
                print('nu='+str(self.nu))
                print('Gamma='+str(self.Gamma))
                print('l='+str(self.l))
