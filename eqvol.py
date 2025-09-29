import numpy as np
import math
import scipy as sc
import scipy.interpolate
from scipy.interpolate import PchipInterpolator
from mpl_toolkits import mplot3d

def getAdjacentIndices(v,t):
    '''
    1. left=-1, right=0 ==> below the first 
    2. left=len(v)-1, right = -1 ==> above the last
    3. left=x, right=x+1 ==> t>=v[left]  
    '''
    if t<=v[0]: return (-1,0)
    if t>=v[-1]: return (len(v)-1, -1)
    for i in range(len(v)):
        if t>=v[i] and t<v[i+1]:
            return (i,i+1)
    raise Exception("Failed to get adjacent index")

class paramvol:
    def __init__(self,dc):
        self.spot = dc.get("spot",3000.0)
        self.strikes = dc.get("strikes",np.linspace(500.0,5000.0,100))
        self.maturities = dc.get("maturities",np.linspace(0.25,10,25))
        self.atmvols = dc.get("atmvols",paramvol.getdefaultvols(self.maturities))
        self.borrow = dc.get("borrow",0.0)
        self.riskfree = dc.get("rfrate",0.0)
        self.paramvoltype = dc.get("paramvoltype","parabolic")
        if self.paramvoltype == "parabolic":
            self.lskew = dc.get("lskew",np.linspace(0.005,0.005,self.maturities.shape[0]))
            self.lconvexity = dc.get("lconvexity",np.linspace(0.01,0.01,self.maturities.shape[0]))
            self.rskew = dc.get("lskew",np.linspace(0.003,0.003,self.maturities.shape[0]))
            self.rconvexity = dc.get("rconvexity",np.linspace(0.003,0.003,self.maturities.shape[0]))
        else:
            self.atminterpolator = PchipInterpolator(self.maturities,self.atmvols*self.atmvols*self.maturities)
            self.eeta = dc.get("eeta",1.6)
            self.lambdav = dc.get("lambdav",0.4)
            self.rho = dc.get("rho",-0.15)
        
    @staticmethod
    def getdefaultvols(maturities,vol0=0.3,voln=0.3):
        fwdvols = np.linspace(vol0,voln,maturities.shape[0])
        vars = np.zeros((fwdvols.shape[0]))
        v = np.zeros((fwdvols.shape[0]))
        varprev = 0.0
        tprev = 0.0
        for i in range(fwdvols.shape[0]):
            vars[i] = ((fwdvols[i]**2)*(maturities[i]-tprev)+varprev)
            varprev = vars[i]
            tprev = maturities[i] 
            v[i]= np.sqrt(vars[i]/maturities[i])
        return v

    def sviVar(self,strike,t):
        y = np.log(strike/(self.spot+(self.riskfree-self.borrow)*t))
        theta = self.atminterpolator(t)
        phi = self.eeta*np.power(theta,-self.lambdav)
        return 0.5*theta*(1.0+self.rho*phi*y+np.sqrt(np.power(phi*y+self.rho,2)+(1.0-self.rho*self.rho)))
    
    def getvol(self,strike,timeindex):
        if self.paramvoltype == "parabolic":
            x = (strike/self.spot-1.0)
            if x<0.:
                if(strike<=self.strikes[0]):
                    x= (self.strikes[0]/self.spot-1.0)
                return self.atmvols[timeindex] + self.lskew[timeindex]*x + self.lconvexity[timeindex]*x*x
            else:
                if(strike>=self.strikes[-1]):
                    x= (self.strikes[-1]/self.spot-1.0)
                return self.atmvols[timeindex] + self.rskew[timeindex]*x + self.rconvexity[timeindex]*x*x
        else:
            atmvol = self.atmvols[timeindex]
            t=self.maturities[timeindex]
            var = self.sviVar(strike,t)
            return np.sqrt(var/t)
    
    def generatesurface(self):
        vols=np.zeros((self.maturities.shape[0],self.strikes.shape[0]))
        for i in range(self.maturities.shape[0]):
            for j in range(self.strikes.shape[0]):
                vols[i,j] = self.getvol(self.strikes[j],i)
        return vols
    
class eqvol:
    def __init__(self,dc):
        self.spot = dc.get("spot",3000.0)
        self.strikes = dc.get("strikes",np.linspace(50.0,25000.0,500))
        self.maturities = dc.get("maturities",np.linspace(0.25,10,25))
        self.isflatvol = dc.get("isflatvol",False)
        self.isparameteric= False
        self.borrow = dc.get("borrow",0.0)
        self.riskfree = dc.get("rfrate",0.0)
        if self.isflatvol:
            self.vols = np.ones((self.maturities.shape[0], self.strikes.shape[0]))*dc.get("vol", 0.3)
        else:
            self.isparameteric = dc.get("isparameteric",True)
            if self.isparameteric:
                self.paramvol = dc.get("paramvol",paramvol({"spot":self.spot,"maturities":self.maturities,"strikes":self.strikes,
                                                            "rfrate":self.riskfree,"borrow":self.borrow,"paramvoltype":"svi"}))
                self.vols = self.paramvol.generatesurface()
            else:
                vol = dc.get("vol", 0.3)
                if isinstance(vol,(float,int)):
                    self.vols = np.ones((self.maturities.shape[0],self.strikes.shape[0]))*vol
                else:
                    self.vols = vol
        self.strikevolinterp= []
        for i in range(self.maturities.shape[0]):
            atmvol = self.getatmvol(i)
            strikestransformed = np.log(self.strikes/self.spot)/(atmvol*math.sqrt(self.maturities[i]))
            self.strikevolinterp.append(sc.interpolate.PchipInterpolator(strikestransformed,self.vols[i,:]))
    
    def getfwd(self,t):
        return self.spot*math.exp((self.riskfree-self.borrow)*t)
    
    def getatmvol(self,ti):
        if self.isparameteric:
            atmvol = self.paramvol.atmvols[ti]
        else:
            if self.isflatvol:
                atmvol = self.vols[0,0]
            else:
                adj = getAdjacentIndices(self.vols[ti,:],self.spot)
                if adj[0]==-1: adj[0]=0
                if adj[1]==-1: adj[1]=self.strikes.shape[0]-1 
                if adj[0]==adj[1]: return self.vols[ti,adj[0]]
                atmvol = self.vols[ti,adj[0]] + (self.vols[ti,adj[1]]-self.vols[ti,adj[0]])*(self.spot-self.strikes[adj[0]])/(self.strikes[adj[1]]-self.strikes[adj[0]])
        return atmvol
    
    def getVol(self,strike,time):
        if self.isflatvol: return self.vols[0][0]
        adj = list(getAdjacentIndices(self.maturities,time))
        if adj[0]==-1: 
            #return self.getVol(strike,self.maturities[0])
            adj[0]=0
        if adj[-1]==-1: 
            #return self.getVol(strike,self.maturities[-1])
            adj[-1]=adj[-2]
        if self.isparameteric:
            if self.paramvol.paramvoltype == "parabolic":
                vol1 = self.paramvol.getvol(strike,adj[0])
                vol2 = self.paramvol.getvol(strike,adj[1])
            else:
                return math.sqrt(self.paramvol.sviVar(strike,time)/time)
        else:
            if strike<self.strikes[0]:
                strike = self.strikes[0]
            elif strike>self.strikes[-1]:
                strike = self.strikes[-1]

            atmvol1 = self.getatmvol(adj[0])
            atmvol2 = self.getatmvol(adj[1])
            scaledstrikes1 = np.log(strike/self.spot)/(atmvol1*math.sqrt(self.maturities[adj[0]]))
            scaledstrikes2 = np.log(strike/self.spot)/(atmvol2*math.sqrt(self.maturities[adj[1]]))
            vol1 = self.strikevolinterp(scaledstrikes1,adj[0])
            vol2 = self.strikevolinterp(scaledstrikes2,adj[1])
        t1=self.maturities[adj[0]]
        t2=self.maturities[adj[1]]
        if(t2==t1): return vol1
        return math.sqrt((vol2*vol2*t2-vol1*vol1*t1)/(t2-t1))

class localvol:
    def __init__(self,eqvol,dc=None):
        if dc is None:
            dc = {}
        self.eq = eqvol
        self.nx = dc.get("nx",500)
        self.nt = dc.get("nt",500)
        self.init()

    def init(self):
        self.times = np.linspace(self.eq.maturities[0],self.eq.maturities[-1],self.nt)
        self.strikes = np.linspace(self.eq.strikes[0]/5,4*self.eq.strikes[-1],self.nx)
        self.forwards = np.asarray([self.eq.getfwd(x) for x in self.times])
        self.lvs = np.zeros((self.strikes.shape[0],self.times.shape[0]))
        #fp=open("c:/temp/lv.csv","w")
        for i in range(self.strikes.shape[0]):
            for j in range(self.times.shape[0]):
                self.lvs[i,j] = self.getLocalVol(i,j)  
                #fp.write("{},".format(self.lvs[i,j]))
            #fp.write("\n")
        #fp.close()
        self.lvinterps = []
        for i in range(self.times.shape[0]):
            self.lvinterps.append(sc.interpolate.PchipInterpolator(np.log(self.strikes/self.eq.spot), self.lvs[:,i]))


    def var(self,k,t):
        v=self.eq.getVol(k,t)
        return v*v*t
    
    def getLocalVol(self,ki,tj):
        delt = 0.0001
        delk = 0.0001
        k=self.strikes[ki]
        t=self.times[tj]
        fwd = self.forwards[tj]
        y = math.log(k/fwd) #dy =fwd/k dk
        w = self.var(k,t)
        if t<=self.times[0] or t>=self.times[-1]:
            v=self.eq.getVol(k,t)
            dwdt = v*v
        else:
            dwdt = (self.var(k,t+delt)-w)/delt
        if k<=self.strikes[0] or k >= self.strikes[-1]:
            dwdy=0.0
            dw2dy2=0.0
        else:
            dwdy = ((self.var(k+delk,t)-w)/delk)*(k/fwd)
            dw2dy2 = ((self.var(k+delk,t)-2.0*w+self.var(k-delk,t))/(delk*delk))*(-k*k/fwd)

        den = 1.0 - (y/w)*dwdy  +0.25*(-0.25-1.0/w + y*y/(w*w))*dwdy*dwdy + 0.5*dw2dy2
        return math.sqrt(max(dwdt/den,1e-8))

    def getLV(self,s,t):
        adj = list(getAdjacentIndices(self.times,t))
        if adj[0]==-1: 
            adj[0]=0
        if adj[-1]==-1: 
            adj[-1]=adj[-2]
        
        vol1 = self.lvinterps[adj[0]](s)
        vol2 = self.lvinterps[adj[1]](s)
        t1=self.times[adj[0]]
        t2=self.times[adj[1]]
        if(t2==t1): return vol1
        return vol1+(t-t1)*(vol2-vol1)/(t2-t1)
        #return np.sqrt(np.maximum((vol2*vol2*t2-vol1*vol1*t1)/(t2-t1),np.ones(s.shape)*1e-8))

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    vol = eqvol({"isflatvol":False})
    
    m,s=np.meshgrid(vol.maturities,vol.strikes)
    vols = vol.vols
    my_cmap = plt.get_cmap('hot')
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(m,s,vols.T,cmap = my_cmap,edgecolor ='none')
    #ax.set(xticklabels=[], yticklabels=[], zticklabels=[])  
    ax.set_xlabel("maturities")
    ax.set_ylabel("strikes")
    plt.show()

    lv = localvol(vol)
    m,s=np.meshgrid(lv.times,lv.strikes)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(m,s,lv.lvs.T,cmap = my_cmap,edgecolor ='none')
    #ax.set(xticklabels=[], yticklabels=[], zticklabels=[])  
    ax.set_xlabel("maturities")
    ax.set_ylabel("spots")
    plt.show()



    time = 8.2
    strikes=np.linspace(100,20000,500)
    vs=[]
    for strike in strikes:
        vs.append(vol.getVol(strike,time))
    plt.plot(strikes,vs)
    plt.show()
    

