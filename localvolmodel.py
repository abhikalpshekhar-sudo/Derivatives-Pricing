import numpy as np
import eqvol 
import math
class lvmodelsimulator:
    def __init__(self,eq,dc=None):
        if dc is None: dc={}
        self.lv = eqvol.localvol(eq,dc)
        self.eq=eq
        self.r = eq.riskfree-eq.borrow
        self.maxsimdates = dc.get("maxsimdate",3.0)
        self.mandatorydates = dc.get("mandatorydates",[])
        self.maxgap = dc.get("maxgap",7.0/365.0)
        self.simdates,self.mdtindices = lvmodelsimulator.simulationdates(self.maxsimdates,self.mandatorydates,self.maxgap)
        self.nsim = dc.get("nsim",100)
        self.seed = dc.get("seed",97)
        self.forwardcorrect = dc.get("forwardcorrect",True)
        self.paths = np.zeros((self.nsim,1+len(self.simdates))) 
        #self.paths[:,0] = np.log(self.eq.spot)
        np.random.seed(self.seed)
        z=np.random.normal(0.0,1.0,(int(self.nsim/2),len(self.simdates)))
        self.rnd = np.vstack([z,-z])

    @staticmethod
    def simulationdates(T,dts,delt):
        mingap = 1e-6
        delt = max(delt,mingap)
        dates = dts
        if len(dts)==0 or dts[-1]<T:
            dates.append(T)
        finaldates=[delt]
        mdates = []
        for i in range(len(dates)):
            while dates[i] - finaldates[-1] > 2.0*delt:
                finaldates.append(delt+finaldates[-1])
            finaldates.append(dates[i])
            mdates.append(len(finaldates)-1)
        return finaldates,mdates
    
    def evolvenext(self,ti):
        if ti>0:
            t=self.simdates[ti-1]
        else:
            t=0.0
        tnext=self.simdates[ti]
        dt = tnext-t
        tmid = 0.5*(t+tnext)
        if abs(t-2.1671232876712336)<1e-6:
            print("check me")
        lvs = self.lv.getLV(self.paths[:,ti],tmid)
        lvs = np.minimum(np.maximum(lvs,np.full(lvs.shape,1e-4)),np.full(lvs.shape,10000.0))
        self.paths[:,ti+1] = self.paths[:,ti] + self.r*dt-0.5*lvs*lvs*dt + lvs*self.rnd[:,ti]*math.sqrt(dt)
        print("time={}\tlv_max={}\tlv_min={}\tspot_max={}\tspot_min={}\tspot_mean={}".format(t,
                                                                               np.max(lvs),
                                                                               np.min(lvs),
                                                                               self.eq.spot*(np.exp(np.max(self.paths[:,ti+1]))),
                                                                               self.eq.spot*(np.exp(np.min(self.paths[:,ti+1]))),
                                                                               self.eq.spot*(np.mean(np.exp(self.paths[:,ti+1])))))
        if self.forwardcorrect:
            corr=self.eq.spot*np.mean(np.exp(self.paths[:,ti+1]))/self.eq.getfwd(tnext)
            self.paths[:,ti+1] = self.paths[:,ti+1] - np.log(corr)

    def evolve(self):
        for i in range(0,len(self.simdates)):
            self.evolvenext(i)
    
import bsutils
class SimpleBarrierOption:
    def __init__(self,dc=None):
        if dc is None:
            dc = {}
        self.strike = dc.get("strike",3000.0)
        self.barrier = dc.get("barrier",3600.0)
        self.option = dc.get("option",bsutils.OPTIONTYPE.CALL)
        self.rebate = dc.get("rebate",0.0)
        self.maturity = dc.get("maturity",2.0)

    def price(self,eq):
        ls=lvmodelsimulator(eq,{"maxsimdate":self.maturity,"nsim":40000})
        ls.evolve()
        paths=eq.spot*np.exp(ls.paths)
        s = 0.0
        vars = 0.0
        for i in range(ls.paths.shape[0]):
            breached = False
            if self.barrier is not None:
                for j in range(ls.paths.shape[1]):
                    if paths[i,j]>self.barrier:
                        s = s+self.rebate
                        breached = True
                        break
            if not breached:
                p = max(paths[i,-1]-self.strike, 0.)
                s = s + p
                vars = vars + p*p
        price = s/ls.paths.shape[0] 
        serr = math.sqrt((vars/ls.paths.shape[0] - price*price)/ls.paths.shape[0])
        return price, serr

if __name__=="__main__":
    
    plotpath=False
    if plotpath:
        #randomly plot ten paths
        eq=eqvol.eqvol({"isflatvol":False,"rfrate":0.05})
        ls=lvmodelsimulator(eq,{"maxsimdate":9.8,"nsim":500})
        ls.evolve()
        import matplotlib.pyplot as plt
       
        for i in range(50):
            plt.plot([0.0]+ls.simdates,ls.paths[i,:])
        plt.show()

    pricevanilla =False
    if pricevanilla:
        eq=eqvol.eqvol({"isflatvol":True,"rfrate":0.05})
        strike = 3500.0
        op=SimpleBarrierOption({"strike":strike,"barrier":None,"maturity":2})
        def pricefunc(strike):
            vol = eq.getVol(strike,op.maturity)
            fwd = eq.spot*math.exp((eq.riskfree-eq.borrow)*op.maturity)
            return bsutils.bspriceundiscounted(fwd,strike,vol,op.maturity,bsutils.OPTIONTYPE.CALL)
        
        print("Price of vanilla = {}".format(op.price(eq)))
        print("Closed form price is {}".format(pricefunc(strike)))
    priceop=True
    if priceop:
        op=SimpleBarrierOption()
        print("Price of barrier = {}".format(op.price(eqvol.eqvol({"isflatvol":False}))))

