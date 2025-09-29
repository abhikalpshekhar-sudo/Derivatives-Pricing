import bsutils
import numpy as np
import math

class SingleVariateDistr:
    def __init__(self,dc=None):
        if dc is None: dc = {}
        self.nx = dc.get('nx',1000)
        self.callspread = dc.get('callspread', 0.0001)
        self.nsigma = dc.get('nsigma',5)
        self.x=None
        self.y=None

    @staticmethod
    def getExpectationNumInt(distr,f):
        s = 0.0
        den  = 0.0
        for i in range(distr.x.shape[0]):
            s += f(distr.x[i])*distr.pdf[i]
            den += distr.pdf[i]
        #Left and right tail
        s /= den 
        return s
    
    def getdistr(self,eq,rfrate,time,borrowrate=0.0,xmin=None,xmax=None):
        spot = eq.spot
        fwd = spot*math.exp((rfrate-borrowrate)*time)
        atmvol = eq.getVol(spot,time)
        def pricefunc(strike):
            vol = eq.getVol(strike,time)
            return bsutils.bspriceundiscounted(fwd,strike,vol,time,bsutils.OPTIONTYPE.PUT)
              
        callspread = self.callspread
        if xmin is None:
            xmin = fwd*math.exp(-0.5*atmvol*atmvol*time-self.nsigma*atmvol*math.sqrt(time))
        if xmax is None:
            xmax = fwd*math.exp(-0.5*atmvol*atmvol*time+self.nsigma*atmvol*math.sqrt(time))
        x = np.linspace(xmin,xmax,self.nx)
        y = np.zeros(x.shape)
        callspread = x[1]-x[0]
        for i in range(x.shape[0]):
            y[i] = ((pricefunc(x[i]+callspread) - pricefunc(x[i])))/(callspread)
            
        
        self.pdf=(y[1:] - y[:-1]) #pdf * delx
        self.x=x[1:]
        self.y=y[1:]
        

def priceEuropean(func,time,eq,distr=None):
    import eqvol
    if distr is None:
        distr = SingleVariateDistr()
        distr.getdistr(eq,0.03,time,borrowrate=0.0,xmin=None,xmax=None)
    return SingleVariateDistr.getExpectationNumInt(distr,func)

    
if __name__ == "__main__":
    import eqvol
    import matplotlib.pyplot as plt
    
    distr = SingleVariateDistr()
    eq = eqvol.eqvol({"isflatvol":True})
    spot = eq.spot
    time = 3.2
    rfrate = 0.0
    borrowrate = 0.0
    strike = 3000.0
    distr.getdistr(eq,rfrate,time,borrowrate,xmin=None,xmax=None)
    plt.plot(distr.x,distr.y)
    plt.show()
    print("Price of [x-K]+ for atm is = {}".format(priceEuropean(
        lambda x:np.maximum(x-strike,0.0), time,eq,distr)))
    def pricefunc(strike):
        vol = eq.getVol(strike,time)
        fwd = spot*math.exp((rfrate-borrowrate)*time)
        return bsutils.bspriceundiscounted(fwd,strike,vol,time,bsutils.OPTIONTYPE.CALL)
    print("Closed form price is {}".format(pricefunc(strike)))