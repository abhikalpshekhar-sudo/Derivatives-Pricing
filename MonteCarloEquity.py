import numpy as np
import math

def correlatedNormals(Cor_Matrix,nsim):
    a = np.transpose(Cor_Matrix)
    if np.array_equal(a,Cor_Matrix) == False:
        raise Exception("Correlations must be symmetric")
    ls = np.random.normal(0.0,1,(nsim,Cor_Matrix.shape[0]))
    tran_chol = np.transpose(np.linalg.cholesky(Cor_Matrix))
    return ls@tran_chol

'''def optPrice(StockP,strikeP,dfs = IRCurve([1.0],[0.05]),optype = 1):
    if optype == 1:
        return max(StockP-strikeP,0.0)*dfs
    return max(strikeP - StockP,0.0)*dfs'''



class IRCurve:
    def __init__(self,times=[],rates=[]):
        self.times=times
        self.discFactors = [math.exp(-r*t) for r,t in zip(rates,times)]

    def discountFactor(self,t):
        if t<=self.times[0]:
            r = -math.log(self.discFactors[0])/self.times[0]
            return math.exp(-r*t)
        elif t>=self.times[-1]:
            if len(self.discFactors)==1:
                r = -math.log(self.discFactors[0])/self.times[0] 
            else:
                r = -math.log(self.discFactors[-1]/self.discFactors[-2])/(self.times[-1] - self.times[-2])
            return self.discFactors[-1]*math.exp(-r*(t-self.times[-1]))
        
        for i,tt in enumerate(self.times):
            if abs(t-tt)<1e-8:
                return self.discFactors[i]
            elif t<tt:
                tprev=self.times[i-1]
                r = -math.log(self.discFactors[i]/self.discFactors[i-1])/(self.times[i] - self.times[i-1])
                return math.exp(-r*(t-tprev))*self.discFactors[i-1]

def optPrice(StockP,strikeP,t,dfs = IRCurve([1.0],[0.05]),optype = 1):
    disc = dfs.discountFactor(t)
    if optype == 1:
        return max(StockP-strikeP,0.0)*disc
    return max(strikeP - StockP,0.0)*disc
       
class Equity:
    def __init__(self,name="",spot=100.0,vol=0.2,r=None):
        if r is None:
            r = IRCurve([1.0],[0.05])
        self.rcurve = r
        self.name = name
        self.spot = spot
        self.vol = vol

    def evolveNext(self,s,t,tnext,randnumbers):
        growth = self.rcurve.discountFactor(t)/self.rcurve.discountFactor(tnext)
        sigma2t = math.exp(-0.5*(tnext-t)*self.vol**2)
        sigmaSqrtT = self.vol*math.sqrt(tnext-t)
        return s*growth*sigma2t*np.exp(randnumbers*sigmaSqrtT)

class EquityBasket:
    def __init__(self,eqs,corr=None):
        self.eqs=eqs
        self.neq = len(self.eqs)
        if corr is None:
            corr = np.eye(self.neq)
        self.corr=corr
        self.slices = []

    def evolve(self,simtimes,nsim):
        s0=np.zeros((nsim,self.neq))
        for i in range(self.neq):
            s0[:,i] = self.eqs[i].spot
        self.slices.append(s0)
        tprev=0.0
        for i,t in enumerate(simtimes):
            rnds = correlatedNormals(self.corr,nsim)
            st = np.zeros(rnds.shape)
            for j in range(self.neq):
                st[:,j]=self.eqs[j].evolveNext(self.slices[-1][:,j],tprev,t,rnds[:,j])
            self.slices.append(st)
            tprev=t

def defaultMarket():
    EQ1Data = {
        "spot"  : 100.0,
        "vol"   : 0.2,
        "growth": IRCurve([1.0],[0.05])
    }
    EQ2Data = {
        "spot"  : 200.0,
        "vol"   : 0.15,
        "growth": IRCurve([1.0],[0.05])
    }
    marketData = {
        "EQ1" : EQ1Data,
        "EQ2" : EQ2Data,
        "corr": np.asarray([[1,0.6],[0.6,1]])
    }
    return marketData

class BasketOption:
    def __init__(self,eqbasket=["EQ1","EQ2"],weights=[0.6,0.4],strike=170,maturity=5, discount=IRCurve([1.0],[0.05])):
        self.eqbasket=eqbasket
        self.weights=np.asarray(weights)
        self.weights.reshape((2,1))
        self.strike=strike
        self.maturity=maturity
        self.discount = discount
        self.basketPrice = 0.0
    
    def price(self,marketData=None,nsim=10000):
        if marketData is None:
            marketData = defaultMarket()
        eqs = []
        for eq in self.eqbasket:
           eqs.append(Equity(eq,marketData[eq]["spot"],marketData[eq]["vol"],marketData[eq]["growth"]))
        eqbasket = EquityBasket(eqs,marketData["corr"]) 
        eqbasket.evolve([self.maturity],nsim)
        self.basketPrice = eqbasket.slices[-1]@self.weights
        optionprice = np.maximum(self.basketPrice-self.strike,0.0)
        df =  self.discount.discountFactor(self.maturity)
        return np.mean(optionprice)*df


class BasketAsianOption:
    def __init__(self,eqbasket=["EQ1","EQ2"],weights=[0.6,0.4],
                 averagingDates=[1.0,2.0],optype=1,strike=140,discount=IRCurve([1.0],[0.05])):
        self.dates = averagingDates
        self.maturity = averagingDates[-1]
        self.strike = strike
        self.eqbasket = eqbasket
        self.weights=np.asarray(weights)
        self.weights.reshape((2,1))
        self.optype = optype
        self.discount=discount

    def price(self,marketData=None,nsim=10000):
        sum = 0
        if marketData is None:
            marketData = defaultMarket()
        eqs = []
        for eq in self.eqbasket:
           eqs.append(Equity(eq,marketData[eq]["spot"],marketData[eq]["vol"],marketData[eq]["growth"]))
        eqbasket = EquityBasket(eqs,marketData["corr"]) 
        eqbasket.evolve(self.dates,nsim)
        for i,t in enumerate(self.dates):
            basketPrice = eqbasket.slices[i]@self.weights
            if i==0:
                sumbasketprices=basketPrice
            else:
                sumbasketprices = sumbasketprices + basketPrice
        sumbasketprices = sumbasketprices/float(len(self.dates))
        optionprice = np.maximum(sumbasketprices-self.strike,0.0)
        df =  self.discount.discountFactor(self.dates[-1])
        return np.mean(optionprice)*df
        

                
          

if __name__=="__main__":

    def testCurve():
        curve = IRCurve([1.0,2.0,3.0,4.0,5.0],[0.05,0.055,0.056,0.058,0.06])
        [print(f"DF[{t}]={curve.discountFactor(t)}") for t in curve.times] 
        times = [0.5,1.5,2.5,3.5,4.5,5.5]
        [print(f"DF[{t}]={curve.discountFactor(t)}") for t in times] 
    def testEq():
        eq = Equity()
        nsim = 10
        randnums=np.random.normal(0.0,1.0,nsim)
        s=np.asarray([100.0]*nsim)
        st1=eq.evolveNext(s,0.0,0.5,randnums)
        print(st1)
        randnums=np.random.normal(0.0,1.0,nsim)
        st2=eq.evolveNext(st1,0.5,1.0,randnums)
        print(st2)
    def testEQBasketSingleAsset():
        eqbasket = EquityBasket([Equity()])
        simtimes = np.linspace(0.1,1.0,10)
        nsims=5
        eqbasket.evolve(simtimes,nsims)
        for s in eqbasket.slices:
            print(f"{s}")

    def testEQBasketSingleAsset():
        eqbasket = EquityBasket([Equity(),Equity(spot=100.0,vol=0.2)],corr=np.asarray([[1.0,0.99999],[0.99999,1.0]]))
        simtimes = np.linspace(0.1,1.0,10)
        nsims=5
        eqbasket.evolve(simtimes,nsims)
        for s in eqbasket.slices:
            print(f"{s}")

    def testEQBasketOption():
        basketOp=BasketOption()
        print(f"Price = {basketOp.price()}")

    def testAsianBasket():
        asian = BasketAsianOption()
        print(f"Asian Option Price is {asian.price()}")

    testEQBasketOption()   
    #testAsianBasket()