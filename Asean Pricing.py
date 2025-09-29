import numpy as np
import math
def SimulateSt(S0,r,sigma,t,nsim=10000):
    randnumbers = np.random.normal(0.0,1.0,(nsim,))*sigma*math.sqrt(t)
    return S0*np.exp((r-0.5*sigma*sigma)*t + randnumbers)

def optionPrice(optype,S0=100.0,k=100.0,r=0.05,mu=0.05,sigma=0.2,t=1.2,nsim=10000):
    St = SimulateSt(S0,mu,sigma,t,nsim)
    df = math.exp(-r*t)
    if optype==1: #call
        return df*np.mean(np.maximum(St-k,0.0))
    else:
        return df*np.mean(np.maximum(k-St,0.0))
    

def Asean(optype=1,S0=100.0,k=100.0,r=0.05,mu=0.05,sigma=0.2,t=1.2,n=12,nsimWhole=100):
    timeList = np.linspace(0,t,n+1)
    Stockmatrix = np.zeros((n+1,nsimWhole))
    Stockmatrix[0] = S0
    Options = np.zeros((nsimWhole,))
    for i in range(Stockmatrix.shape[0]-1):
        Stockmatrix[i+1][:] = SimulateSt(Stockmatrix[i][:],r,sigma,t/n,nsimWhole)
    Stockmatrix = np.transpose(Stockmatrix)
    for i in range(nsimWhole):
        avg = sum(Stockmatrix[i][1:])/n
        if optype == 1:
            Options[i] = max(avg - k,0.0)
        else:
            Options[i] = max(k-avg,0.0)
    return np.mean(Options)*math.exp(-r*t)
    '''CumulativeOpPrice = 0
    for i in range(nsimWhole):
        GrandSum = 0
        for j in timeList:
            StockP = SimulateSt(S0,r,sigma,j,nsim)
            GrandSum = GrandSum + StockP
        AvgPrice = GrandSum/(n+1)
        if optype == 1:#Call
            CumulativeOpPrice = CumulativeOpPrice + max(AvgPrice - k,0)
        if optype == 0:#Put
            CumulativeOpPrice = CumulativeOpPrice + max(k-AvgPrice,0)
    return CumulativeOpPrice/nsimWhole
    '''
a = Asean(optype=0,S0=100.0,k=100.0,r=0.05,mu=0.05,sigma=0.2,t=1.2,n=12,nsimWhole=100)
print(a)

b = optionPrice(optype=0,S0=100.0,k=100.0,r=0.05,mu=0.05,sigma=0.2,t=1.2,nsim=10000)
print(b)



