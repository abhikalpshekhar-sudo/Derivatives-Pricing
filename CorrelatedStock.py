import numpy as np
import math


def optionPrice(optype,S0=100.0,k=100.0,r=0.05,mu=0.05,sigma=0.2,t=1.2,nsim=10000):
    St = SimulateSt(S0,mu,sigma,t,nsim)
    df = math.exp(-r*t)
    if optype==1: #call
        return df*np.mean(np.maximum(St-k,0.0))
    else:
        return df*np.mean(np.maximum(k-St,0.0))

def Correlated_Normals(Cor_Matrix,nsim,t):
    a = np.transpose(Cor_Matrix)
    if np.array_equal(a,Cor_Matrix) == False:
        raise Exception("Correlations must be symmetric")
    ls = np.random.normal(0.0,t,(nsim,Cor_Matrix.shape[0]))
    tran_chol = np.transpose(np.linalg.cholesky(Cor_Matrix))
    return ls@tran_chol
'''def simul_correlatedStocks(Cor_Matrix,Bases,r=0.05,sigmas,t=1.2,nsim=1000):
    RandMat = Correlated_Normals(Cor_Matrix,nsim)*sigmas*math.sqrt(t)
    n = Cor_Matrix.shape[0]
    sum = np.zeros(n)
    for i in range(n):
        temp_sum = 0
        for j in range(nsim):
            tempsum = sum(Bases[i]*math.exp(((r-(sigmas[i]*sigmas[i])/2))*t + RandMat[i][j]))
        sum[i] = temp_sum
    sum = sum*1/nsim
    return sum
'''
def SimulateSt(Cor_Matrix,Bases,sigmas,r=0.05,t=1.2,nsim=10000):
    corr = Correlated_Normals(Cor_Matrix,nsim)
    weighted = np.zeros(corr.shape)
    for i in range(len(corr)):
        weighted[:][i] = corr[:][i]*sigmas[i]*math.sqrt(t)
        

    
    #randnumbers = np.random.normal(0.0,1.0,(nsim,))*sigma*math.sqrt(t)
    #return Bases[i]*np.exp((r-0.5*sigma*sigma)*t + randnumbers)

def optionsPortfolioReturn(Bases,sigmas,CorMat,portfolio,optype=1,k=150.0,r=0.05,t=1.2,nsim=1000):
    stockPrices = [SimulateSt(Bases,j,CorMat,portfolio,optype,k,r,t,nsim) for j in sigmas]
    pfLenAdjusted = []
    for i in range(len(stockPrices)):
        if i < len(portfolio):
            pfLenAdjusted.append(portfolio[i])
        else:
            pfLenAdjusted.append(0)
    OverallStock = np.dot(stockPrices,pfLenAdjusted)*math.exp(-r*t)
    return max(optype*(OverallStock - k),0)
    
alpha = np.array([[1,0.6],[0.6,1]])
print(optionsPortfolioReturn([100,200],[0.07,0.06],alpha,[0.5,0.5]))


    
        



