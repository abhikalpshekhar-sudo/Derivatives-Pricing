import math
import numpy as np

class Node:
    def __init__(self,s,pu,pd):
        self.s=s
        self.pu=pu
        self.pd=pd
        self.v = 0.0

class BinTree:
    def __init__(self,optype,S0=100.0,r=0.05,sigma=0.2,T=1.2,k=100.0,nsteps=10):
        self.optype = optype
        self.S0=S0
        self.r = r
        self.sigma = sigma
        self.T=T
        self.k = k
        self.nsteps = nsteps
        self.timelines = np.linspace(0,T,self.nsteps+1)#Added +1 in self.nsteps, also why is T = 1.2 
 
    def createTree(self):
        nt = self.timelines.shape[0] 
        delt = self.timelines[1]-self.timelines[0]
        sqrtDelt = math.sqrt(delt)
        up = math.exp(self.sigma*sqrtDelt)
        dn = 1.0/up 
        self.dfs = np.ones(self.timelines.shape)*math.exp(-self.r*delt)
        self.slices = []
        pu = (math.exp(self.r*delt)-dn)/(up-dn) 
        pd = 1.0 - pu
        for i in range(nt): #time direction
            cols = []
            for j in range(i+1): #number of spots at each time step
                S = self.S0*math.pow(up,i-2*j) 
                cols.append(Node(S,pu,pd))
            self.slices.append(cols)
        
    def price(self):
        self.createTree()
        nt = self.timelines.shape[0]
        for i in range(nt-1,-1,-1):
            for j in range(i+1):
                if i==nt-1:
                    if self.optype==1: #call option
                        self.slices[i][j].v = max(self.slices[i][j].s-self.k,0.0)
                    else:
                        self.slices[i][j].v = max(-self.slices[i][j].s+self.k,0.0)
                else:
                    self.slices[i][j].v = self.dfs[i]*(self.slices[i][j].pu*self.slices[i+1][j].v +
                                                       self.slices[i][j].pd*self.slices[i+1][j+1].v)
        return self.slices[0][0].v
    
if __name__ == "__main__":
    tr = BinTree(1)
    print(f"Price Of Call Option is {tr.price()}")
    tr.optype=0    
    print(f"Price Of Put Option is {tr.price()}")            
