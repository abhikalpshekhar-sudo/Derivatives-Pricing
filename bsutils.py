from enum import Enum
import math
import numpy as np
from scipy.stats import norm

class OPTIONTYPE(Enum):
    CALL = 0
    PUT  = 1
    FORWARD = 2

def bsprice(spot, strike, vol, time, optype, riskfreerate,borrow=0.0 ):
	df = math.exp(-riskfreerate*time)
	return bspriceundiscounted(spot*math.exp((riskfreerate-borrow)*time),strike, vol, time, optype )*df

def bspriceundiscounted(fwd, strike, vol, time, optype ):
	if type==OPTIONTYPE.FORWARD:  return fwd-strike
	if abs(time)<1e-6:
		cp = max([fwd - strike,0.])
	else:
		sqrtT = math.sqrt(time)
		d1 = math.log(fwd/strike)/(vol*sqrtT) + 0.5*vol*sqrtT
		d2  = d1-vol*sqrtT
		cp = fwd*norm.cdf(d1)-strike*norm.cdf(d2)
	
	if optype==OPTIONTYPE.CALL: return cp
	if optype==OPTIONTYPE.PUT:  return cp - (fwd-strike)

def lognormal_vega(fwd,strike,sigma,time, optype ):
	if abs(time)<1e-6: return 0.0001
	if fwd<0. or strike <0. or sigma<0 : raise Exception("fwd/strike/vol cannot be negative")
	sqrtT = math.sqrt(time)
	d1 = math.log(fwd/strike)/(sigma*sqrtT) + 0.5*sigma*sqrtT
	vega = fwd * sqrtT * norm.pdf(d1)
	return vega

def blackLnImpliedVol(price,fwd, strike, time, optype):
	#price -> undiscounted option price
	niter = 100
	tol = 1e-7
	guess = 10.0/10000.0
	for i in range(niter):
		p0 = bspriceundiscounted(fwd, strike, guess, time, optype)
		if abs(p0-price)<tol:
			return guess
		vega = lognormal_vega(fwd, strike, guess, time, optype)
		guess = max([guess - (p0-price)/vega,0.0001])
	return guess


