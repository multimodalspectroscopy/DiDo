from logging import raiseExceptions
import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, exp, pi, absolute
from scipy.optimize import minimize

class dcsfit:
    def __init__(self, tau, g2, beta):
        self.tau = tau
        self.g2 = g2
        self.beta = beta
        self.rho = 2
        self.mua = 0.01
        self.mus = 0.30

 
    def fitting(self):
        self.g1meas = sqrt(absolute((self.g2-1)/np.nanmean(self.beta)))
        output = minimize(self.chi2, x0 = 1e-9, method='Nelder-Mead', tol=1e-10)
        if output.success == False:
            raise Warning(" did not converge ")
        db = output.x
        return db
    
    def chi2(self,x):          
        chi2 = np.nansum((self.dcsfun(x) - self.g1meas)**2)
        return chi2

    def dcsfun(self,x):
        n0 = 1; n1 = 1.4; S0 = 1; c = 3*10**10
        k02 = (2*pi*n1/785e-7)**2
        Reff = -1.440 * n1**(-2) + 0.710 * n1**(-1) + 0.668 + 0.0636 * n1
        alpha = 1
        zb =  2*(1+Reff) / (3*self.mus*(1-Reff))
        D =  c/(3*self.mus*n1)
        K   = lambda x : sqrt(c/n1/D*(self.mua + 2*alpha*x*self.tau*self.mus*k02))
        r1  = sqrt(self.rho**2 + (-1/self.mus)**2)
        r2  =  sqrt(self.rho**2 + (1/self.mus + 2*zb)**2)

        G1 = lambda x: c/n1/(4*pi*D)*(exp(-K(x)*r1)/r1 - exp(-K(x)*r2)/r2)    
        G1norm = G1(x)/G1(0)
      
        return G1norm

   
