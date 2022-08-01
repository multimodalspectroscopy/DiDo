from logging import raiseExceptions
import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, exp, pi, absolute
from scipy.optimize import minimize

tau = np.loadtxt("code_development\\tautest.txt", delimiter=',')
rho = 2
mua = 0.01
mus = 0.3


n0 = 1; n1 = 1.4; S0 = 1; c = 3*10**10
k02 = (2*pi*n1/785e-7)**2
Reff = -1.440 * n1**(-2) + 0.710 * n1**(-1) + 0.668 + 0.0636 * n1
alpha = 1
zb =  2*(1+Reff) / (3*mus*(1-Reff))
D =  c/(3*mus*n1)
    
K   = lambda x : sqrt(c/n1/D*(mua + 2*alpha*x*tau*mus*k02))
r1  = sqrt(rho**2 + (-1/mus)**2)
r2  =  sqrt(rho**2 + (1/mus + 2*zb)**2)

G1 = lambda x: c/n1/(4*pi*D)*(exp(-K(x)*r1)/r1 - exp(-K(x)*r2)/r2)
        
G1norm = lambda x: G1(x)/G1(0)

xx = np.linspace(-2e-9,-0.5e-9,100)
print(G1norm(-1e-9))