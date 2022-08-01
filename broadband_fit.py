import numpy as np
import scipy as sp
from scipy.optimize import minimize, Bounds
from scipy.interpolate import splrep, splev
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

class bb:
    def __init__(self, reflectance, extinction_int, wavelengths):
        self.reflect = reflectance
        self.coefs = np.array(extinction_int)
        self.separation = 30
        self.boundaries = [
        [0.8,20,20,1,3],
        [0.8,0,0,0,0],
        [0.8,100,100,2,4]
        ]
        self.wavelengths = wavelengths
    

    def Reflectance_function(self,args):
        # mu_a = WF*mua + ln(10)*e*c(HHb) + ln(10)*e*c(HbO2)
        mua =  args[0]*self.coefs[:,2] + np.log(10)*(args[1]*self.coefs[:,0]+args[2]*self.coefs[:,1])
        mus = args[3]*(self.wavelengths*0.001)**(- args[4])
        rho = 30
        z0 = 1/mus
        D = 1/(3*(mua+mus))
        zb = (1+0.493)/(1-0.493)*2*D
        r1 = rho**2
        r2 = (z0 + 2*zb)**2 + r1
        mueff = np.sqrt(3*mua*mus)
        Reflectance =  1/(4*np.pi)*(z0*(mueff + 1/np.sqrt(r1)) * (np.exp(-mueff*np.sqrt(r1))/r1) + (z0 + 2*zb)*(mueff + 1/np.sqrt(r2))* (np.exp(-mueff*np.sqrt(r2))/r2))
        #self.Reflectance_model_diff = np.diff(Reflectance, self.order)
        spl = splrep(self.wavelengths,Reflectance,k=3)
        self.Reflectance_model_diff = splev(self.wavelengths,spl,der=self.order)
        
        Difference = -(
        self.refdif[:,np.searchsorted(self.wavelengths,self.lowerlam):np.searchsorted(self.wavelengths,self.upperlam)]
        -self.Reflectance_model_diff[np.searchsorted(self.wavelengths,self.lowerlam):np.searchsorted(self.wavelengths,self.upperlam)]
        )
        self.L= np.sum(np.square(Difference),axis=1)
        return self.L

    def run(self):
        self.refl1diff = np.diff(self.reflect)     
        self.refl2diff = np.diff(self.refl1diff)
        #output = self.apply(825, 850, refl2diff, self.boundaries[0], self.boundaries[1], self.boundaries[2])
        self.order = 2 
        self.lowerlam = 825
        self.upperlam = 850
        self.refdif = self.refl2diff
        output = minimize(self.Reflectance_function, x0 = self.boundaries[0], method='COBYLA', tol=1e-10)
        self.lowerlam = 710
        self.upperlam = 800
        self.order = 2
        output2 =  minimize(self.Reflectance_function, x0 = output.x, method='COBYLA', tol=1e-10)
        #self.visual()
        self.lowerlam = 710
        self.upperlam = 845
        self.order = 1
        self.refdif = self.refl1diff
        self.output = minimize(self.Reflectance_function, x0 = output2.x, method='COBYLA', tol=1e-10)
        self.mua, self.mus, self.SbO2 = self.calcs()
        self.dmua = self.mua - np.mean(self.mua)
        return self.dmua, self.mua, self.mus, self.SbO2

    def calcs(self):
        # last step
        mua =  self.output.x[0]*self.coefs[:,2] + np.log(10)* (self.output.x[1]*self.coefs[:,0] + self.output.x[2]*self.coefs[:,1])
        mus = self.output.x[3]*(self.wavelengths*0.001)**(- self.output.x[4])
        
        SbO2 = self.output.x[3]/(self.output.x[2]+self.output.x[1]) * 100
        """fig, axs = plt.subplots()
        plt.plot(self.wavelengths, mua, color='aqua')
        plt.ylabel('$\mu_a$ / ')
        plt.xlabel('wavelength / nm')
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.savefig('code_development\mua.png', dpi=400)
        #plt.show()
        fig, axs = plt.subplots()
        plt.plot(self.wavelengths, mus, color='lime')
        plt.ylabel('$\mu_s$ / ')
        plt.xlabel('wavelength / nm')
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.savefig('code_development\mus.png', dpi=400)
        #plt.show()"""
        return mua, mus, SbO2
    
    def getdpf(self,mua,mus,d):
        y = 0.5*((3*mus/mua)**0.5)*(1-1/(1+d*(3*mua*mus)**0.5))
        return y


