import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
from scipy.stats import linregress
from scipy.integrate import trapezoid
from scipy.fft import fft,  fftfreq, fftshift
from scipy.interpolate import splrep, splev
import math

class mystats:
    def __init__(self, indataset, measure):
        self.datac = indataset
        self.name = measure
        self.cs = ['r', 'g', 'b']

    def describe(self, array):
        array = array[~np.isnan(array)]
        mu = np.mean(array)
        std = np.std(array)
        medi = np.median(array)
        #print(self.name, "  mu: %s  sigma: %s  median: %s" %(round(mu,4),round(std,4),round(medi,4)))
        colors = ['r']
        ax = plt.subplots(figsize=(8,5))
        bplot = plt.boxplot(array, vert=False,patch_artist=True, labels=[self.name])
        
      
        #ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        #ax.legend([bplot["boxes"][0]],self.name)
        imagename =  self.name  + "_boxplot.png"
        plt.savefig(imagename)
        plt.show()
        
        return [mu,std, medi]

    def ddy(self, *args):
        
        fig, axs = plt.subplots()
       
        spl = splrep(self.wave[0,:],self.spectra[0,:],k=3)
        ddy = splev(self.wave,spl,der=2)
        
        if len(args) == 1:    
            l1 = axs.plot(self.wave[0,:],-self.spectra[0,:], label='data', color='purple')
            axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax2 = axs.twinx()
            l2 = ax2.plot(self.wave[0,:], ddy[0,:], alpha=0.6, label='2nd derivative', color='aqua')
            ax2.set_ylim(-1.5e4, 1.5e4)
            axs.set_xlabel(" wavelength / nm")
            axs.set_ylabel(" photon count / ")
            ax2.set_ylabel(" $d^{2}/ dλ^{2}$")
            lns = l1 + l2
            labels = [l.get_label() for l in lns]
                
            plt.tight_layout()
            imagename = 'statgraphs/' + self.filename + "_comparing_curvature.png"
            #plt.savefig(imagename)

            plt.show()
        return ddy
    
    def regress(self, data):
        #perform regression on one dataset
        r2, p, sigma = [], [], []
       
        for i in range(0,data.shape[0]):
            for j in range(i+1,data.shape[0]):
                data1 = data[i,:]
                data1 = data1[~np.isnan(data1)]
                data2 = data[j,:]
                data2 = data2[~np.isnan(data2)]
                slope, intercept, r_value, p_value, std_err = linregress(data1, data2)
                r2.append(r_value**2); p.append(p_value); sigma.append(std_err)
        #print(r2)
        return r2, p, sigma
    
    def auc(self, array_y):
        array = array_y[~np.isnan(array_y)]
        area = np.abs(trapezoid(array))
        #print(area, "area")
        return area

    def runstats(self, types):
        statdata = self.datac
        if types == 'area':
            output = self.auc(statdata)
        if types == 'r2':
            #zip your 3 arrays togrther snd pass in here
            output = self.regress(statdata)
            output = output[0]
        if types == 'describe':
            output = self.describe(statdata)
        return output
    
    def fourier(self):
        fs = []
        
        #padding
        length = len(self.datac[:,0])
        array = list(self.datac[:,0])
        print(length)
        """while math.log(length, 2).is_integer() == False :
            array.append(0)
            length = len(array)
        empty = np.zeros((length,3))
        for i in range(empty.shape[0]-len(self.datac[:,0])):
            for j in range(empty.shape[1]):
                empty[i,j] = self.datac[i,j]
        self.datac = empty.copy()"""
        print(len(self.datac[:,0]))
        import scipy
        a = scipy.fft.next_fast_len(length, real=True)
        #print(a)
        N = a
        T = 1 /(2*N)
        for i in range(3):
            four = fft(self.datac[:,i],a)
            xf = fftfreq(N, T)
            xf = fftshift(xf)
            #plt.plot(xf, four.real, xf, four.imag)
            #fs.append(four)
        return fs
    
    #def showanim(self, dt, action):
    #    f = runanim(dt, action, self.datac, self.time, self.auc, self.regress)
    #    f.run()
    #    return plt.show()
    
    def writing(self, file):
        file = open(file, "w")
        header = (" enzyme  mean   median    std    area  R2")
        file.write(header)
        for i in range(3):
            s  = self.fetchdata(i)
            file.write("\n"+ s)
        return file.close()
    
    def fetchdata(self,enzyme):
        statdata = self.describe()
        stats = statdata[enzyme]
        aucdata = self.auc(self.datac)
        auc = aucdata[enzyme]
        r2data = self.regress(self.datac)
        r2 = r2data[enzyme]
        stri = str(stats) + "  "+ str(auc) + "   "+ str(r2)
        return stri
    
    def envelope(self):
        signal = self.ddy()
        ref = signal[0,0]
        signal -= ref
        plt.subplots(figsize=(10,8))
        plt.plot(self.wave[0,:], signal.flatten(), '-b')
        plt.show()
        
        return print("hi")

   

