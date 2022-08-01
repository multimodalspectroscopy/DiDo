import numpy as np
from numpy import nanmean, nanstd
from opendcs import dcsload
from fitfundcs import dcsfit
import matplotlib.pyplot as plt
import copy, os
from sideprocess import SideProc
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message') 


class dcs:
    def __init__(self, pathin):
        self.path = pathin
        self.load_dcsdata()
        self.helpfun1 = SideProc('dcs')

    def load_dcsdata(self):
        #print("loading data")
        p = dcsload(self.path)
        data = p.united()
        data = np.array(data,dtype='object')
        self.tau = np.array(data[:,0][0])
        self.g2in = np.concatenate(data[:,1],axis=0)
        self.N2 = self.g2in.shape[0]
        
        return self.N2

    def clean(self):
        # params
        overfloats = np.argwhere((self.g2in > 1e20) | (self.g2in < 1e-20))
        self.g2in[overfloats] = np.nan
        g2clean=self.g2in.copy()

        # no extrema or mean over threshold
        tau0 = np.argwhere((self.tau<10**-1) & (self.tau>10**-2))
        meanattau = nanmean(self.g2in[:,:,tau0].reshape(self.g2in[:,:,tau0].shape[0], self.g2in[:,:,tau0].shape[1],self.g2in[:,:,tau0].shape[2]),axis=2).flatten()
        a1_rej = np.argwhere(np.abs(meanattau)>1.02).flatten()
        g2clean[a1_rej,:,:]=np.nan
        
        # beta stat calcs
        beta0 = (np.nansum(self.g2in[:,:,2:4],axis=2)-3)/3
        betamean0 = nanmean(beta0)
        betasigma0 = nanstd(beta0)

        # remove outisde threshold
        rejbeta = np.argwhere((beta0 > 0.54) & (beta0 < 0.40 ))
        g2clean[rejbeta]=np.nan
        
        # clean g1 calcs
        g1 = np.sqrt(np.absolute((self.g2in-1))/betamean0)
       # print(self.g2in[~np.isnan(self.g2in)])
        g1mean = nanmean(g1,axis=2)      
        g1sigma = nanstd(g1,axis=2)

        # remove rej g1 from g2
        g2clean[np.argwhere(nanmean(g1[:,:,tau0],axis=2) > 0.08)] = np.nan
        g2mean = nanmean(g2clean, axis=2)
        g2sigma = nanstd(g2clean, axis = 2)

        self.g2 = g2clean
        self.beta0 = beta0
        
        return self.g2

    def runfitting(self):
        # fitting      
        fv0=np.zeros((self.N2))
        beta0 = (np.nansum(self.g2,axis=1)-1)/1
        betamean0 = nanmean(beta0)
        betasigma0 = nanstd(beta0)
        self.beta0 = beta0
        BFI = np.zeros((self.N2))

        for j in range(0,self.N2):
            if np.isnan(self.g2[j,0,0])==False:   
                BFI[j] = dcsfit(self.tau,self.g2[j,0,:].flatten(),self.beta0[j,:]).fitting() 
        plt.plot(np.arange(0,2*len(BFI),2),BFI, color='blue')
        plt.savefig("C:\\Users\\user\\mres\\code_development\\bfiiii.jpg")
        self.bfi = np.array(BFI.copy())[np.newaxis,:]

    def run(self, name_out, target_dir):      
        self.clean()
        self.runfitting()
        self.helpfun1.savedata(self.bfi, self.path, name_out, target_dir)
        print("ready")
        
        return self.bfi


if __name__ == '__main__':
    

    data_dir = 'D:\\Raw_NIRS_Data'
    p = dcs(data_dir)
    p.load_dcsdata()
    p.run()
                
