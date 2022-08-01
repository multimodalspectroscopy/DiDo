import numpy as np
import matplotlib.pyplot as plt
from import_nirs import importing
from broadband_fit import bb
import time
import matplotlib.ticker as ticker
from tkinter import filedialog
import copy, os
from sideprocess import SideProc
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message') 



class nirs:
    def __init__(self,path):
        self.path = path
        self.helpfun = SideProc('nirs')
    
    def importbroad(self):
        self.a = importing(self.path)
        length = self.a.setupexp()

        return length

    def run(self, target_dir):

        b = bb(self.a.rfl, self.a.coefs1, self.a.wav)
        dmua, mua, mus, sbo2 = b.run()
        #print("mua", np.mean(mua), "mus", np.mean(mus))
       
        self.mua = mua.copy()
        self.mus = mus.copy()
       
        coefinv = np.linalg.pinv(self.a.coefs1)
        dpf = b.getdpf(mua,mus,self.a.optsep)
        coefbydpf = np.zeros((coefinv.shape))
        for i in range(coefinv.shape[0]):
            for m in range(coefinv.shape[1]):
                coefbydpf[i,m] = coefinv[i,m]/dpf[m]
        DeltaC_DPFdynamic=np.matmul(coefbydpf,self.a.atn.T)/(self.a.optsep)
        DeltaC_HBT = DeltaC_DPFdynamic[0,:] + DeltaC_DPFdynamic[1,:]
        DeltaC_HBdiff = DeltaC_DPFdynamic[0,:] - DeltaC_DPFdynamic[1,:]

        plt.close('all')
        cols = ['g','b','r']
        fig, axs = plt.subplots()
        for k in range(3):
            plt.plot(DeltaC_DPFdynamic[k,:],color=cols[k])
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.xlabel('wavelength / nm')
        plt.ylabel('Concentration / μΜ')
        #plt.savefig('code_development\\pythontest.jpg', dpi=400)
        #plt.show()
        plt.close()
        self.waves = self.a.wav
        self.data_save = np.vstack([DeltaC_DPFdynamic, DeltaC_HBT, DeltaC_HBdiff])
        nameout, namefile = self.helpfun.savedata(self.data_save, self.path, target_dir)
        self.session = self.helpfun.session
        return DeltaC_DPFdynamic,nameout, namefile
        

if __name__ == '__main__':
    
    filepaths = []
    times = []
    lengths = []
    data_dir = 'D:\\Raw_NIRS_Data'
    for file in os.listdir(data_dir)[:1]:
        if file.endswith(".mat") and len(os.path.join(data_dir, file)) > 40 and file.count('BBS') == 1:        
            #nirspath = filedialog.askopenfilename(title='Please choose NIRS data file')  
            p = nirs(os.path.join(data_dir, file))
            p.importbroad()
            tic = time.perf_counter()
            dconc, nn = p.run()
            toc = time.perf_counter()
            if dconc.shape[1] > 10:
                times.append(toc-tic)
                lengths.append(dconc.shape[1])
