import numpy as np
import matplotlib.pyplot as plt
from prel_data import preanalysis
from scanfiles import scanning

directory  = "D:\\uclh_data"

class load:
    def __init__(self,action,*args):
        self.action = action
        self.args = args

    def loaddata(self,m):
        directory  = "D:\\uclh_data"
        self.q = self.args[m]
        f = scanning(directory)
        self.groups = f.scannow(self.args[m]) 
        flistin = f.getfiles(self.args[m])       
        return flistin
        

    def runsetup(self):
        flist = self.loaddata(0)
        flist = flist[::3]
        files = [file[0:25] for file in flist]
        return np.array(self.groups), files


    def runload(self):
        for m in range(len(self.args)):
            flist = self.loaddata(m)
            if self.args[0] == ['[HbO2]_Deltamua', '[HHb]_Deltamua', '[oxCCO]_Deltamua']:
                flist = flist[::3]
            elif(len(self.args)>1):
                flist = flist[::int(len(self.args))]
        #patient number
        for k,file in enumerate(flist):
            if self.action == 'plot' and len(self.args) > 1:
                fig, axs = plt.subplots(int(len(self.args)), 1, figsize = (5,int(len(self.args)*2)))                     
            for j in range(len(self.args)):
                print(file[0:25])
                dataset = self.groups[k,:][0]
                # truncating signal             
                if dataset.shape[0] >= 0:
                    # broadcasting data  to be plotted
                    p = preanalysis(dataset, self.q)    
                    if self.action == 'plot' and len(self.args) > 1:
                        p.plot(axs,j,k)
                    if self.action == 'plot' and len(self.args) == 1:
                        p.plot1(dataset,file[0:25])
        plt.close('all')        
        return flist, self.groups

    
if __name__=='__main__':
    quantity1 = ['[HbO2]_Deltamua', '[HHb]_Deltamua', '[oxCCO]_Deltamua']

    quantity2 = ['Pulse-Pulse(bpm)_1']

    quantity3 = ['ART-MEAN(mmHg)']

    quantity4 = ['NBP-MEAN(mmHg)']

    quantity5 = ['RefitDCSaDB']

    quantity6 = ['SpO2-O2(%)']

    quantity7 = ['HR-HR(bpm)']

    quantities = [quantity1,quantity2,quantity3,quantity4,quantity5,quantity6,quantity7]
    quantities = [['Pulse-Pulse(bpm)_1'], ['NBP-MEAN(mmHg)'], ['CO2-ET(kPa)'], ['tcpCO2-TCUT(kPa)']]
    quantities = ["ART-SYS(mmHg)","CO2-ET(kPa)", "CO2-MIN(kPa)", "RR-RR(rpm)","SpO2 l-LEFT(%)","SpO2 r-RIGHT(%)","SpO2-O2(%)","StO_2"]
    quantities = [quantity1]
    for quantity in quantities:
        p = load('plot', quantity)
        p.runload()

