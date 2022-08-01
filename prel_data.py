import matplotlib.pyplot as plt
import numpy as np
from windows import mystats
import matplotlib.ticker as ticker
import glob, os
from itertools import repeat
import matplotlib


#os.chdir("C:\\Users\\user\\mres\\data")
#filelist = []
#for file in glob.glob("*.txt"):
#    filelist.append(file)

#dataset = np.array(filelist)

#dataset = 'C:\\Users\\user\\mres\\data\\BS007_Session1_SynchronizedRaw.txt'

class preanalysis:
    """"
    This class is aimed to facilitate the preliminary analysis of the dataset

    """
    def __init__(self, inputdata, quantity):
        # inputdata is the whole csv
        self.inputarray = inputdata
        #self.inputarray = np.loadtxt(inputdata)
        self.name = quantity[0]
        self.q = quantity
        self.data = np.array(self.inputarray)
        self.statdata = self.inputarray


    def stats(self,*args):
        """
        area, r2, describe
        """
        #d is an array of lists
        statype = args[0]
        datarr = []
       
        for i in range(len(self.q)):
            fd = np.argwhere(str(self.q[i]) == np.array(self.q))
    
            if statype == 'r2':
                datarr.append(self.statdata[:,fd])
            else:
                if len(self.q)==1:
                    data = np.array(self.statdata[:]).flatten()
                else:
                    data = self.statdata[:,fd]
                s = mystats(data, self.q[i])
                output = s.runstats(statype)
        if statype == 'r2':
            data = np.array(list(zip(datarr)))
            s = mystats(data, self.q[i])
            output = s.runstats(statype)
        return output
   
    
    def sety(self):
        
        self.ylabel = None
        if self.q[0] in ('[HbO2]_Deltamua','[HHb]_Deltamua', '[oxCCO]_Deltamua', '[HBT]_Deltamua', '[HBdiff]_Deltamua'):
            self.ylabel =  'Hb Concentration / μM'
        if self.q[0] in ( 'BFIch_1', 'BFIch_2', 'BFIch_3', 'BFIch_4', 'BFi'):
            self.ylabel = 'BFI / cm $s^{-1}$'
        if self.q[0] in ( 'Pulse-Pulse(bpm)_1', 'Pulse-Pulse(bpm)_2'):
            self.ylabel = 'Pulse-Pulse / bpm'
        if self.q[0] in ('ABP-DIA(mmHg)', 'ABP-MEAN(mmHg)', 'ABP-SYS(mmHg)'):
            self.ylabel = 'ABP / mmHg'
        if self.q[0] == 'SpO2-O2(%)':
            self.ylabel = self.q[0]
        if self.q[0] == 'CO2-MIN(kPa)':
            self.ylabel = 'CO2-MIN / kPa'
        if self.q[0] == 'RR-RR(rpm)':
            self.ylabel = 'RR-RR / rpm'
        if self.q[0] == 'ART-MEAN(mmHg)':
            self.ylabel = 'ART-MEAN(mmHg)'
        if self.q[0] == 'NBP-MEAN(mmHg)':
            self.ylabel = 'NBP-MEAN(mmHg)'
        if self.q[0] == 'DeltaSpO2-DIFF(%)':
            self.ylabel = 'DeltaSpO2-DIFF(%)'
        if self.q[0] == 'RefitDCSaDB':
            self.ylabel = 'BFI / $cm^2 s^-1$'
        return self.ylabel

    def findq(self, q):
        self.sety(q)
        if q == '[HbO2]':
            found = 0
        if q == '[HHb]':
            found = 1
        if q == '[oxCCO]':
            found = 2
        return found
    
    def broadcast(self, q):
        a = self.findq(q)
        data = self.statdata[a]
        return data

    def plot(self, axs, m, patientn):
        """
        ['Time', '[HbO2]', '[HHb]', '[oxCCO]', 'Start', '[HBT]', '[HBdiff]',
       'BFIch_1', 'BFIch_2', 'BFIch_3', 'BFIch_4', 'BFi', 'ABP-DIA(mmHg)',
       'SpO2-O2(%)', 'ABP-MEAN(mmHg)', 'ABP-SYS(mmHg)', 'Pulse-Pulse(bpm)',
       'Pulse-Pulse(bpm)_2']
        """
       # fig, ax = plt.subplots()
        d = [[] for i in repeat(None, len(self.q))]; labels = []; idc = []
        argarr = np.array(self.q)
        
        for i in range(len(self.q)):
            found = np.argwhere(self.q[i] == argarr)
            if len(self.q) == 1:
                i = 0; found = 0
            self.sety()
            data = self.data[:,found].flatten()
            d[i].append(data)
            idc.append(found)
            labels.append(self.q[i])
        
        d = np.array(d)
        cfs = np.ones(len(d))
        if self.q[0] == '[HbO2]_Deltamua':
            self.cs = ['red', 'blue', 'green']
        else:
            self.cs = plt.cm.jet(np.linspace(0,1,20))

        for i in range(int(min(idc)),int(max(idc))+1):
            if self.q[0] == '[HbO2]_Deltamua':
                axs[m].plot(cfs[i] * d[i,0], label = labels[i], color=self.cs[i])
            else:
                axs[m].plot(cfs[i] * d[i,0], label = labels[i], color=self.cs[np.random.randint(1,20)])
  
        axs[m].legend()
        
        axs[m].set_ylabel(self.sety())
        axs[m].set_xlabel("Time / s")

        axs[m].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[m].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        #plt.legend()
        fn = "C:\\Users\\user\\mres\\data_analysis\\uclh_data\\variableplots\\" + str(self.name) + "_patient_" + str(patientn) + ".png"
        
        
        return fn

    def plot1(self, dataset,patientn):
        """
        ['Time', '[HbO2]', '[HHb]', '[oxCCO]', 'Start', '[HBT]', '[HBdiff]',
       'BFIch_1', 'BFIch_2', 'BFIch_3', 'BFIch_4', 'BFi', 'ABP-DIA(mmHg)',
       'SpO2-O2(%)', 'ABP-MEAN(mmHg)', 'ABP-SYS(mmHg)', 'Pulse-Pulse(bpm)',
       'Pulse-Pulse(bpm)_2']
        """
        fig, ax = plt.subplots()
        
        labels = []; idc = []
        argarr = np.array(self.q)
      
        d = dataset
        d = np.array(d)
   
        cfs = np.ones(len(d))
        if self.q[0] == '[HbO2]_Deltamua':
            self.cs = ['red', 'blue', 'green']
            labels = ['HbO2', 'HHb', 'oxCCO']
        colors = plt.cm.jet(np.linspace(0,1,10))

        for i in range(d.shape[1]):
            if self.q[0] == '[HbO2]_Deltamua':
                if i == 2:
                    ax2 = ax.twinx()
                    ax2.plot(d[:,i], color=self.cs[i], alpha=1)
                    ax2.set_ylabel('Concentration oxCCO / / μM', color=self.cs[i])
                    
                else:
                    ax.plot(d[:,i], color=self.cs[i], alpha=1)
            else:
                plt.plot(d[:,i], color=colors[np.random.randint(0,len(colors))], label=f'{self.q} {i}')
        #plt.legend()
        #plt.ylabel(' z signal /  ')
       # plt.ylabel('Concentration / μM')
        markers = [  '--', '--', '--']
        forcecolors = ['red', 'blue', 'green']   
        columns = [plt.plot([], [], markers[m], color=forcecolors[m])[0] for m in range(0,3)]

        ax2.legend(handles=columns, labels=[ '$HbO_{2}$', 'HHb', 'oxCCO'], loc='best')
        ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        checklab = self.sety()
        if checklab != None:
            ax.set_ylabel(self.ylabel)
        else:
            ax.set_ylabel(self.q)
       
        plt.xlabel("Time / s")

        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.xlim(0,7000)
        #plt.legend()
        # fn = "C:\\Users\\user\\mres\\data\\" + str(self.name) + "_patient_" + str(patientn) + ".png"
        filename = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\variableplots\\"+str(self.name)+"_"+str(patientn)+".png"
        #filename = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\testplots\\chromophoreplots\\chromoph"+"_"+str(patientn)+".png"
        plt.tight_layout()
        plt.savefig(filename, dpi=400)
        
        #plt.show()
        #plt.close()
        return filename

    


#plt.tight_layout()
#plt.savefig(fn)
#plt.show()