import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import time
import scipy as sp
from scipy.stats import zscore
import os, glob
from load_prel_data import load
from scanfiles import scanning
from scipy import stats

#rootdirectory  = "D:\\uclh_data"
#p = scanning(rootdirectory)
#groups = p.scannow([ '[HbO2]',  '[HHb]', '[oxCCO]'])
#print(groups)

datadirectory = "C:\\Users\\user\mres\\data_analysis\\data\\uclh_data"


class pwrcst:
    def __init__(self, datadirectory):
        self.datadir = datadirectory
        self.datalist = self.openfiles()
        self.M = np.eye(2)
        self.M[0,0] = 1/np.sqrt(2); self.M[0,1] = 1/np.sqrt(2); self.M[1,0] = -1/np.sqrt(2); self.M[1,1] = 1/np.sqrt(2)

    def openfiles(self):
        
        os.chdir(self.datadir)
        filearray = []
        for file in glob.glob("*.txt"):
            filearray.append(file)

        # list of files with chromophore data
        filelist = np.array(filearray)
        return filelist
    
    def getdata(self,arrayname):
        file = 'C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\' + arrayname
        array = np.loadtxt(file)
        # array is one dataset
        data1 = np.array(array[:,0]); data2 = np.array(array[:,1]); data3 = np.array(array[:,2])  
        self.imagename = 'statgraphs/relatecbf' + "_" + ".png"
        self.time = np.arange(len(data1))
        self.datac = np.array(list(zip(data1, data2, data3)))   
        return self.datac
    
    def mybeta(self,array):
        mz = zscore(array, nan_policy='propagate')
        medi = np.nanmedian(mz)
        return medi

    def zscoring(self,t0,t1):
        data = [[], [], []]
        for i in range(3):
            array = self.datac[:,i]
            data[i] = array[t0:t1]
        zhbo = self.mybeta(data[0])
        zhhb = self.mybeta(data[1])
        zcco = self.mybeta(data[2])
        zs1 = np.array([[zhbo], [zcco]]); zs2 = np.array([[zhhb], [zcco]])
        return zs1, zs2
    
    def eventrun(self):
        # necessary step to get event times
        q1 = ['[HbO2]', '[HHb]', '[oxCCO]']; q2 = ['Pulse-Pulse(bpm)_1']
        k = load('observe', q2, q1)
        output = k.runload()
        eventout = np.array(output[0], dtype='object')
        events = eventout[:,0]; intervals = eventout[:,1]; pulses = eventout[:,2]

        # iterate over all patients
        for number in range(len(events)):
            # read chromophore data
            self.getdata(self.datalist[number])
            #calculate every event in each patient
            magnhbo = [] ; magnhhb = []; pulsedata = []
            colors = plt.cm.plasma(np.linspace(0,1,len(events[number])))

            if len(events[number]) not in (0, None,1,2):
                for j in range(len(events[number])):
                    
                    fullid1 = []; fullid2 = []; pulsedata.append(pulses[number][j])

                    id1, id2 = self.runcals(int(events[number][j]-intervals[number][j]), intervals[number][j])
                    fullid1.append(id1); fullid2.append(id2)
                    id1, id2 = self.runcals(events[number][j], intervals[number][j])
                    fullid1.append(id1); fullid2.append(id2)
                    id1, id2 = self.runcals(int(events[number][j]+intervals[number][j]), intervals[number][j])
                    fullid1.append(id1); fullid2.append(id2)

                    d1 = np.array(fullid1, dtype = 'object') ; d2 = np.array(fullid2, dtype = 'object')
                    
                    rpwrhbo = d1[:,0]; rcsthbo = d1[:,1] ; rpwrhhb = d2[:,0];  rcsthhb = d2[:,1]
                    
                    magnhbo.append(np.sqrt((rpwrhbo[1]**2 + rcsthhb[1]**2)[0]))
                    magnhhb.append(np.sqrt((rpwrhhb[1]**2 + rcsthhb[1]**2)[0]))
                    #self.plotting(rpwrhbo, rcsthbo, rpwrhhb, rcsthhb, number, j) 
                
                magnhbo = np.array(magnhbo); pulsedata = np.array(pulsedata)
                # calculations on each patient over all events
                self.lsfit(pulsedata, magnhbo, magnhhb)

                

    def lsfit(self, eventnumber, depvar1, depvar2, indepvar1, indepvar2):
        depvarname = 'rpwr / '; indepvar1name = 'rcst / '; indepvar2name = 'rcst / '
        params = []
        """idx1 = np.where(~np.isnan(np.array(depvar1[:,0] + indepvar1[:,0], dtype=np.float64)))
        depvar1 = depvar1[idx1]; indepvar1 = indepvar1[idx1]
        idx2 = np.where(~np.isnan(np.array(depvar2[:,0] + indepvar2[:,0], dtype=np.float64)))
        depvar2 = depvar2[idx2]; indepvar2 = indepvar2[idx2]"""
        if len(depvar1) > 40:

            # for feature fits, scan excel file, row is file column is variable, if value == 1 then calculate, otherwise pass

            fig, axs = plt.subplots(1,2, figsize=(10,5))
            colors = plt.cm.plasma(np.linspace(0,1,len(depvar1)))
            axs[0].scatter(depvar1, indepvar1, color=colors)
            axs[1].scatter(depvar2,indepvar2,  color=colors)
            m1, b1 = np.polyfit(list(depvar1[:,0]), list(indepvar1[:,0]), 1)
            axs[0].plot(depvar1, m1*depvar1 + b1, color='gray')
            m2, b2 = np.polyfit(list(depvar2[:,0]), list(indepvar2[:,0]), 1)
            axs[1].plot(depvar2, m2*depvar2 + b2, color='gray')
            axs[0].set_xlabel(depvarname); axs[0].set_ylabel(indepvar1name)
            axs[1].set_xlabel(depvarname); axs[1].set_ylabel(indepvar2name)
            name  = 'C:\\Users\\user\mres\\data_analysis\\data\\uclh_data\\testplots\\linear_fit_' + str(eventnumber[0:25]) + "_interval_200_pulse_rpwr.png"
            plt.tight_layout()
            plt.savefig(name, dpi = 400)
            plt.close()
            params.append([m1,b1,m2,b2])
        return params
        
    def runcals(self, event, interval):
        #perform calcs within interval
        zs1, zs2 = self.zscoring(int(event-interval/2), int(event+interval/2))
        calc1 = np.dot(self.M,zs1); calc2 = np.dot(self.M,zs2)
        data1 = [calc1[0], calc1[1]]; data2 = [calc2[0], calc2[1]]
        return data1, data2
    
    def runcalsuniversal(self, t0,t):
        #perform calcs within interval
        zs1, zs2 = self.zscoring(t0, t)
        calc1 = np.dot(self.M,zs1); calc2 = np.dot(self.M,zs2)
        data1 = [calc1[0], calc1[1]]; data2 = [calc2[0], calc2[1]]
        return data1, data2
    


    def calcuniversal(self, interval):
        # for each patient
        params = []
        pwrcstdata = []
        
        for i in range(len(self.datalist)): 
            fulldata1 = []; fulldata2 = [] ; rpwrhbodata = []; rcsthhbdata = []; rpwrhhbdata = []; rcsthbodata = []
            self.getdata(self.datalist[i])
            t0 = 0
            # convert to minutes
            self.time = np.arange(interval/60,int(len(self.datac))/60, interval/60)
             
            for t in range(interval,len(self.datac), interval):
                data1, data2 = self.runcalsuniversal(t0, t)
                t0 += interval
               
                fulldata1.append(data1); fulldata2.append(data2)
            d1 = np.array(fulldata1, dtype = 'object') ; d2 = np.array(fulldata2, dtype = 'object')
            rpwrhbo = d1[:,0]; rcsthbo = d1[:,1] ; rpwrhhb = d2[:,0];  rcsthhb = d2[:,1]
          
            #self.plotting(rpwrhbo, rcsthbo, rpwrhhb, rcsthhb, i, self.datalist[i]) 
            self.interval = interval
            self.timeplot(rpwrhbo, rcsthbo, rpwrhhb, rcsthhb, i, self.datalist[i]) 
            
          
            pwrcstdata.append([rpwrhhb, rcsthbo, rpwrhhb, rcsthhb])
            #print(f"rpwr rcst cals completed for dataset {i}")
           
            # take care, you are saving the floats not the "expanded" files
            name = 'D:\\uclh_data\\' + self.datalist[i][0:25] + "_" + 'rpwrhbo' + ".txt"
            np.savetxt(name, rpwrhbo, fmt='%s', header='rpwrhhbo')
            name = 'D:\\uclh_data\\' + self.datalist[i][0:25] + "_" + 'rpwrhhb' + ".txt"
            np.savetxt(name,rpwrhhb, fmt='%s', header='rpwrhbo')
            name = 'D:\\uclh_data\\' + self.datalist[i][0:25] + "_" + 'rcsthbo' + ".txt"
            np.savetxt(name, rcsthbo, fmt='%s', header='rcsthbo')
            name = 'D:\\uclh_data\\' + self.datalist[i][0:25] + "_" + 'rcsthhb' + ".txt"
            np.savetxt(name, rcsthhb, fmt='%s', header='rcsthhb')

            name =  "C:\\Users\\user\\mres\\data_analysis\\pwrcst\\data_varying_interval_test_week3t2\\" + self.datalist[i][0:25] + "_" + 'rpwrhbo' + "_" + str(interval)+".txt"
            np.savetxt(name, rpwrhbo, fmt='%s', header='rpwrhhbo')
            name =  "C:\\Users\\user\\mres\\data_analysis\\pwrcst\\data_varying_interval_test_week3t2\\" + self.datalist[i][0:25] + "_" + 'rpwrhhb' + "_" + str(interval)+".txt"
            np.savetxt(name,rpwrhhb, fmt='%s', header='rpwrhbo')
            name =  "C:\\Users\\user\\mres\\data_analysis\\pwrcst\\data_varying_interval_test_week3t2\\" + self.datalist[i][0:25] + "_" + 'rcsthbo' + "_" + str(interval) +".txt"
            np.savetxt(name, rcsthbo, fmt='%s', header='rcsthbo')
            name =  "C:\\Users\\user\\mres\\data_analysis\\pwrcst\\data_varying_interval_test_week3t2\\" + self.datalist[i][0:25] + "_" + 'rcsthhb' +"_" + str(interval) +".txt"
            np.savetxt(name, rcsthhb, fmt='%s', header='rcsthhb')
        return rpwrhbo, rpwrhhb, rcsthhb
    
    def timeplot(self, rpwrhbo, rcsthbo, rpwrhhb, rcsthhb, num, eventnumber):
        fig, axs = plt.subplots(2,1,figsize=(15,10))
        axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator()); axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator()); axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].set_ylabel(' $\Delta HbO_2$ coupling / '); axs[1].set_ylabel(' $\Delta HHb$ coupling / ')
        axs[0].set_xlabel(' time / mins'); axs[1].set_xlabel(' time / mins')
        axs[0].plot(self.time, rpwrhbo, marker='o', color='red', alpha=0.8, label='rpwr'); axs[0].plot(self.time, rcsthbo,  marker='o', color='orange', alpha=0.8, label='rcst')
        axs[1].plot(self.time, rpwrhhb, marker= 'o', color='blue', alpha=0.8, label='rpwr'); axs[1].plot(self.time, rcsthhb,  marker='o', color='aqua', alpha=0.8, label='rcst')
        name = 'C:\\Users\\user\mres\\data_analysis\\data\\uclh_data\\testplots\\'+ str(eventnumber[0:25]) + "_reference_"+ str(self.interval) +'.png'
        plt.tight_layout()
        axs[1].legend(); axs[0].legend()
        plt.savefig(name, dpi=400)
        return plt.close()

    def plotting(self, rpwrhbo, rcsthbo, rpwrhhb, rcsthhb, num, eventnumber):
        
        colors = plt.cm.viridis(np.linspace(0,1,len(rpwrhbo)))
        fig, axs = plt.subplots(1,2,figsize=(8,4))
        name = "C:\\Users\\user\\mres\\data_analysis\\data\\pwrpct_hbocco_eventsintervals_patient_" + str(num)  + str(eventnumber) +".png" 
        
        axs[0].set_xlabel('rPWR'); axs[0].set_ylabel('rCST'); axs[1].set_xlabel('rPWR'); axs[1].set_ylabel('rCST')

        axs[0].scatter(rpwrhbo,rcsthbo, color=colors, marker='+')
        axs[1].scatter(rpwrhhb,rcsthhb, color=colors, marker='+')
        axs[0].set_xlim(-1.5,1.5), axs[0].set_ylim(-1.5,1.5); axs[1].set_xlim(-1.5,1.5), axs[1].set_ylim(-1.5,1.5)
        axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator()); axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator()); axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.tight_layout()
       # plt.show()
        name = 'C:\\Users\\user\mres\\data_analysis\\data\\uclh_data\\testplots\\' + str(eventnumber[0:25]) + '_interval_200.png'
        #plt.savefig(name, dpi=400)
        return plt.close()



"""p = pwrcst(datadirectory)

tic = time.perf_counter()
        
p.calcuniversal(200)

toc = time.perf_counter()
print("code ran in:", toc-tic, "s")"""

