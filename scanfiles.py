import numpy as np
import pandas as pd
import glob, os
from itertools import repeat
from collections import Counter
import matplotlib.pyplot as plt

class scanning:
    def __init__(self, directory):
        if os.getcwd() != directory:
            os.chdir(directory)
        self.dir = directory
    
    def create_table(self, instrument, output_name):
        self.instrument = instrument

        filelist = self.datavis()
        nameref = filelist[0][0:25]
        name = filelist[0][0:25]

        namearr = [[] for f in repeat(None,len(filelist))]
        metrics = []
        metrics1 = [[] for f in repeat(None,len(filelist))]
        namearr1 = []

        for i in range(len(namearr)):
            for j in range(len(filelist)-1):
                name = filelist[j][0:25]
                if name == nameref:
                    namearr[i].append(name)
                    metrics1[i].append([nameref, filelist[j][26:len(filelist[j])-4]])
                    metrics.append(filelist[j][26:len(filelist[j])-4])
                    checkname = filelist[j+1][0:25]
                    if checkname != name and j < len(filelist):
                        nameref = filelist[j+1][0:25]
                        namearr1.append(nameref)
                        break

        des = Counter(metrics)
        vars1 = []; reps = []
        for k,v in des.items():
            vars1.append(k)
            reps.append(round(v/len(namearr),3))

        acceptvars = ['HR-HR(bpm)', 'StO_2', '[HBT]_Deltamua', '[HBdiff]_Deltamua',
                '[oxCCO]_Deltamua', 'RefitDCSaDB',  'ART-MEAN(mmHg)','RR-RR(rpm)']
        if self.instrument == 'cyril':
            acceptvars.remove('RefitDCSaDB')

        datadict = pd.DataFrame(np.array(reps).reshape(1,len(reps)), columns=vars1)

        reps = datadict[acceptvars].flatten().copy()
        positions = np.linspace(0,len(acceptvars), len(acceptvars))

        plotting = 'no'
        if plotting == 'yes':
            fig, ax = plt.subplots(figsize=(15,8))
            #colors = cm.terrain(positions / float(max(positions)*1.1))
            plt.bar(positions, reps, color='orange', alpha=0.7)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(bottom=False, left=False)
            ax.set_axisbelow(True)
            ax.yaxis.grid(True)
            ax.xaxis.grid(False)
            plt.xticks(positions, acceptvars)

            plt.xlabel("variables / ")
            plt.ylabel("dataset frequency / ")
            fig.tight_layout()
            plt.savefig("C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\datainfo\\datafileinfo_featselect.png", dpi=500)
            plt.show()

        met = np.concatenate(metrics1,axis=0)
        df = pd.DataFrame({'Files':np.array(met[:,0]), 'Variables':np.array(met[:,1])})
        tab = pd.crosstab(index=df['Files'], columns=df['Variables'], margins=True)
        tab.to_csv(f'C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\datainfo\\{output_name}')

    def datavis(self):
        self.filelist = []
        for file in glob.glob("*.txt"):
            if file[0] == 'B' and file[8]=='P':
                self.filelist.append(file)
        self.filelist = np.array(self.filelist)
        return self.filelist
    
    def getfiles(self, variables):
        self.cleanlist = self.grouping(variables)
        return self.cleanlist
        
    def grouping(self, types):
        self.filelist = []
        for file in glob.glob("*.txt"):
            if file[0] == 'B' and file[8]=='P':
                for m in range(len(types)):
                    typename =  types[m] + '.txt'
                    if file.count(str(typename)) == 1:    
                        self.filelist.append(file)
        self.filelist = np.array(self.filelist)
        
        return self.filelist

    def filechoice(self,filename, variable):
        files = self.datavis()
        for file in files:
            if file.count(variable) == 1 and file.count(filename) == 1:
                return 1 
            
        return 0 
                
    def getspo2(self,filename):
        files = self.datavis()
        spo2files = []
        lengths = []
        for file in files:
            if file.count('SpO2') == 1 and file.count(filename) == 1:
                spo2files.append(file)
                lengths.append(len(pd.read_csv(file)))
        arrays = [x for _,x in sorted(zip(lengths,spo2files))]
        longestspo2 = pd.read_csv(arrays[-1])
        return longestspo2

    
    def signaldata(self, filename, variable):

        name =  filename + "_" + variable + ".txt"
        data = pd.read_csv(name, delimiter=';')
        if data.ndim == 2:
            data = np.average(data,axis=1)
        return data


    def individual_data(self, filename, interval,*variables):
        total = self.grouping([*variables])
        unpack = []
        datas = []  
        for i in range(len(total)):
            if total[i][0:25] == filename:
                for j in range(len(variables)):
                    name = 'D:\\uclh_data\\' + filename + "_" + variables[j] + ".txt"
                    if variables[j] not in ('rpwrhbo', 'rpwrhhb', 'rcsthbo', 'rcsthhb'):
                        data = np.array(pd.read_csv(name, delimiter=';'))  
                        if data.ndim == 2:
                            data = np.average(data,axis=1)
                        data = np.mean(data[:(len(data)//interval)*interval].reshape(-1,interval), axis=1)
                    else: 
                        data = np.loadtxt(name)
                    unpack.append(data)
        unite = list(zip(variables, unpack))
        diction = dict(unite)
        df = pd.DataFrame(diction)
        
        return unpack,df

    def scannow(self, quantities):
        # group hbo2, hhb, cco, hbt, dbdif from each file
        self.labels = quantities
        self.grouping(self.labels.copy())
        name = self.filelist[0][0:25]
        d = [[] for f in repeat(None,int(len(self.filelist)/len(self.labels)))]
        groups = []
        # iterating over files with chromophores
        for j in range(int(len(self.filelist)/len(self.labels))):
            # iterating over list of names to group all chromophores for each patient
            for i in range(len(self.filelist)):
                if name == self.filelist[i][0:25]:
                    d[j].append((self.filelist[i]))
                    
                if len(d[j]) == len(self.labels):
                    filename = "C:\\Users\\user\mres\\data_analysis\\data\\uclh_data\\" + name + "_chromophores.txt"
                    data = pd.read_csv(d[j][0],delimiter=';')
                    if len(d[j]) != 1: 
                        datarr = []
                        
                        for m in range(len(d[j])):
                            datarr.append(pd.read_csv(d[j][m], delimiter=';'))
                       
                        data = np.hstack((datarr))
              
                    #np.savetxt(filename, data)
                    groups.append([data])
                    if i < len(self.filelist) - 1:
                        lab = self.filelist[i+1]
                        name = lab[0:25]
                      
                
        d = np.array(d, dtype='object')
        groups = np.array(groups, dtype='object')
        
        return groups