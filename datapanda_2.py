import numpy as np
import matplotlib.pyplot as plt
from scanfiles import scanning
import pandas as pd
import os, sys
import time 
from scipy.stats import zscore
from scipy.integrate import trapezoid
import copy
from lines import myaxis, rotation
from hypoxicevents import findhypoxic, extract_corr_features
from structure_setup import DataStructure
from sidepy import DataHelper
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message') 


infopath = "C:\\Users\\user\\mres\\data_analysis\\variabletable.csv"
datapath = "D:\\uclh_data"
instrument = 'live'

class varian:
    def __init__(self, instrument, datapath, *args):
        self.datapath =  datapath
        self.instrument = instrument
        self.sideprocess = DataHelper(instrument) 
        self.variables = self.sideprocess.getparamlist()
        self.structure = DataStructure(self.instrument)
        
        if self.instrument == 'florence' or self.instrument == 'cyril':
            infopath = args[0]

            self.finfo = pd.read_csv(infopath, delimiter=';')
            self.allfiles = self.finfo['Files'] 
            self.quality = np.array(self.finfo['optical_data'])
        else:
            self.allfiles = [args[0]]
            self.quality = np.ones(len(self.allfiles))
    
    def filechoice(self, n):
        counter = 0
        chosenvariables = copy.deepcopy(self.variables)
        chosenvariables.remove('SpO2')
        for variable in chosenvariables:
            if self.gg1.filechoice(self.allfiles[n], variable) == 1:
                counter += 1
        
        if counter != len(chosenvariables):
            conv = False
        else:
            conv = True
        return conv
    
    def qualitychoice(self,n):         
        if self.quality[n] == 1:
            conv = True
        else:
            conv = False
        return conv

    def printProgressBar(self,iteration):
        prefix = 'Progress:'
        suffix = 'Complete'
        fill = '█'
        percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(len(self.allfiles))))
        filledLength = int(50 * iteration // len(self.allfiles))
        bar = fill * filledLength + '-' * (50 - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = "\r")
        # Print New Line on Complete
        if iteration == len(self.allfiles): 
            print()
    
    def acquisition_setup(self):
        print('assessing the optical data quality')
        #self.printProgressBar(0)
        filearr = []  
        # find datasets that have all variables
        self.gg1 = scanning(self.datapath)
        self.has_spo2 = []
        
        if self.instrument != 'live':
            for n in range(len(self.allfiles)):
                #self.printProgressBar(n + 1)
                if self.filechoice(n) == True and self.qualitychoice(n) == True:
                    filearr.append(self.allfiles[n])  
                    self.has_spo2.append(self.gg1.filechoice(self.allfiles[n], 'SpO2'))
                            
            if len(filearr) < 4:
                raise Warning('too few points')   
            self.has_spo2 = np.array(self.has_spo2)
        else:
            self.has_spo2 = np.zeros(np.array(self.allfiles).shape)
            filearr = copy.deepcopy(self.allfiles)
        combs = self.sideprocess.checkcombs()

        print(f"using {len(self.variables)} variables and {len(filearr)} files")
        return combs, filearr

 
    def zfeature_acquisition(self):
    
        combs, filearrin = self.acquisition_setup()
        
        alreadysaved = []; total_signal_feats = []
        for vv in range(len(combs)):   
            print(f"examining the pair: {combs[vv][0]} and {combs[vv][1]}")

            vals1 = []; vals2 = [];  sessionid = []; eventvals1 = []; eventvals2 = []
            tablerefs = []; tablesess = []; tableevents = []; tabledurs = []; tableventlocs = []; tabletotaldurs = []; tablerefval = []; tablemins = []
            # iterate over files found
            changes_before = []; changes_during = []; changes_after = []
            timings = ['before', 'during', 'after']
            
            for m in range(len(filearrin)):  
                if self.instrument == 'live':
                    filearr = filearrin[0]
                else:
                    filearr = copy.deepcopy(filearrin)
                
                sessionid.append(int(filearr[m][-1]))
                if int(filearr[m][5]) == 0 :
                    patnum = filearr[m][6:7]
                elif(int(filearr[m][5]) != 0):
                    patnum = filearr[m][5:7]
           
                # get data
                gg = scanning(self.datapath)

                data1 = gg.signaldata(filearr[m],combs[vv][0])
                data2 = gg.signaldata(filearr[m],combs[vv][1])
                
                if self.has_spo2[m] == 1:
                    dataspo2 = gg.getspo2(filearr[m])
                    pulses, self.events, self.steps,minsvals,refs = findhypoxic(filearr[m], dataspo2)
                else:
                    pulses, self.events, self.steps,minsvals,refs = [], [], [], [], []
                
                gettable = 'on'
                if gettable == 'on':
                    tablesess, tablerefs, tableevents,tablerefval, tablemins, tableventlocs, tabledurs, tabletotaldurs = self.sideprocess.tableadd(filearr[m], tablesess, tablerefs, tableevents,tablerefval, tablemins, tableventlocs, tabledurs, tabletotaldurs, refs, minsvals,data1,self.events,self.steps)
                
                featout1 = self.run(data1)
                vals1.append(featout1[0])
                eventvals1.append(featout1[1])
                dataa = copy.deepcopy(self.signal)
                featout2 = self.run(data2)
                vals2.append(featout2[0])
                eventvals2.append(featout2[1])
                datab = copy.deepcopy(self.signal)
                if self.instrument == 'live':
                    dataspo2 = np.empty(np.array(dataa).shape)
                    dataspo2[:] = np.nan
                # get p values of combination for every type
                if combs[vv][0] == '[HBdiff]_Deltamua' and combs[vv][1] == '[oxCCO]_Deltamua' :
                    change_datasets = [changes_before, changes_during, changes_after]
                    for timing, change_data in zip(timings, change_datasets):
                        change_data = extract_corr_features(timing, change_data, patnum, sessionid[m], self.events, self.steps, dataspo2, dataa, datab)

            if combs[vv][0] == '[HBdiff]_Deltamua' and combs[vv][1] == '[oxCCO]_Deltamua':
                change_datasets = [changes_before, changes_during, changes_after]
                for timing, change_data in zip(timings, change_datasets):
                    np.save(f'C:\\Users\\user\\mres\\data_analysis\\changeshypoxia_all_{timing}_{self.instrument}', change_data)
            #print(tabledurs)
            tables = pd.DataFrame(list(zip(tablerefs, tablesess, tabletotaldurs, tableevents, tableventlocs, tabledurs, tablerefval, tablemins)))
            tables.columns = ['patient number', 'session', 'duration / hr', 'event numbers', 'event time / hr' , 'event duration / s', 'ref value / %', 'min value / %']
            total_signal_feats, alreadysaved = self.sideprocess.datasaver(tables, np.array(vals1), np.array(vals2), eventvals1, eventvals2, tablerefs, tablesess, vv, total_signal_feats, alreadysaved, combs)        
            
        # data acquisition - total features
        total_signal_feats_out = self.sideprocess.save_total_feats(total_signal_feats, alreadysaved)
        
        # data acquisition - z features
        files = [file for file in os.listdir("C:\\Users\\user\\mres\\data_analysis") if file.count('eventtest_')==1 and file.count(self.instrument)==1]
        zdatadf = self.structure.zload(files)
        self.zdata = self.structure.zcombine(zdatadf, total_signal_feats_out,f'zdata_uclh_{self.instrument}.csv')

        # data acquisition - z features with correlation information
        corrdf = self.structure.loadcorr()
        zzz = pd.read_csv(f'zdata_uclh_{self.instrument}.csv')
        self.zdata =  self.structure.zcombine(zzz,corrdf, f'zdata_incl_corr_{self.instrument}.csv')
        
        return self.zdata
    
    def rpwrrcst_acquisition(self):
        combinations = self.sideprocess.checkcombs()
        self.zdata = pd.read_csv(f'C:\\Users\\user\\mres\\data_analysis\\zdata_incl_corr_{self.instrument}.csv')
        R = rotation(np.pi/4)
        features = ['min','medi','max','auc','rt']
        rpwr_features_raw = []
        rcst_features_raw = []
        name_definitions = []
        timing_calibration = ['before', 'during', 'after']
        for feature in features:
            for combination in combinations:  
                for timing_index in timing_calibration:              
                    dftest2 = [
                    col for col in self.zdata.columns if (col.count(feature)==1) and 
                    (col.count(combination[0])==1 or col.count(combination[1])==1) and 
                    (col.count(timing_index)==1)
                    ]
                    zsignal = self.zdata[dftest2]

                    if np.array(zsignal).shape[1] == 2:
                        rpwr_rcst_features = np.dot(R,np.array(zsignal).T)
                        rpwr_features_raw.append(rpwr_rcst_features[0,:])
                        rcst_features_raw.append(rpwr_rcst_features[1,:])
                        name_definitions.append(f'{combination}_{feature}_{timing_index}')

        # extracting features
        rpwr_features_raw_0 = np.hstack([self.zdata[['patient','session']], np.array(rpwr_features_raw).T])
        rcst_features_raw_0 = np.hstack([self.zdata[['patient','session']], np.array(rcst_features_raw).T])   
        rpwr_features_out = pd.DataFrame(np.array(rpwr_features_raw_0), columns=np.hstack(['patient','session',[f'{nameitem}_rpwr' for nameitem in name_definitions]]))
        rej_cols = [col for col in rpwr_features_out.columns if col.count('2')==1 or col.count('1')==1]
        rpwr_features_out.drop(rej_cols, axis=1)
        rpwr_features_out.to_csv(f"C:\\Users\\user\\mres\\data_analysis\\rpwr_uclh_{self.instrument}.csv")
        rcst_features_out = pd.DataFrame(np.array(rcst_features_raw_0), columns=np.hstack(['patient','session',[f'{nameitem}_rcst' for nameitem in name_definitions]]))
        rej_cols = [col for col in rcst_features_out.columns if col.count('2')==1 or col.count('1')==1]
        rpwr_features_out.drop(rej_cols, axis=1)
        rcst_features_out.to_csv(f"C:\\Users\\user\\mres\\data_analysis\\rcst_uclh_{self.instrument}.csv")
        
        # uniting with z data
        zrpwrrcst_raw = np.hstack([self.zdata, np.array(rpwr_features_raw).T, np.array(rcst_features_raw).T])
        labeling = np.hstack([self.zdata.columns,[f'{nameitem}_rpwr' for nameitem in name_definitions], [f'{nameitem}_rcst' for nameitem in name_definitions]])
        self.zpwrcst = pd.DataFrame(zrpwrrcst_raw, columns=labeling)
        self.zpwrcst.to_csv(f"C:\\Users\\user\\mres\\data_analysis\\z_rpwr_rcst_uclh_{self.instrument}.csv")
        
        print("relative power and relative cost features successfully extracted")
        


            
    def region_overlap(self, i):
        if i < len(self.events)-1 and i>0:
            return (self.events[i]+self.steps[i] < self.events[i+1] or self.events[i]-self.steps[i] > self.events[i-1]+self.steps[i-1])
        if i == len(self.events)-1:
            return self.events[i]-self.steps[i] > self.events[i-1]+self.steps[i-1]
        if i ==0: 
            return self.events[i]+self.steps[i] < self.events[i+1]    
    
    def run(self, datasetin):
        # this is going to return the z score of 1 trunacted signal, using the interval
        length = 3600
        # accept signal
        if len(datasetin) >= length:               
            #datasetin = datasetin[:-int(datasetin.shape[0]-length)]       
            dataset = zscore(datasetin, axis=0, nan_policy='omit')        
            totalvalmax = np.nanmax(dataset)
            totalvalmin = np.nanmin(dataset)
            totalvalmedi = np.nanmedian(dataset)
            self.signal = dataset.copy()
            # events
            interval = 200
            if len(self.events) > 0:
                # 1 event               
                if len(self.events) == 1:
                    # features 
                    regions = [(-interval,0), (0,self.steps[0]), (self.steps[0], int(self.steps[0]+interval))]
                    feats = []
                    durations = [np.nan, self.steps[0], np.nan]
                    self.events[0] = int(self.events[0])
                    for duration,r in zip(durations, regions):
                        dataregion = dataset[self.events[0]+r[0]:self.events[0]+r[1]]
                   
                        if len(dataregion[~np.isnan(dataregion)]) > 0:
                            valmin = np.nanmin(dataregion)    
                            if regions.index(r) == 1:
                                valmax = np.nan
                                valmedi = np.nan
                                valrtout = self.events[0] + self.steps[0] - np.argwhere(np.abs(dataset-valmin)<1e-3)[0]
                                valrt = valrtout[0]
                            else:
                                valmax = np.nanmax(dataregion)
                                valmedi = np.nanmedian(dataregion)
                                valrt = np.nan
                            valauc = trapezoid(dataregion[~np.isnan(dataregion)])/self.steps[0]
                            
                            feats.append([valmin, valmedi, valmax, np.abs(valauc), valrt, duration])
                # multiple events
                if len(self.events) > 1:
                    dataauc = []
                    feats = []
                    for i in range(len(self.events)):
                        # features 
                        if self.events[i] != 0:
                            regions = [[-interval,0], [0,self.steps[i]], [self.steps[i], int(self.steps[i]+interval)]]
                            durations = [np.nan, self.steps[i], np.nan]
                            singlefeats = []
                            for duration,r in zip(durations,regions):        
                                #print(self.events)
                                dataregion = dataset[self.events[i]+r[0]:self.events[i]+r[1]]
                                
                                if self.region_overlap(i) == True:
                                    if len(dataregion[~np.isnan(dataregion)]) > 0:
                                        valmin = np.nanmin(dataregion)  
                                        if regions.index(r) == 1:
                                            valmax = np.nan
                                            valmedi = np.nan
                                            valrtout = self.events[i] + self.steps[i] - np.argwhere(np.abs(dataset-valmin)<1e-3)[0]
                                            valrt = valrtout[0]
                                        else:
                                            valmax = np.nanmax(dataregion)
                                            valmedi = np.nanmedian(dataregion)
                                            valrt = np.nan
                                        valauc = trapezoid(dataregion[~np.isnan(dataregion)])/self.steps[i]
                                        
                                        singlefeats.append([valmin, valmedi, valmax, np.abs(valauc), np.nan, duration])  
                                    else:
                                        singlefeats.append([np.nan, np.nan, np.nan, np.nan, np.nan,np.nan])   
                                if self.region_overlap(i) == False:
                                    singlefeats.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])           
                        feats.append(singlefeats)
            else:             
                feats = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                valauc = np.nan              
        # reject signal
        else:
            self.signal = []
            totalvalmax = np.nan ; totalvalmin = np.nan ; totalvalmedi = np.nan
            feats = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        
        totalfeats = [totalvalmin, totalvalmedi, totalvalmax]
        
        return [totalfeats, feats]

if __name__ == '__main__':
    
    tic = time.perf_counter()
    
    datapath = 'D:\\test'

    files = [
        'BBS_022_RTSynchedSession2'
    ]
    p = varian(instrument, datapath, files)

    p.zfeature_acquisition()
    p.rpwrrcst_acquisition()


    toc = time.perf_counter()
    print(f"code ran in {(toc-tic)/60} minutes")

