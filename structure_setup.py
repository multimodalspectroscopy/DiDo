import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os 
from sklearn.preprocessing import StandardScaler

class DataStructure:
    def __init__(self, instrument):
        self.instrument = instrument

    def loadoutcome(self,plotoption):
        xls = pd.ExcelFile("NEOLIGHT_COHORT.xlsx")
        outcomes_uclh = pd.read_excel(xls, "HIE cohort")
        outcomes_uclh['STUDY_ID'] = [patient[4:] for patient in outcomes_uclh['STUDY_ID']]
        outcomes_class = []
        precols = list(outcomes_uclh.columns)
        precols.append('CLASS')
        for score in outcomes_uclh['LACT/NAA']:
            score = float(score)
            if score < 0.18:
                outcomes_class.append(0)
            if score >= 0.18:
                outcomes_class.append(1)
            if np.isnan(score) == True:
                outcomes_class.append(np.nan)
        outcome_data = pd.DataFrame(np.hstack([np.array(outcomes_uclh), np.array(outcomes_class)[:,np.newaxis]]), columns=precols).drop(index=0, axis=0)
        outcomes = np.array(outcome_data['CLASS'])
        outcome_data = outcome_data.astype({'STUDY_ID':'int'}).dropna()
        data = list(outcome_data['LACT/NAA'])
        if plotoption == 'yes':
            fig, axs = plt.subplots(figsize=(8,5))
            plt.hist(data, color='olivedrab', bins=8)
            plt.ylabel('population / ')
            plt.xlabel('LACT/NAA score / ')
            axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            #plt.savefig('outcomes.png', dpi=400)
        for classif in [0,1]:
            #print(f'class {classif} was found to be:')
            info = outcome_data.loc[outcome_data['CLASS']==classif]
            infogroup = np.array(info['LACT/NAA'])
            #print(infogroup)
            #print(f'with mean {round(np.mean(infogroup),3)} and sigma {round(np.std(infogroup),3)}')
        return outcome_data

    def unpack(self, data, patients, sessions):
        modified = []
        eventtracker = []
        patientid = []
        sessionid = []
        for patient,session,item in zip(patients,sessions,data):
            item = np.array(item)
            #print(item.shape)
            if item.ndim == 1:
                item = np.empty([3,6])
                item[:] = np.nan
                modified.append(item)
                eventtracker.append(np.nan)
                patientid.append(patient); sessionid.append(session)
            else:
                if item.ndim == 3:
                    for dime in range(item.shape[0]):
                        modified.append(item[dime,:])
                        eventtracker.append(dime)
                        patientid.append(patient); sessionid.append(session)
                if item.ndim == 2:
                    modified.append(item)
                    eventtracker.append(0)
                    patientid.append(patient); sessionid.append(session)
        modified = np.array(modified)
        eventtracker = np.array(eventtracker)
        output = np.concatenate(modified.T, axis=0)
        return output, patientid, sessionid

    def name_initialization(self):
        features = ['min', 'medi', 'max', 'auc', 'rt', 'duration']
        domains = ['before', 'during', 'after']
        pnames = ['patient', 'session']
        names = []
        for feature in features:
            names.append([f'{feature}_{timing}' for timing in domains])
        names = np.hstack([pnames, np.array(names).flatten()]).flatten()
        return names

    def zload(self,files):
        # unite data
        os.chdir("C:\\Users\\user\\mres\\data_analysis")
        zdata = []; zdata_names = []
        for datafile in files:
            reffile = f'reftest_{datafile[10:]}'
            zdata_variable = datafile[10:len(datafile)-(4+len(self.instrument)+1)]
            data = np.load(datafile, allow_pickle = True)
            patientdata = np.load(reffile, allow_pickle = True)
            patients = patientdata[0,:]; sessions = patientdata[1,:]
            output, patientid, sessionid = self.unpack(data, patients, sessions)
            zdata.append(output)
            feat_names = [f'{zdata_variable}_{feat}' for feat in self.name_initialization()[2:]]
            zdata_names.append(feat_names)
        zdata = np.array(zdata, dtype='object'); zdata_names = np.array(zdata_names, dtype='object').flatten()
        zdata = np.vstack([patientid, sessionid, np.concatenate(zdata, axis=0)])
        zdata_names = np.hstack(['patient', 'session', zdata_names])
        zdatadf = pd.DataFrame(zdata.T, columns=zdata_names).astype('float')
        duplicates = [k for k in zdatadf.columns if k.count('_duration_') == 1]
        for duplicate in duplicates:
            if len(np.array(zdatadf[duplicate])[~np.isnan(np.array(zdatadf[duplicate]))]) > 0:
                zdatadf.rename({duplicate: 'duration'}, axis=1, inplace=True)
                break
        duplicates = [k for k in zdatadf.columns if k.count('_duration_') == 1]
        for duplicate in duplicates:
            zdatadf = zdatadf.drop(duplicate,axis=1)
        zdatadf.to_csv(f"zdata_plain_{self.instrument}.csv",index=False, na_rep='nan')
        zdatadf_out = pd.read_csv(f"zdata_plain_{self.instrument}.csv")
        return zdatadf_out

    def zcombine(self, zdatadf, additional_data, destfilepath):
        os.chdir('C:\\Users\\user\\mres\\data_analysis')
        total_data = copy.deepcopy(additional_data)
        zdatafull = []
        k = 0
        for patient in np.unique(total_data['patient']):
            for session in np.unique(zdatadf.loc[zdatadf['patient'] == patient]['session']):
                zdatain = zdatadf.loc[zdatadf['patient'] == patient]
                zdatain = zdatain.loc[zdatain['session'] == session]
                fulldatain = total_data.loc[total_data['patient'] == patient]
                fulldatain = fulldatain.loc[fulldatain['session'] == session]
                fulldatain = fulldatain.drop(['patient', 'session'],axis=1)
                if np.array(fulldatain).shape[0] > 0:
                    totalextended = []
                    for eventid in range(zdatain.shape[0]):
                        totalextended.append(fulldatain)
                    totalextended = np.array(totalextended)
                
                    totalextended = totalextended.reshape(totalextended.shape[0],totalextended.shape[2])
                else:
                    totalextended = np.empty((1,fulldatain.shape[1]))
                    totalextended[:] = np.nan
                
                allz = np.hstack([zdatain, totalextended])
                if k == 0:
                    zdatafull = copy.deepcopy(allz)
                if k > 0:
                    zdatafull = np.vstack([allz,saveallz])
                saveallz = copy.deepcopy(zdatafull)
                k += 1
        fullnames = np.hstack([zdatadf.columns,total_data.columns[2:]])
        zdata_uclh_df = pd.DataFrame(zdatafull[::-1], columns=fullnames)
        zdata_uclh_df.to_csv(destfilepath, index=False, na_rep='nan')
        zdata_uclh_df = pd.read_csv(destfilepath)
        return zdata_uclh_df

    def loadcorr(self):   
        timings = ['before', 'during', 'after']
        datagather = []; labelgather = []
        for axidx, timing in enumerate(timings):
            dataload = np.squeeze(np.load(f'changeshypoxia_all_{timing}_{self.instrument}.npy', allow_pickle=True))
            if dataload.ndim == 1:               
                data = dataload[np.newaxis,:]
            else:
                data = copy.deepcopy(dataload)
            df = pd.DataFrame(data, dtype='float32')
            df.columns = ['patient', 'session','SpO2', 'HBdiff', 'oxCCO', 'errspd', 'errHBdiff', 'erroxcco', 'corr SpO2HBdiff', 'corr HBdiffoxCCO', 'corr SpO2oxCCO']
            clin = df[['patient','session']]
            df = df[df.columns[2:]]
            labels = [f'{name}_{timing}' for name in np.array(df.columns)]
            labelgather.append(labels)
            datagather.append(df)
        datagather = np.concatenate(np.array(datagather),axis=1); labelgather = np.array(labelgather).flatten()
        corrdf = pd.DataFrame(np.hstack([clin, datagather]),columns=np.hstack(['patient', 'session', labelgather]))
        corrdf.to_csv(f'zdata_corr_{self.instrument}.csv',index=False, na_rep='nan')
        corrdf = pd.read_csv(f'zdata_corr_{self.instrument}.csv')
        return corrdf

    def pearsonr_pval(self, x,y):
        return pearsonr(x,y)[1]

    def preprocess(self, path):
        df = pd.read_csv(path,delimiter=',')
        df = df.dropna(axis=1, how='all')
        removert = [coll for coll in df.columns if coll.count('rt')==1 or coll.count('Unnamed') == 1]
        df = df.drop(removert, axis=1)
        scaler = StandardScaler()
        channels = list(np.array(df.columns[2:]))
        for channel in channels:
            if channel != 'duration':
                df[channel] = scaler.fit_transform(np.array(df[channel]).reshape(-1, 1))
        singular_data = []
        session_weights = []
        params = ['RefitDCSaDB' , '[HBdiff]_Deltamua', '[HBT]_Deltamua', '[oxCCO]_Deltamua']
        if self.instrument == 'cyril':
            params.remove('RefitDCSaDB')

        for session in np.unique(df['session']):
            weights = []
            for patient in np.unique(df['patient']):
                data = df.loc[df['session'] == session].loc[df['patient'] == patient]          
                if len(data) > 0:      
                    weights.append([weight/np.nanmax(data['duration']) for weight in data['duration']])
                    currentweights = [weight/np.nanmax(data['duration']) for weight in data['duration']]  
                    if len(np.array(data['duration'])[~np.isnan(np.array(data['duration']))]) > 0:     
                        data = data.drop(index=data.index[list(np.argwhere(np.isnan(currentweights).flatten()==True).flatten())],axis=0)
                        currentweights = np.array(currentweights)[~np.isnan(np.array(currentweights))]                   
                        datatoprocess = np.average(np.array(data.drop(['patient','session'],axis=1)),axis=0,weights=np.array(currentweights))               
                    else:        
                        if len(data['duration']) > 1:
                            datatoprocess = np.array(data.drop(['patient','session'],axis=1))[0,:].flatten()
                        if len(data['duration']) == 1:
                            datatoprocess = np.array(data.drop(['patient','session'],axis=1)).flatten()
                    #print(datatoprocess)
                    databroadcast = np.hstack([np.array(data['patient'])[0], np.array(data['session'])[0], datatoprocess]) 
                    singular_data.append(databroadcast)
            if len(weights) > 1:
                weights = np.concatenate(np.array(weights,dtype='object'),axis=0)
                weights[np.array(np.argwhere(np.isnan(list(weights))==True))] = 0
            else:
                weights = np.array(weights,dtype='object').flatten()
                if weights[0] == np.nan:
                    weights[0] = 0
            session_weights.append(weights)
        zdata_patients = pd.DataFrame(singular_data, columns=df.columns).sort_values(by='patient').reset_index(drop=True).dropna(axis=1, how='all').dropna(axis=0, how='all')
        zdata_patients_summary = copy.deepcopy(zdata_patients)
        for zcol in zdata_patients.columns:
            for patient, session in list(zip(zdata_patients['patient'], zdata_patients['session'])):
                refvalue = np.nanmean(np.array(zdata_patients[zcol]))
                itemtestted = zdata_patients[zcol].loc[(zdata_patients['patient']==patient) & (zdata_patients['session']==session)]
                if np.isnan(np.array(itemtestted)) == True:
                    zdata_patients[zcol].loc[(zdata_patients['patient']==patient) & (zdata_patients['session']==session)] = copy.deepcopy(refvalue)

        return df, zdata_patients, zdata_patients_summary

