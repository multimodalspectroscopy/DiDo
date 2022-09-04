import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import scipy
import pickle
from sklearn.metrics import rand_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest
from scipy.stats import ttest_ind
from functionsforclusters import plotcluster, plothisto, distribution_accuracy, parameter_setup, identify_severity_regions, derive_decision_boundary, extract_distributions
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message')


def runcluster(df, outcome_data, tune_state, plotoption, saveoption, *args):
    if tune_state not in ['active', 'tuned']:
        raise Exception('please choose from active, tuned')
    if tune_state == 'active':
        if len(args) == 0:
            raise Exception('pass tune factor')
        else:
            kvalue = args[0]
    if tune_state == 'tuned':
        parameter_data = parameter_setup(df,'auto')

    keepingnoevents = 'yes'
    if keepingnoevents == 'no':
        df = df.dropna()
    session_data = []
    physio = []
    session_sig = []
    session_error = []
    classificiation_data = []
    parameters_selected_features = []
    parameters_decision_bound = []
    for kindex, session in enumerate(np.unique(df['session'])[:4]):
        # initialization
        n_cl = 2
        cls = KMeans(n_clusters=n_cl, random_state=3425)
        if tune_state == 'active':
            selector = SelectKBest(k=kvalue)
        if tune_state == 'tuned':
            selector = SelectKBest(k=parameter_data[kindex])
    

        # unpacking arrays
        X = df.loc[df['session'] == session].drop(['patient','session'],axis=1)
        patients = df.loc[df['session'] == session]['patient']
        y = np.array(outcome_data[outcome_data['STUDY_ID'].isin(np.array(patients,dtype='int'))]['CLASS'])
        outcome_unpacked = np.array(outcome_data[outcome_data['STUDY_ID'].isin(np.array(patients,dtype='int'))]['LACT/NAA'])
        
        # feature selection
        for k,pat_test in enumerate(df.loc[df['session'] == session]['patient']):
            if pat_test not in np.array(outcome_data[outcome_data['STUDY_ID'].isin(np.array(patients,dtype='int'))]['STUDY_ID']):
                X = X.drop(k, axis=0)
                patients = patients.drop(k, axis=0)
        zdata = selector.fit_transform(
                    X = X, 
                    y = y
                    )
        
        namelist = selector.get_feature_names_out()
        #parameters_selected_features.append(namelist)
        
        # dimensionality reduction
        physio.append([
            np.sum([1 for name in namelist if name.count('HBT')==1 or name.count('HBdiff')==1]), 
            np.sum([1 for name in namelist if name.count('RefitDCS')==1]), 
            np.sum([1 for name in namelist if name.count('oxCCO')==1])
            ])
        pca_2 = PCA(n_components=2)
        
        pca_2_forward = pca_2.fit(zdata)
        update_fit = 'yes'
        if update_fit == 'yes':
            pca_2_base = PCA(n_components=2)
            colsdrop = [col for col in df.columns if (col.count('before')==1 or col.count('during')==1 or col.count('after')==1)]
            zdata_live = df.drop(colsdrop,axis=1)
            
            colkeep = [
                "[HBT]_Deltamua_total_minimum","[HBT]_Deltamua_total_median","[HBT]_Deltamua_total_maximum", 
                "RefitDCSaDB_total_minimum","RefitDCSaDB_total_median"," RefitDCSaDB_total_maximum"
                ]
            colsdrop2 = [col for col in zdata_live.columns if col not in colkeep]
            zdata_live = zdata_live.drop(colsdrop2,axis=1)
            zdata_live = zdata_live.dropna(axis=1, how='all')
            
            parameters_selected_features.append(zdata_live.columns)
            
            pca_2_forlive = pca_2_base.fit(zdata_live)
            ## only run if you are using the variables (bNIRS, DCS, system) as will be available in real time, i.e. dont update when checking variable dependencies or tuning
            config_dictionary = {f'item_save_{kindex}': pca_2_forlive}
            with open('config.dictionary', 'wb') as config_dictionary_file:
                pickle.dump(config_dictionary, config_dictionary_file)
            

        pca_2_result = pca_2_forward.transform(zdata)

        # statistical significance
        session_sig.append(ttest_ind(np.array(pca_2_result)[np.argwhere(y==0)],np.array(pca_2_result)[np.argwhere(y==1)], equal_var=False)[1])

        # clustering
        labels = cls.fit_predict(pca_2_result)
        centers = cls.cluster_centers_
        outcomes_patients = np.array(outcome_data[outcome_data['STUDY_ID'].isin(np.array(patients,dtype='int'))]['CLASS'])          
        scaler = MinMaxScaler()
        pca_2_resultx = scaler.fit_transform(np.array(pca_2_result[:,0]).flatten().reshape(-1, 1))
        pca_2_resulty = scaler.fit_transform(np.array(pca_2_result[:,1]).flatten().reshape(-1, 1))
        pca_2_result = copy.deepcopy(np.hstack([pca_2_resultx,pca_2_resulty]))
        
        # classification and threshold identification
        x_db, y_db, xx, yy, Z = derive_decision_boundary(pca_2_result, labels, define_decision=False)
        m, b = np.polyfit(x_db, y_db, 1)
        

        # distribution fitting
        data0 = np.array(outcome_unpacked[np.argwhere(labels==0)])
        data1 = np.array(outcome_unpacked[np.argwhere(labels==1)])
        mu0, sigma0, mu1, sigma1 = extract_distributions(data0, data1) 

        # result plotting
        updated_labels, regions = identify_severity_regions(x_db, y_db, pca_2_result, labels, patients, np.unique(df['patient']), mu0, mu1, centers)
        parameters_decision_bound.append([m,b, regions])
        plotcluster(pca_2_result, updated_labels, patients, session, keepingnoevents, plotoption, saveoption, xx, yy, Z)
        updated_labels_inuse = updated_labels[~np.isnan(updated_labels)]
        data0 = np.array(outcome_unpacked[np.argwhere(updated_labels_inuse==0)])
        data1 = np.array(outcome_unpacked[np.argwhere(updated_labels_inuse==1)])
        if n_cl == 2:
            mu0, sigma0, mu1, sigma1 = plothisto(data0, data1, session, keepingnoevents, saveoption, plotoption) 
        plt.close()
        classificiation_data.append(updated_labels)
        updated_labels = updated_labels[~np.isnan(updated_labels)]
        
        # accuracy calculation  
        accuracy, uncertainty = distribution_accuracy(data0, data1, updated_labels, outcome_unpacked, sigma0, sigma1)
        session_data.append(accuracy)
        session_error.append(uncertainty*100)
        if tune_state == 'tuned':
            print(f'accuracy: {accuracy} %')
    
    # comparing the classifications assigned 
    classificiation_data = np.array(classificiation_data)
    flags = []
    sum_checks = np.nansum(classificiation_data,axis=0)
    for idx, patient_test in enumerate(sum_checks):
        tent = classificiation_data[:,idx]
        if patient_test != len(tent[~np.isnan(tent)]) and patient_test != 0:
            flags.append(int(1))
        else:
            flags.append(int(0))
    classificiations = pd.DataFrame(
                    np.vstack([np.unique(np.array(df['patient'],dtype=int))[:,np.newaxis].T, classificiation_data, flags]).T, 
                    columns=np.hstack(['patient', [f'session {n} classification' for n in range(1,int(classificiation_data.shape[0]+1))], 'flag'])
                    )
    #classificiations.to_csv('classification_results.csv')
    potential_seizures = np.array(classificiations.loc[classificiations['flag']==1]['patient'])
    actual_seizures = np.array([13,14,16,33,34])
    score = sum([1 for candidate in potential_seizures if candidate in actual_seizures])*100/len(actual_seizures)
    sign = ttest_ind(potential_seizures, actual_seizures, equal_var=False)[1]
    #print(score, sign, potential_seizures,actual_seizures )
    return session_data, session_sig, session_error, physio, parameters_decision_bound, parameters_selected_features
   
        