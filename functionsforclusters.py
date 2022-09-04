import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from myclustering import clustering
import matplotlib.ticker as ticker
from sklearn.svm import SVC
import matplotlib
import copy
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import rand_score
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message')


def parameter_setup(df, tuning_state):
    # identify missing variables
    if tuning_state != 'auto' and tuning_state !='default':
        raise Exception('please choose tuning state auto or default')
    if tuning_state == 'auto':
        reference = list(['oxCCO', 'HBT', 'HBdiff', 'RefitDCSaDB'])
        missing_variables = []
        for var in reference:
            identified = [item for item in df.columns[2:] if item.count(var) == 0]
            if len(identified) == len(df.columns[2:]) and var not in missing_variables:
                missing_variables.append(var)
        for phase in missing_variables:
            reference.remove(phase)
        tuned_parameters = pd.read_csv('tuning//mean.csv')
        tuned_parameters = tuned_parameters.drop('Unnamed: 0',axis=1)
        tuned_parameters.columns = ['RefitDCSaDB', 'oxCCO', 'HBT', 'HBdiff', 'RefitDCSaDB and oxCCO', 'RefitDCSaDB and HBT', 'oxCCO and HBT', 'all']
        # find column
        if len(missing_variables) > 0:
            col_found = []
            for col in tuned_parameters.columns:
                if col.count(reference[0]) == 1 and len(reference) == 1:
                    col_found.append(col)
                    break
                if len(reference) == 2:
                    if col.count(reference[0]) == 1 and col.count(reference[1]) == 1:
                        col_found.append(col)
                        break
        else:
            col_found = 'all'
        parameter_data = np.array(tuned_parameters[col_found], dtype=int).flatten()
    if tuning_state == 'default':
        parameter_data = [10,10,10,10]
    return parameter_data

def filter_variables(df, *args):
    reference = list(['oxCCO', 'HBT', 'HBdiff', 'RefitDCSaDB'])
    if len(*args) != 0:
        if len(*args) == 1:
            acccols = [col for col in df.columns if (col.count('err')==0 and col.count(*args[0])>0) or col.count('session') == 1 or col.count('patient') == 1]
            item = args[0]
            reference.remove(item[0])
            rejcols = []      
            for refitem in reference:
                if not ((refitem == 'HBdiff' and item == ['HBT']) or (refitem == 'HBT' and item == ['HBdiff'])):
                    rejcols.append([col for col in acccols if col.count(refitem) == 1])
        if len(*args) == 2:
            acccols = [col for col in df.columns if (col.count('err')==0 and col.count(args[0][0])>0 and col.count(args[0][1])>0) or col.count('session') == 1 or col.count('patient') == 1]
            reference.remove(args[0][0])
            reference.remove(args[0][1])
            rejcols = []
            for refitem in reference:            
                rejcols.append([col for col in acccols if col.count(refitem) == 1])       
        df = df[acccols]      
        rejcols = np.concatenate(np.array(rejcols),axis=0)
        df = df.drop(rejcols,axis=1)
    return df




def distribution_accuracy(data0, data1, predicted_labels, lactnaascores, sigma0, sigma1):
    predicted_labels = copy.deepcopy(predicted_labels)  
    if len(data0) >= 2:
        data_ref = copy.deepcopy(data0)
    else:
        if len(data1) >= 2:
            data_ref = copy.deepcopy(data1)
    
    true_labels = np.zeros(np.array(predicted_labels).shape) 
    true_labels[(lactnaascores <= np.max(data_ref.flatten())) & (lactnaascores >= np.min(data_ref.flatten()))] = 1
    accuracy = rand_score(true_labels,predicted_labels) * 100
    if len(data0) > 2 and len(data1) > 2:
        uncertainty = sigma0 + sigma1
    else:
        if len(data0) <= 2:
            uncertainty = sigma1
        if len(data1) <= 2:
            uncertainty = sigma0       
    #print(f'accuracy: {accuracy} +- {round(uncertainty*100,4)}')
    return accuracy, uncertainty

def visualise_corrs(dfout, combination, combinations):
    savecolumns = dfout.columns
    if combination == ['RefitDCSaDB']:        
        dfout.columns = np.hstack([dfout.columns[0:2], [f'BFI{colname[11:]}' for colname in dfout.columns[2:]]])
    if combination in combinations[1:4]:        
        cnames = []
        for dflab in dfout.columns:       
            if dflab.count('mua') == 1:
                cnames.append(f'{dflab[1:len(combination[0])+1]}{dflab[len(combination)+15:len(dflab)]}')
            else:
                cnames.append(dflab)
        dfout.columns = cnames
        """plt.figure(combinations.index(combination))
        sns.heatmap(dfout.drop(['patient', 'session'], axis=1).corr()**2)
        plt.show()"""
    dfout.columns = savecolumns
    return dfout.corr()**2

def plotbar(inputtable, outputname):
    inputtable.plot.bar(cmap='winter')
    plt.legend(loc=(1.02, 0))
    plt.ylabel('accuracy % / ')
    plt.tight_layout()
    #plt.savefig(outputname, dpi=500)

def fitting(data, bins, color, heights, plotoption):
    bins = np.array(bins.flatten(),dtype='float32')
    mu, sigma = np.mean(bins), np.std(bins)
    if plotoption == 'yes':
        xx = np.linspace(0.1,0.3,500)
        best_fit_line = scipy.stats.norm.pdf(xx, mu, sigma)
        if len(heights[heights > 0]) > 1:
            plt.plot(xx, best_fit_line, color=color)   
    else:
        best_fit_line = []
    return best_fit_line, mu, sigma

def plotcluster(pca_2_result, labels, patients, session, keepingnoevents, plotoption, saveoption, xx, yy, Z):
    pca_2_resultx, pca_2_resulty = pca_2_result[:,0], pca_2_result[:,1]
    fig, axs = plt.subplots(figsize=(6,5.5))
    if plotoption == 'yes':
        # plotting contour
        plt.contour(xx, yy, Z,  cmap='winter')  
        font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
        matplotlib.rc('font', **font)
        clustering(pca_2_resultx, pca_2_resulty, patients, labels, axs)
    if saveoption == 'yes':
        plt.savefig(f'cluster_{session}_{keepingnoevents}.png', dpi=500)
    if plotoption == 'yes':
        plt.show()
    else:
        plt.close()

def derive_decision_boundary(pca_2_result, labels, define_decision):
    # optimizing model
    if define_decision == True:
        dec_bound_opt = {'C': [ 10, 100, 1000],
                'gamma': [1, 0.1, 0.01, 0.001, 0.0001],
                'kernel': ['linear', 'rbf']}
        dec_grid = GridSearchCV(SVC(), dec_bound_opt, refit = True, verbose = 1, cv=3)
        dec_grid.fit(pca_2_result, labels)
        model = dec_grid.best_estimator_
        print('model used: ', model)
    else:
        model = SVC(C=10, gamma=1, kernel='linear')
    
    # fitting optimal model
    clf = model
    X = pca_2_result
    clf.fit(X=X, y=labels)
    xx, yy = np.meshgrid(np.arange(X[:, 0].min() - 1, X[:, 0].max() + 1, 0.01),np.arange(X[:, 1].min() - 1, X[:, 1].max() + 1, 0.01))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)  

    # extracting decision boundary
    m0 = np.diff(Z, axis=1) != 0
    m1 = np.diff(Z, axis=0) != 0
    x_db = np.concatenate((xx[:, 1:][m0], xx[1:, :][m1]))
    y_db = np.concatenate((yy[:, 1:][m0], yy[1:, :][m1]))
    
    return x_db, y_db, xx, yy, Z

def identify_severity_regions(x_boundary, y_boundary, pca_2_result, labels, patients_session, patients_all, distr0, distr1, centres):
    pca_regions = np.ones(labels.shape)
    pca_2_resultx, pca_2_resulty = pca_2_result[:,0], pca_2_result[:,1]
    decision_position = []
    # find position relative to decision boundary - for every subject
    threshold = 1e-2
    cluster_regions = []
    centre = centres[0]
    ycoord_indices = np.argwhere(np.abs(centre[1] - y_boundary) < 1e-2).flatten()
    if centre[0] < x_boundary[ycoord_indices[0]]:
        cluster_regions = [0,1]
    else:
        cluster_regions = [1,0]

    for case in range(pca_2_result.shape[0]):
        conv = False
        while conv == False:
            ycoord_indices = np.argwhere(np.abs(pca_2_resulty[case] - y_boundary) < threshold).flatten()
            if len(ycoord_indices) > 0:
                threshold = 1e-2
                if pca_2_resultx[case] <= x_boundary[ycoord_indices[0]]:
                    # LHS
                    decision_position.append(0)
                else:
                    # RHS
                    decision_position.append(1)
                conv = True
                break
            else:
                threshold *= 2
    decision_position = np.array(decision_position)
    # define left or right position for each cluster, 0 for left, 1 for right
    classifications_0 = np.mean(decision_position[np.array(np.argwhere(labels==0).flatten(),dtype=int)])
    classifications_1 = np.mean(decision_position[np.array(np.argwhere(labels==1).flatten(),dtype=int)])
    if classifications_0 >= 0.5 and classifications_1 <= 0.5:
        cl0 = 1
        cl1 = 0
    if classifications_0 <= 0.5 and classifications_1 >= 0.5:
        cl0 = 0
        cl1 = 1
    pca_regions[np.array(np.argwhere(labels==0).flatten(),dtype=int)] = cl0
    pca_regions[np.array(np.argwhere(labels==1).flatten(),dtype=int)] = cl1

    # saving the positions of the patients
    pca_classification = pd.DataFrame(np.array(np.vstack([np.array(patients_session),np.array(labels), pca_regions])).T, columns=['patient', 'cluster', 'position'])
    classification_general = []
    for patient in patients_all:
        if len(np.array(pca_classification.loc[pca_classification['patient']==patient])) > 0:
             # relate the position to the severity
            if distr0 <= distr1:
                if np.array(pca_classification.loc[pca_classification['patient']==patient]['cluster']).flatten()[0] == 0:
                    classification_general.append(0)
                if np.array(pca_classification.loc[pca_classification['patient']==patient]['cluster']).flatten()[0] == 1:
                    classification_general.append(1)
            else:
                if np.array(pca_classification.loc[pca_classification['patient']==patient]['cluster']).flatten()[0] == 0:
                    classification_general.append(1)
                if np.array(pca_classification.loc[pca_classification['patient']==patient]['cluster']).flatten()[0] == 1:
                    classification_general.append(0)
        else:
            classification_general.append(np.nan)
    classification_general = np.array(classification_general)
    return classification_general, cluster_regions

def extract_distributions(data0, data1, *args):
    n1, bins1, _ = plt.hist(data0, color='blue', density=True, label='cluster 0')
    n2, bins2, _ = plt.hist(data1, color='lime', alpha=0.6, density=True , label='cluster 1') 
    best_fit_line, mu0, sigma0 = fitting(data0, bins1, color='aqua', heights=n1, plotoption='no')
    best_fit_line, mu1, sigma1 = fitting(data1, bins2, color='green', heights=n2, plotoption='no')  
    if len(args) > 0:
        data2 = args[0]
        n3, bins3, _ = plt.hist(data2, color='red', alpha=0.6, density=True , label='cluster 2')
        best_fit_line, mu2, sigma2 = fitting(data2, bins3, color='orange', heights=n3, plotoption='no')
        plt.close()
        return mu0, sigma0, mu1, sigma1, mu2, sigma2
    else:
        plt.close()
        return mu0, sigma0, mu1, sigma1

def plothisto(data0, data1, session, keepingnoevents, saveoption, plotoption, *args):
    if saveoption == 'yes' and plotoption != 'yes':
        raise Warning('modify plot option')
    
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}
    matplotlib.rc('font', **font)
    fig, axs = plt.subplots(figsize=(6,2))
    n1, bins1, _ = plt.hist(data0, color='blue', density=True, label='cluster 0')
    n2, bins2, _ = plt.hist(data1, color='lime', alpha=0.6, density=True , label='cluster 1') 
    best_fit_line, mu0, sigma0 = fitting(data0, bins1, color='aqua', heights=n1, plotoption=plotoption)
    best_fit_line, mu1, sigma1 = fitting(data1, bins2, color='green', heights=n2, plotoption=plotoption)
    if len(args) > 0:
        data2 = args[0]
        n3 = args[1]
        bins3 = args[2]
        best_fit_line, mu2, sigma2 = fitting(data2, bins3, color='orange', heights=n3, plotoption=plotoption)
    if plotoption == 'yes':
        plt.ylabel('Amplitude / ')
        plt.xlabel('LACT/NAA score / ')
        plt.legend()
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.tight_layout()
        if saveoption == 'yes':
            plt.savefig(f'threshold_{session}_{keepingnoevents}.png', dpi=500)
        plt.show()
    else:
        plt.close()
    return mu0, sigma0, mu1, sigma1      



