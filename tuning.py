import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pipeline_function import runcluster
from functionsforclusters import filter_variables
import pandas as pd
import copy

def runtuning(dfin, outcome_data, test_density, reps, plotfluctuations, plotoption, saveoption):
    table1 = []
    j = 0
    combinations = [['RefitDCSaDB'], ['oxCCO'], ['HBT'], ['HBdiff'], ['RefitDCSaDB','oxCCO'], ['RefitDCSaDB', 'HBT'], ['oxCCO', 'HBT'], []]
    table1 = []
    for rep in range(reps):
        optimal_rep = []
        for i,combination in enumerate(combinations):
            df = filter_variables(dfin, combination)
            optimal_k, j = tuning_number_features(df, outcome_data, j, int(test_density * len(combinations) * reps), test_density, plotoption, saveoption)
            optimal_rep.append(optimal_k)          
        table1.append(optimal_rep)
    table1 = np.array(table1)
    np.save(f'data_tuning_x{reps}.npy', table1)
    if plotfluctuations == 'yes':
        plot_fluctuations(table1)   
    return table1

def tuning_number_features(df, outcome_data, j, maxlength,test_density, plotoption, saveoption):
    test_range = np.linspace(0.2,1,test_density)
    tuning_data = []
    if df.shape[1] < 30:
        test_factors = np.unique(np.array(test_range * (df.shape[1]-2), dtype=int)).flatten()
    else:
        test_factors = np.unique(np.array(test_range * 30, dtype=int)).flatten()
    tuning_data = []
    for m,tuningfactor in enumerate(test_factors):
        printProgressBar(j+m+1,maxlength)
        run_i = runcluster(df, outcome_data, 'active', plotoption, saveoption, tuningfactor)
        tuning_data.append(run_i[0])
    j += (m + 1)
    tuning_data = np.array(tuning_data)
    optimal_factors = []   
    for session in range(0,4):
        # only keep n features >= 10
        indices = np.array(np.argwhere(test_factors >= 10).flatten())
        k_ten = test_factors[indices]
        tuning_tent = tuning_data[indices,session]
        # sort by accuracy
        k_optimal = k_ten[np.argwhere(np.abs(tuning_tent-np.max(tuning_tent))<1e-5)].flatten()[0]
        optimal_factors.append(k_optimal)
    optimal_factors = np.array(optimal_factors)
    if plotoption == 'yes':
        fig, axs = plt.subplots(figsize=(8,6))
        colors = ['blue', 'red', 'purple', 'lime']
        for session in range(0,4):
            plt.plot(test_factors, tuning_data[:,session], color=colors[session], label=f'session {session+1}', marker='.')
        plt.legend()
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.ylabel('accuracy / %')
        plt.xlabel('feature inclusion / %')
        plt.show()
    return optimal_factors, j

def printProgressBar(iteration,maxlength):
    prefix = 'Progress:'
    suffix = 'Complete'
    fill = '█'
    percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(maxlength)))
    filledLength = int(50 * iteration // maxlength)
    bar = fill * filledLength + '-' * (50 - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = "\r")
    # Print New Line on Complete
    if iteration == maxlength: 
        print()


def plot_fluctuations(results):
    reps = results.shape[0]
    combinations = [['RefitDCSaDB'], ['oxCCO'], ['HBT'], ['HBdiff'], ['RefitDCSaDB','oxCCO'], ['RefitDCSaDB', 'HBT'], ['oxCCO', 'HBT'], []]
    colnames = ['BFI', 'oxCCO', 'HBT', 'HBDiff', 'BFI and oxCCO', 'BFI and HBT', 'oxCCO and HBT', 'all']
    results_std = pd.DataFrame(np.std(results, axis=0).T, columns=colnames, index=[f'session {n+1}' for n in range(4)])
    fig, axs = plt.subplots(2,4, figsize=(24,10))
    ind = [[0,0], [0,1], [0,2], [0,3], [1,0], [1,1], [1,2], [1,3]]
    for m, combination in enumerate(combinations):
        i = ind[m][0]
        j = ind[m][1]
        colors = ['yellow', 'orange', 'red', 'purple']
        axs[i,j].set_ylim(0,50)
        axs[i,j].set_title(colnames[m])
        for session in range(0,4):
            errval = np.array(results_std)[session,i]
            err = np.linspace(errval, errval, reps)
            axs[i,j].plot(np.arange(reps), results[:,m,session], color=colors[session], label=f'session {session+1}', marker='o')
            axs[i,j].errorbar(np.arange(reps), results[:,m,session],yerr=err, fmt='None', ecolor=colors[session], capsize=2)
        axs[i,j].legend()
        axs[i,j].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[i,j].set_ylabel('optimal number of features / ')
        axs[i,j].set_xlabel('repetition')
    plt.tight_layout()
    return plt.show()