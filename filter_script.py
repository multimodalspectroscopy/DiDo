import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pipeline_function import runcluster
from functionsforclusters import filter_variables, visualise_corrs
import pandas as pd
import copy


def decorrelate(dfin, corr_threshold):
    corrdata = dfin.corr()**2
    tested = []
    df_mask = copy.deepcopy(dfin)
    for k,variable1 in enumerate(corrdata.columns[2:]):
        correlated_variable = []
        for variable2 in corrdata.columns[k:]:
            if variable1 != variable2 and corrdata[variable1][variable2] > corr_threshold and variable2 not in tested:
                correlated_variable.append(variable2)
        if len(correlated_variable) > 0:             
            weighting = 1/(len(correlated_variable)+1)
            df_mask[variable1] *= weighting
            for corr_var2 in correlated_variable:
                df_mask[corr_var2] *= weighting
        tested.append(variable1)
    return df_mask

def examine_correlations(dfin, combinations, colnames, threshold_array, reps):

    rep_data = []
    for rep in np.arange(reps):
        ref_arr = []; test_arr = []   
        for filtervals in threshold_array:
            data1, data2, sig1, sig2 = filter(dfin, filtervals, combinations, colnames)
            test_arr.append(data2)
            ref_arr.append(data1)
        test_arr2 = np.array(test_arr)
        ref_arr2 = np.array(ref_arr)
        rep_data.append(test_arr2)
    rep_data = np.array(rep_data)

    return rep_data

