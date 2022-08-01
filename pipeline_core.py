from distutils.command.config import config
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import pickle
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from structure_setup import DataStructure
import matplotlib.ticker as ticker
from pipeline_function import runcluster
from tuning import runtuning
from filter_script import decorrelate
from functionsforclusters import filter_variables
import os
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message')

instrument = 'florence'

class Pipeline:
    def __init__(self, instrument, state, plotoption, saveoption):
        if state not in ['tune', 'examine variables', 'run simple', 'run custom', 'run live']:
            raise Exception('choose from: tune, examine variables, run simple, run custom', 'run live')
        os.chdir('C:\\Users\\user\\mres\\data_analysis')
        self.state = state
        self.instrument = instrument
        self.structure = DataStructure(self.instrument)
        self.outcome = self.structure.loadoutcome('no')
        self.plotoption = plotoption
        self.saveoption = saveoption
    
    def learn_machine(self, *args):
        if self.state != 'run live':
            self.load_data()
            if self.state == 'tune':
                plotfluctuations = 'no'
                runtuning(self.dfin, self.outcome, 5, 2, plotfluctuations, self.plotoption, self.saveoption)
            if self.state == 'examine variables':
                self.examine_variable_accuracy()
            if self.state == 'run simple':
                self.run_simple()
        else:
            self.fast_calcs(feats_keep=args[0], session=args[1])
    
    def fast_calcs(self, feats_keep, session):
        zdata = pd.read_csv(f'z_rpwr_rcst_uclh_{self.instrument}.csv')

        with open('config.dictionary', 'rb') as config_dictionary_file:
            config_dictionary = pickle.load(config_dictionary_file)
            pca_2 = config_dictionary[f'item_save_{session}']
        
        rejcols = [col for col in zdata.columns if col not in feats_keep]
        zdata = zdata.drop(rejcols,axis=1)
        pca_2_result = pca_2.transform(zdata)
        print('data loaded and preprocessed successfully')


    def load_data(self):
        zdata_uclh, zdata_patients, summary = self.structure.preprocess(f'z_rpwr_rcst_uclh_{self.instrument}.csv')
        self.dfin = copy.deepcopy(zdata_patients)

    def initialize_variables(self):
        self.combinations = [['RefitDCSaDB'], ['oxCCO'], ['HBT'], ['HBdiff'], ['RefitDCSaDB','oxCCO'], ['RefitDCSaDB', 'HBT'], ['oxCCO', 'HBT'], []]
        self.colnames = ['BFI', 'oxCCO', 'HBT', 'HBDiff', 'BFI and oxCCO', 'BFI and HBT', 'oxCCO and HBT', 'all']
        if self.instrument == 'cyril':
            for comb, name in list(zip(self.comblab, self.colnames)):
                if comb.count('DCS') == 1:
                    self.combinations.remove(comb)
                    self.colnames.remove(name)

    def examine_variable_accuracy(self):
        errors = []
        table = []
        self.initialize_variables()
        for i,combination in enumerate(self.combinations):
            df = filter_variables(self.dfin, combination)
            #df = decorrelate(df,0.9)
            a = runcluster(df, self.outcome, 'tuned', self.plotoption, self.saveoption)
            table.append(a[0])
            errors.append(a[2])
        table = pd.DataFrame(np.array(table).T, columns=self.colnames, index = [f'session {n+1}' for n in range(4)])
        table_errors = pd.DataFrame(np.array(errors).T, columns=self.colnames, index = [f'session {n+1}' for n in range(4)])
        data = table.T
        fluctuations = np.load('data_tuning_x100.npy')
        stds = np.std(fluctuations,axis=0)
        stds = pd.DataFrame(stds, columns=data.columns, index=data.index)
        error_data = table_errors.T + stds
        plotresults = 'yes'
        if plotresults == 'yes':
            fig, axs = plt.subplots(figsize=(8,5))
            data.plot(kind='bar', yerr=error_data, error_kw=dict(ecolor='k', capsize=3), cmap='winter',ax=axs)
            plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.16))
            axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            plt.ylabel('accuracy / % ')
            plt.tight_layout()
            saveplot = 'yes'
            if saveplot  == 'yes':
                plt.savefig('accuracies.png', dpi=500)
            plt.show()
        return data

    def run_simple(self,*args):
        if len(args) == 0:
            var = []
        else:
            var = args[0]
        df = filter_variables(self.dfin, var)
        #df = decorrelate(df,0.8)
        cluster_output = runcluster(df, self.outcome, 'tuned', self.plotoption, self.saveoption)
        update_params = 'yes'
        if update_params == 'yes':
            np.save('decision_boundaries.npy', cluster_output[4])
            np.save('selected_features.npy', cluster_output[5])

if __name__ == '__main__':
    plot_state = 'no'
    save_state = 'no'
    p = Pipeline(instrument, 'run simple', plot_state, save_state)
    p.learn_machine()

