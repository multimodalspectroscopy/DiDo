import numpy as np
import matplotlib.pyplot as plt
from import_nirs import importing
from broadband_fit import bb
import time
import matplotlib.ticker as ticker
from tkinter import filedialog
import copy, os
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message') 


class SideProc:
    def __init__(self,datatype):
        self.datatype = datatype

    def savedata(self, data, filepath, target_dir, *args):
        self.data_save = data
        self.path = filepath
        if self.datatype == 'nirs':
            chromophores = ['[HbO2]_Deltamua', '[HHb]_Deltamua', '[oxCCO]_Deltamua', '[HBT]_Deltamua', '[HBdiff]_Deltamua']
            id_loc = self.path.index('B')
            id_out = self.path[id_loc:int(id_loc+7)]
            session_id = self.path[int(id_loc+8):int(id_loc+9)]     
            self.session = copy.deepcopy(session_id)
            name_root = f'{target_dir}\\{id_out}_RTSynchedSession{session_id}_'
            name_file = f'{id_out}_RTSynchedSession{session_id}'
        if self.datatype == 'dcs':
            chromophores = ['RefitDCSaDB']
            name_root = args[0]
            name_file = ''
         
        for i,chromo in enumerate(chromophores):
            np.savetxt(f'{name_root}{chromo}.txt', self.data_save[i,:])
           
            
        return name_root, name_file