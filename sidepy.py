import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
warnings.warn('DelftStack')
warnings.warn('Do not show this message') 

class DataHelper:
    def __init__(self, instrument):
        self.instrument = instrument
    
    def save_total_feats(self,total_signal_feats, alreadysaved):
        total_signal_feats_out = pd.DataFrame(np.array(total_signal_feats).T)
        total_featnames = ['total_minimum', 'total_median', 'total_maximum']
        total_labels = []
        for variablename in alreadysaved:
            total_labels.append([f'{variablename}_{featname}' for featname in total_featnames]) 
        total_labels = np.array(np.hstack([['patient', 'session'],  np.array(total_labels).flatten()])) 
        total_signal_feats_out.columns = total_labels
        total_signal_feats_out.to_csv(f"C:\\Users\\user\\mres\\data_analysis\\total_signal_data_{self.instrument}.csv", index=False, na_rep='nan')
        total_signal_feats_out = pd.read_csv(f"C:\\Users\\user\\mres\\data_analysis\\total_signal_data_{self.instrument}.csv")
        return total_signal_feats_out
    
    def datasaver(self, tables, vals1, vals2, eventvals1, eventvals2, tablerefs, tablesess, vv, total_signal_feats, alreadysaved, combs):
        updatetable = 'no'                 
        if vv == 0 and updatetable == 'yes':
            tables.to_csv('C:\\Users\\user\\mres\\data_analysis\\datasummary.csv', index=False, sep =',')  
        totalvalues = [vals1,vals2]
        eventvalues = [eventvals1, eventvals2]
        refvals1 = np.vstack([tablerefs, tablesess])                        
        for item in range(2):
            if vv == 0 and item == 0:
                total_signal_feats.append(tablerefs); total_signal_feats.append(tablesess)
            if combs[vv][item] not in alreadysaved:
                np.save(f'C:\\Users\\user\\mres\\data_analysis\\eventtest_{combs[vv][item]}_{self.instrument}', eventvalues[item])
                np.save(f'C:\\Users\\user\\mres\\data_analysis\\reftest_{combs[vv][item]}_{self.instrument}', refvals1)
                for dataobj in range(totalvalues[item].shape[1]):
                    total_signal_feats.append(np.array(totalvalues)[item][:,dataobj])
                alreadysaved.append(combs[vv][item])
        return total_signal_feats, alreadysaved
    
    def tableadd(self, *args):
        filename, tablesess, tablerefs, tableevents,tablerefval, tablemins, tableventlocs, tabledurs, tabletotaldurs, refs, minsvals,data1,events,steps = args
        tablesess.append(filename[len(filename)-1:])
        tablerefs.append(f'{filename[5:7]}')
        tableevents.append(len(events))
        tablerefval.append(refs)
        tablemins.append(minsvals)
        tableventlocs.append(np.round(np.array(events)/3600,2))
        tabledurs.append(steps)
        tabletotaldurs.append(round(len(data1)/3600,2))
        return tablesess, tablerefs, tableevents,tablerefval, tablemins, tableventlocs, tabledurs, tabletotaldurs

    def getparamlist(self):
        systemicdata = 'no'
        if systemicdata == 'yes':
            paramlist = ['ART-MEAN(mmHg)', 'HR-HR(bpm)', 'StO_2', '[HBT]_Deltamua', '[HBdiff]_Deltamua',
            '[oxCCO]_Deltamua', 'RefitDCSaDB']
        else:
            paramlist = ['[HBT]_Deltamua', '[HBdiff]_Deltamua', '[oxCCO]_Deltamua', 'RefitDCSaDB', 'SpO2']
        
        if self.instrument == 'cyril':
            paramlist.remove('RefitDCSaDB')
        if self.instrument == 'live':
            paramlist.remove('SpO2')


        print(f"variables used: {paramlist}")
        return paramlist


    def checkcombs(self):
        # hard coding combinations 
        systemicdata = 'no'

        if systemicdata == 'yes':
            choices = [
            ['[HBT]_Deltamua','[oxCCO]_Deltamua'], ['[HBdiff]_Deltamua','[oxCCO]_Deltamua'], ['RefitDCSaDB','[oxCCO]_Deltamua'], ['RefitDCSaDB','[HBT]_Deltamua'],
            ['RefitDCSaDB','[HBdiff]_Deltamua'], ['StO_2','[HBdiff]_Deltamua'], ['StO_2','[HBT]_Deltamua'], ['StO_2','[oxCCO]_Deltamua'],
            [ 'StO_2','RefitDCSaDB'], ['[HBT]_Deltamua','ART-MEAN(mmHg)'],[ '[HBdiff]_Deltamua','ART-MEAN(mmHg)'], ['[oxCCO]_Deltamua','ART-MEAN(mmHg)'], 
            ['RefitDCSaDB','ART-MEAN(mmHg)'],[ '[HBT]_Deltamua','HR-HR(bpm)'], ['[HBdiff]_Deltamua','HR-HR(bpm)'], ['[oxCCO]_Deltamua','HR-HR(bpm)'], 
            ['RefitDCSaDB','HR-HR(bpm)']
            ]

            self.comblab = [
            'HBT - oxCCO', 'HBdiff-oxCCO', 'BFi-oxCCO', 'BFi-HBT',
            'BFi-HBdiff', 'StO2-HBdiff', 'StO2-HBT', 'StO2-oxCCO',
            'StO2-BFi', 'HBT-ARTmean', 'HBdiff-ARTmean', 'oxCCO-ARTmean', 
            'BFi-ARTmean', 'HBT-HR-HR', 'HBdiff-HR', 'oxCCO-HR', 
            'BFi-HR'
            ]

        else:
            choices = [
            ['[HBT]_Deltamua','[oxCCO]_Deltamua'], ['[HBdiff]_Deltamua','[oxCCO]_Deltamua'], ['RefitDCSaDB','[oxCCO]_Deltamua'], ['RefitDCSaDB','[HBT]_Deltamua'],
            ['RefitDCSaDB','[HBdiff]_Deltamua']]

            self.comblab = [
            'HBT - oxCCO', 'HBdiff-oxCCO', 'BFi-oxCCO', 'BFi-HBT',
            'BFi-HBdiff']

        if self.instrument == 'cyril':
            for choice, combination in list(zip(choices, self.comblab)):
                if choice.count('DCS') == 1:
                    choices.remove(choice)
                    self.comblab.remove(combination)

        return choices 
