import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
import os, glob
import matplotlib.ticker as ticker

#

# In this script I am plotting the effect of varying the averaging interval with every other variable. 

# The results are given by file. 

#
os.chdir("C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest")
names = []
for file in glob.glob("*.csv"):
    if file.count('r_vals') == 1:
        names.append(file[0:25])
names = list(set(names))

names.remove('BBS_014_RDSynchedSession1')
cm_subsection = np.linspace(0, 1, len(names)) 
colors = [ cm.viridis(x) for x in cm_subsection ]





def minimize():
    
    variables = ['ART-MEAN(mmHg)', 'HR-HR(bpm)', 'Pulse-Pulse(bpm)', 'RR-RR(rpm)', 'SpO2-O2(%)', '[HHb]', 
    '[HbO2]', '[oxCCO]', 'rcsthbo', 'rcsthhb', 'rpwrhbo', 'rpwrhhb', 'tcpCO2-TCUT(kPa)', 'tcpO2-TCUT(kPa)']
 

    inters = [10,50,100,  200,  300,  400, 500, 750, 1000]
    x = inters
    for var in variables:
        fig, ax = plt.subplots(figsize=(10,8))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0*2, box.width, box.height])
        for i in range(len(names)):
            cm_subsection = np.linspace(0, 1, len(inters)) 
            colors = [ cm.jet(x) for x in cm_subsection ]
            mu = []
            for k in range(len(inters)):
                name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\"+ names[i] + "_" + str(inters[k]) +"_stats.csv"
                data = pd.read_csv(name)
                
                mu.append(data[var].mean()/data[var].std()) 
            
            plt.plot(inters,mu, c=colors[i], label=names[i], marker='+')
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=3, fancybox=True, shadow=False)
            plt.xlabel("Averaging width / s")
            ylabel =  "z " + str(var) 
            plt.ylabel(ylabel)
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator()); ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\results\\"+ str(var)  + "_zscores.png"
        plt.savefig(name, dpi=400)
       # plt.show()



output = minimize()
   



