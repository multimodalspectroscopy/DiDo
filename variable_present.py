import numpy as np
import pandas as pd
import os

os.chdir("C:\\Users\\user\\mres\\data_analysis")
data = pd.read_csv('variabletable.csv', delimiter=';')
cols = data.columns[1:len(data.columns)-2]
cols = np.array(cols)[2:]
pops = np.array(data.loc[data['Files']=='All'])[0]
pops = pops[1:len(pops)-2][2:]
print(pops)
print(len(pops), len(cols))

#pops = pops.reshape(3,int(len(pops)/3))
#cols = cols.reshape(3,int(len(cols)/3))

comp = np.vstack([cols,pops]).T
comp = comp.reshape(3,24,2)
comp = np.concatenate(comp,axis=1)
print(comp.shape)
#print(comp)

dataf = pd.DataFrame(comp)


#print(dataf)

dataf.to_csv('dissertation_variableprep.csv')
