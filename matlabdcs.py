import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from main_dcs import dcs
import os, glob
import pandas as pd

os.chdir("D:\\uclh_data")
namelist = []
files = []
        
for file in glob.glob("*"):
    if file.count('Refit') == 1 and float(file[5:7]) <= 15:
        data = np.loadtxt(file, skiprows=1)
        if (len(data[~np.isnan(data)])) != 0:
            namelist.append(file[5:])
            files.append(file)

comparison = np.loadtxt("D:\\uclh_data\\BBS_013_PDSynchedSession4_RefitDCSaDB.txt", skiprows=1)
print(len(comparison))
comparison /= np.nanmean(comparison)
plt.plot(comparison, label='UCLn', color='tomato')
tic = perf_counter()
p1 = dcs("D:\\uclh_dcs_data\\patient13session4\\DCS")
bfi = p1.run()
bfi /= np.nanmean(bfi)
plt.plot(np.arange(0,2*len(bfi),2), bfi, color='blue', label='my DCS')
toc = perf_counter()
plt.legend()
plt.show()

print(f'results \n UCL: mean: {np.nanmean(comparison)}, sigma: {np.nanstd(comparison)} \n Mine: mean: {np.nanmean(bfi)}, sigma: {np.nanstd(bfi)}')
print(f"code ran in {np.round((toc-tic)/60,3)} mins")

