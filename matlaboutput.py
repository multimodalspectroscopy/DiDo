import wave
import numpy as np
import matplotlib.pyplot as plt
from main_nirs import nirs
from plotting import plotcomp
import pandas as pd

patients = ["21-Jul-2021 213204 LWP_718 Patient 1", "06-Feb-2021 153753 BBS_017_01 Patient 1",
 "20-Mar-2021 113838 BBS_020_2 Patient 1", 
"05-May-2021 103356 BBS_025_3 Patient 1"]

patients = ["06-Feb-2021 153753 BBS_017_01 Patient 1"]

def tune(patient):
    p = nirs(f"D:\\uclh_code\\{patient}")
    p.run()
    plt.close('all')
    wavelengths = p.waves
    mymua = p.mua
    mymus = p.mus

    mua = np.loadtxt(f"D:\\uclh_code\\mua{patient}.txt",delimiter=',')
    mus = np.loadtxt(f"D:\\uclh_code\\mus{patient}.txt",delimiter=',')

    resmua = mua - mymua
    resmua2 = np.sum(np.square(resmua),axis=1)
    
    idx = np.argwhere(resmua2 == np.nanmin(resmua2)).flatten()[0]

    return wavelengths, [mymua, mymus], [mua[idx,:],mus[idx,:]]

datameans = []
datastds = []
mymeans = []
vals730mua = []
vals850mua = []
vals730mus = []
vals850mus = []

for patient in patients:
    wavelength, mydata, data = tune(patient)
    low = np.searchsorted(wavelength, 730)
    high = np.searchsorted(wavelength, 850)
    
    mymeans.append(np.mean(mydata,axis=1))
    datameans.append(np.mean(data,axis=1))
    datastds.append(np.std(data,axis=1))
    data = np.array(data); mydata = np.array(mydata)
    vals730mua.append([mydata[0,low],data[0,low]])
    vals850mua.append([mydata[0,high],data[0,high]])
    vals730mus.append([mydata[1,low],data[1,low]])
    vals850mus.append([mydata[1,high],data[1,high]])
    #plotcomp(patient,wavelength, low, high, mydata,data)


mymeans = np.array(mymeans); datameans = np.array(datameans); datastds = np.array(datastds)
vals730mua = np.array(vals730mua); vals850mua = np.array(vals850mua); vals730mus = np.array(vals730mus); vals850mus = np.array(vals850mus)

print("\n \n ------------------------------------------------------ \n")
dictionmua = dict({'patients':patients,'my means': mymeans[:,0], 'data means': datameans[:,0], 'std': datastds[:,0],
'my vals 730nm':vals730mua[:,0], 'data 730nm':vals730mua[:,1]})
dfstatmua = pd.DataFrame(dictionmua)
df = dfstatmua.round(4)
print(df)

print("\n ------------------------------------------------------ \n")

dictionmus = dict({'patients':patients,'my means': mymeans[:,1], 'data means': datameans[:,1], 'std': datastds[:,1],
'my vals 730nm':vals730mus[:,0], 'data 730nm':vals730mus[:,1]})
dfstatmua = pd.DataFrame(dictionmus)
df = dfstatmua.round(4)
print(df)

print("\n ------------------------------------------------------")
