from importlib.metadata import files
import numpy as np
import matplotlib.pyplot as plt
from load_prel_data import load


k = load('plot', ['[HbO2]_Deltamua', '[HHb]_Deltamua', '[oxCCO]_Deltamua'])
data, files = k.runsetup()


def snr(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m/sd)


for k,file in enumerate(files):
    print(data[k,:][0].shape)
    dataset = data[k,:][0]
    val = np.diff(dataset,n=1,axis=1)
