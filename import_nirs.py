import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from scipy.interpolate import interpn
import os,glob
from scipy.io import loadmat

class importing:
    def __init__(self,path):
        self.dpf = 4.999
        self.optsep = 30
        self.path = path
    
    def setupexp(self):
        self.getcoefs()
        self.pickvalues()
        return self.spec.shape[1]
 
    def pickvalues(self):
        self.nirsraw = loadmat(self.path)
        self.waves = self.nirsraw['Wavelengths']
        self.spec = self.nirsraw['Spectra']
        self.ref = self.nirsraw['Ref']
        self.spec = np.array(self.spec)-np.array(self.ref)
        r,a = self.interaction()
        self.rfl = self.adjustreso(r)
        self.atn = self.adjustreso(a)
    
    def getcoefs(self):
        dataw = pd.DataFrame(pd.read_csv('code_development\extinction_int.txt', delimiter=';'))
        self.wav = dataw['extinction_int_1']
        self.coefs1 = dataw.drop('extinction_int_1', axis=1)

        datac = pd.DataFrame(pd.read_csv('code_development\ExtCoef_OD_M_CM_HBO2_HHb_oxCCO.txt', delimiter=';'))
        self.wav2 = datac['ExtCoef_OD_M_cm_1']
        self.coefs2 = datac.drop('ExtCoef_OD_M_cm_1', axis=1)
        return self.coefs1, self.coefs2

    def interaction(self):
        reflectance = np.log10(np.abs(self.ref/self.spec))
        attenuation = np.log10(np.abs(np.nanmean(self.ref)/self.spec))
        return reflectance, attenuation
    
    def adjustreso(self, arradjust):
        # interpolate from experimental wavelengths to known coef wavelengths
        arr = np.array([np.interp(self.wav, self.waves.flatten(), arradjust[i,:]) for i in range(arradjust.shape[0])])
        return arr

