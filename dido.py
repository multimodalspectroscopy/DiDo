import numpy as np
import matplotlib.pyplot as plt
from main_dcs import dcs
from main_nirs import nirs
import tkinter as tk
from datapanda import varian
from tkinter import filedialog
from output import clinicwindow
from pipeline_core import Pipeline
import time, os

# DiDo: Data Interpreting Diagnostic Optics

class dido:
    def __init__(self):
        self.base = os.getcwd()
        self.barstat = 0      
        self.win = clinicwindow(tk.Tk())  
        
    def run(self):
        self.win.setwin()
        self.root = self.win.getroot()
        self.root.after(2000, self.start)
        self.root.mainloop()
    
    def start(self):
        self.win.startup(self.destination_folder_initial)

    def destination_folder_initial(self):
        self.win.dest_dir_choice()
        self.target_dir = filedialog.askdirectory(title='Please choose the destination folder')
        print(self.target_dir)  
        self.processing_nirs()

    def nirs_initial(self):
        conv = False
        while conv == False:
            self.win.nirsdataload()
            nirspath = self.getnirsdata()
            if len(nirspath) > 0:
                conv = True
                break
            else:
                self.root.destroy()
                break
        nnirs = self.getnirslength(nirspath)
        self.win.nirsproc(nnirs)
        self.runnirs()
        return nirspath
    
    def dcs_initial(self):
        conv = False
        while conv == False:
            self.win.dcsdataload()
            dcspath = self.getdcsdata()
            if len(dcspath) > 0:
                break
            else:
                self.root.destroy()
                break
        ndcs = self.getdcslength(dcspath)
        self.win.dcsproc(ndcs)
        self.rundcs()
        return dcspath

    def processing_nirs(self):       
        self.nirs_initial()
        plt.plot(self.nirsdata.T)
        plt.show()
        self.win.nirsdataaccept(self.processing_dcs, self.processing_nirs)

    def processing_dcs(self):       
        self.dcs_initial()
        plt.plot(self.dcsdata)
        plt.show()
        self.win.dcsdataaccept(self.analyse, self.processing_dcs)

    
    def analyse(self):
        self.win.mlproc()
        self.runml()
        self.broadcastresult()
        self.win.showres()

    def getstatus(self):
        return self.status

    def getnirsdata(self):
        self.status = 'warming up'
        nirspath = filedialog.askopenfilename(title='Please choose NIRS data file')  
        return nirspath
    
    def getdcsdata(self):
        dcspath = filedialog.askdirectory(title='Please choose DCS data folder')
        return dcspath

    def getdcslength(self,dcspath):
        self.status = "running dcs data processing"
        self.p1 = dcs(dcspath)
        n = self.p1.load_dcsdata()
        return n

    def getnirslength(self,nirspath):
        self.status = "running nirs data processing"
        self.p = nirs(nirspath)
        n = self.p.importbroad()
        return n

    def runnirs(self):
        self.nirsdata, self.namefile, self.namebase = self.p.run(self.target_dir)
        self.sess = self.p.session
        return self.status

    def rundcs(self):
        self.dcsdata = self.p1.run(self.target_dir, self.namefile)
        return self.status
    
    def runml(self):
        os.chdir(self.base)
        boundaries = np.load('decision_boundaries.npy')
        features = np.load('selected_features.npy', allow_pickle='True')
        current_boundaries = boundaries[int(self.sess),:]
        current_features = list(features[int(self.sess)])

        files = [
            self.namebase
            ]


        self.p2 = varian('live', self.target_dir, files)
        self.p2.zfeature_acquisition()
        self.p2.rpwrrcst_acquisition()

        self.pipeline_base = Pipeline('live', 'run live', 'no', 'no')
        self.pipeline_base.learn_machine(current_features,int(self.sess))

        self.status = 'running ML magic'

        return self.status
    
    def broadcastresult(self):
        a = 'mild'
        b = '90%'
        return a,b

if __name__ == '__main__':
    p0 = dido()
    p0.run()