from audioop import avg
import numpy as np
import glob, os
from time import perf_counter, sleep


class dcsload:
    def __init__(self,pathin):
        self.pathin = pathin
        os.chdir(self.pathin)
    
    def united(self):
        data = []
        directors = self.findfolers(self.pathin)
        
        for dirs in directors:
            datafolder = self.findfolers(dirs)
            if len(datafolder) > 0:
                os.chdir(datafolder[0])
                data.append(self.opendcsfiles(f'{self.pathin}\\{dirs}\\{datafolder[0]}'))
                os.chdir(self.pathin)
            else:
                os.chdir(self.pathin)
        return data

    def findfolers(self,pathused):
        os.chdir(pathused)
        folders = []
        [folders.append(name) for name in os.listdir(".") if os.path.isdir(name)]
        dirsort = np.sort(folders)
        return dirsort


    def opendcsfiles(self,infolderpath):
        
        namelist = []
        files = []
        
        for file in glob.glob("*"):
            namelist.append(file[5:])
            files.append(file)
    
        namelist = np.array(namelist, dtype='float32')
        filessorted = [x for _, x in sorted(zip(namelist, files))]
    
        output = []
        dims = []

        ch_index = []
        ch_intensity = []
        ch_mark = [] 
        ch_beta = []
        ch_aDb = []
        ch_g2 = []
        ch_fitted_g2 = []

        avg_index = []
        avg_intensity = []
        avg_mark = []
        avg_beta = []
        avg_aDb = []
        avg_g2 = []
        avg_fitted_g2 = []
        timings = []

        tic = perf_counter()
        for idx in range(len(filessorted)):
            file = filessorted[idx]
            data = open(f"{infolderpath}\\{file}", "rb")
            numCh = np.fromfile(data,count=1,dtype='uint32')
            numAverages = np.fromfile(data,count=1,dtype='uint32')
            numBins = np.fromfile(data,count=1,dtype='uint32')
            
            timings.append(perf_counter()-tic)
            if len(numBins) > 0:
                tau = np.fromfile(data,count=int(numBins),dtype='float64')*1e-9

                
                ch_indexdata = []
                ch_intensitydata = []
                ch_markdata = []
                ch_betadata = []
                ch_aDbdata = []
                ch_g2data = []
                ch_fitted_g2data = []


                for k in range(0,numCh[0]):
                    ch_indexdata.append(np.fromfile(data,count =1,dtype='uint32'))
                    ch_intensitydata.append(np.fromfile(data,count=1,dtype='float64'))
                    ch_markdata.append(np.fromfile(data,count=1,dtype='uint32'))
                    ch_betadata.append(np.fromfile(data,count=1,dtype='float64'))
                    ch_aDbdata.append(np.fromfile(data,count=1,dtype='float64'))
                    ch_g2data.append(np.fromfile(data,count=int(numBins),dtype='float64'))
                    ch_fitted_g2data.append(np.fromfile(data,count=int(numBins),dtype ='float64'))

                ch_index.append(ch_indexdata)
                ch_intensity.append(ch_intensitydata)
                ch_mark.append(ch_markdata)
                ch_beta.append(ch_betadata)
                ch_aDb.append(ch_aDbdata)
                ch_g2.append(ch_g2data)
                ch_fitted_g2.append(ch_fitted_g2)

                avg_indexdata = []
                avg_intensitydata = []
                avg_markdata = []
                avg_betadata = []
                avg_aDbdata = []
                avg_g2data = []
                avg_fitted_g2data = []
            
                for k in range(0,numAverages[0]):
                    avg_indexdata.append(np.fromfile(data,count=1,dtype='uint32'))
                    avg_intensitydata.append(np.fromfile(data,count=1,dtype='float64'))
                    avg_markdata.append(np.fromfile(data,count=1,dtype='uint32'))
                    avg_betadata.append(np.fromfile(data,count=1,dtype='float64'))
                    avg_aDbdata.append(np.fromfile(data,count=1,dtype='float64'))
                    avg_g2data.append(np.fromfile(data,count=int(numBins[0]),dtype='float64'))
                    avg_fitted_g2data.append(np.fromfile(data,count=int(numBins[0]),dtype='float64'))  
                
                
                avg_index.append(avg_indexdata)
                avg_intensity.append(avg_intensitydata)
                avg_mark.append(avg_markdata)
                avg_beta.append(avg_betadata)
                avg_aDb.append(avg_aDbdata)
                avg_g2.append(avg_g2data)
                avg_fitted_g2.append(avg_fitted_g2data)  
        
        
        output = list(zip(timings,ch_index, ch_intensity, ch_mark, ch_beta, ch_aDb, ch_g2, ch_fitted_g2, avg_index, avg_intensity, avg_mark, avg_beta, avg_aDb, avg_g2, avg_fitted_g2))
        
        output = np.array(output,dtype='object')
        g2 = avg_g2
        tau = np.stack(tau)
        #if g2.ndim > 1:
        #    g2 = g2.reshape(g2.shape[0],g2.shape[1])
        return [tau, g2]




