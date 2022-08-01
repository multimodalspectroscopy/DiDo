from turtle import pos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from numpy import round
from scipy.stats import pearsonr, ttest_ind
import copy

def findcorrelation(*args):
        dcin = []
        dataset1, dataset2, events, steps = args[0], args[1], args[2], args[3]            
        if len(events) > 0:
            for k,hyptime in enumerate(events):
                datain1 = np.array(dataset1[int(hyptime):int(hyptime+steps[k])]).flatten()
                datain2 = np.array(dataset2[int(hyptime):int(hyptime+steps[k])]).flatten()
                if len(datain1[~np.isnan(datain1)]) > 0 and len(datain2[~np.isnan(datain2)]) > 0:
                    indices = np.array(np.logical_not(np.logical_or(np.isnan(datain1), np.isnan(datain2))))
                    dc = pearsonr(datain1[indices], datain2[indices])[0]
                    #dc = ttest_ind(datain1, datain2, nan_policy='omit', equal_var=False)[0]
                else:  
                    dc = np.nan
                dcin.append(dc)
            dcin = np.array(dcin)
            dcout = copy.deepcopy(dcin)
        else:
            dcout = np.nan
        return dcout



def findchange(*args):
        dq = []; err = []
        data, events, steps = args[0], args[1], args[2]
        if len(args) == 4:
            turningpoints = np.array(args[3], dtype=int)
        if len(events) > 0:
            positions = []
            for k,hyptime in enumerate(events):
                datain = data[int(hyptime):int(hyptime+steps[k])]
                if len(datain[~np.isnan(datain)]) > 1:
                    if type(hyptime) is not np.int32:
                        if len(hyptime) == 1:
                            hyptime = hyptime[0]
                    calc_min = np.nanmean(data[hyptime-50:hyptime])-np.nanmin(datain)
                    if len(args) == 3:
                        calc = copy.deepcopy(calc_min)
                    if len(args) == 4:
                        calc_max = np.abs(np.nanmean(data[hyptime-50:hyptime])-np.nanmax(datain))
                        if np.abs(calc_min) >= calc_max:
                            calc = copy.deepcopy(calc_min)
                        else:
                            calc = copy.deepcopy(calc_max)
                    dq.append(calc)
                    points = np.argwhere(np.array(datain).flatten()-np.nanmin(datain)<1e-3).flatten()
                    if len(points) > 1:
                        points = points[0]
                    point = int(points)
                    positions.append(point)
                    if len(args) == 3:
                        err.append(np.nanmean(datain[point-2:point+2]) - np.nanmin(datain))
                    if len(args) == 4:
                        err.append(np.nanvar(datain[turningpoints[k]-2:turningpoints[k]+2]))
                else:             
                    err.append(np.nan)
                    dq.append(np.nan)
         
            positions = np.array(positions) 
            dqout = np.array(dq) 
            errout = np.array(err)    
        else: 
             
            errout = np.nan
            dqout = np.nan
            positions = np.nan
        if len(args) == 3:
            return dqout,errout, positions
        if len(args) == 4:
            return dqout, errout


def findhypoxic(filename, data):
        data = np.array(data)
        """fig, axs = plt.subplots()
        colors = plt.cm.jet(np.linspace(0,1,10))
        plt.plot(data, color =colors[np.random.randint(0,len(colors))])
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.ylabel('$SpO_2 - O_2$ / %')
        plt.xlabel(' time / s')"""
        #plt.show()
        ref = data[~np.isnan(data)][0]
        step = 1; event = 0; events = []; steps = []; pulses = []; mins = []; refs = []
        for i in range(len(data)):
            #print(i,event)
            if i >= 50:
                ref = np.nanmean(data[i-50:i])
                #print(ref, val)
            if np.abs(data[i] - ref) <= 0.5:
               startp = i
              
            # loop until you satisfy these conditions
            if (data[i] < 90  or data[i] > 1.5 * ref) and i > event + step:
                #print("oi")
                step = 1; convstep = False
                while convstep == False:
                    if data[int(i+step)] <  ref:
                        if step == 1:
                            event = i
                            pulse = data[event]
                        step += 1    
                    else:
                        convstep = True
                        break
  
                if convstep == True and step >= 30:
                    startp = event - startp
                    endp = 0
                    #plt.fill_between(x=np.arange(event - startp,event+step+endp),y1=np.linspace(100,100,step+int(endp+startp)), y2=np.min(data[event-startp:event+step+endp]), color='tomato', alpha=0.2)
                    #print("event at %s with duration %s" % (event, step))
                    events.append(event); steps.append(step); pulses.append(ref - pulse); mins.append(np.min(data[event-startp:event+step+endp]))
                    refs.append(ref)
        if len(events) >= 1:
            cent = np.mean(events)
            diffs1 = np.diff(events,1)
            if all(item < 200 for item in diffs1) == True:
                if cent + 500 < len(data) and cent - 500 > 0:
                    plt.xlim(cent-500,cent+500)
               # plt.savefig(f"C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\variableplots\\hypoxicevents\\method_0\\SpO2_{filename}_event_0")
               # plt.show()
            else:
                #print("here")
                for eventid in events:
                    cent = eventid
                    if cent + 500 < len(data) and cent - 500 > 0:
                        plt.xlim(cent-500,cent+500)
                    #plt.show()
                    #plt.savefig(f"C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\variableplots\\hypoxicevents\\method_0\\SpO2_{filename}_event_{events.index(eventid)}")
        plt.close('all')
        return pulses, round(events,3), round(steps,3),round(mins,3), round(refs,3)

def extract_corr_features(timing, changes, patnum, session, events, steps, dataspo2, dataa, datab):
    if len(events) > 0 and len(dataa) >= 3600 :
        if events[-1] + 200 < len(dataa):
            #eventuse, stepsuse = checkoverlap()
            if timing == 'before':
                eventuse = np.array(events-200).flatten()
                stepsuse = np.linspace(200,200,len(steps)).flatten()
            if timing == 'during':
                eventuse = copy.deepcopy(events)
                stepsuse = copy.deepcopy(steps)
            if timing == 'after':
                eventuse = np.array(events+steps).flatten()
                stepsuse = np.linspace(200,200,len(steps)).flatten()
            eventsin = copy.deepcopy(eventuse)
            stepsin = copy.deepcopy(stepsuse)
                        
            # calculating errors
            spd, errspd, turningpoints = np.array(findchange(dataspo2, eventsin,stepsin))
            hbd,errhbd = np.array(findchange(dataa,eventsin,stepsin, turningpoints))
            oxcd,erroxcd = np.array(findchange(datab,eventsin,stepsin, turningpoints))
            # calculating correlations
            SHcorr = findcorrelation(dataspo2, dataa, eventsin, stepsin)
            HOcorr = findcorrelation(dataa, datab, eventsin, stepsin)
            SOcorr = findcorrelation(dataspo2, datab, eventsin, stepsin)

            dataout = np.vstack([spd,hbd,oxcd,                                        
                                        errspd, errhbd, erroxcd,
                                        SHcorr, HOcorr, SOcorr])
                           
            if len(eventuse) > 1:                                
                xx = np.argwhere(np.isnan(dataout))
                currentweights = [weight/np.nanmax(stepsuse) for weight in stepsuse]
                if len(np.unique(xx[:,1])) < dataout.shape[1]:
                    for nanidx in np.unique(xx[:,1]):
                        dataout[np.isnan(dataout)] = 0
                        currentweights[nanidx] = 0
                datatoprocess = np.average(np.array(dataout),axis=1,weights=np.array(currentweights).flatten()) 
            else:
                datatoprocess = copy.deepcopy(dataout.flatten())
            changes.append([
                          np.hstack([patnum,session,datatoprocess])
                                        ])

    if len(dataa) < 3600:
        changes.append([
                        patnum,session,
                            np.nan, np.nan, np.nan,                        
                            np.nan, np.nan, np.nan,
                            np.nan, np.nan, np.nan
                            ])
    return changes