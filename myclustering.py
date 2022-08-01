import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.pylab as pl
from sklearn.cluster import KMeans


def clustering(x, y, patients, labs, axs):

    patients = np.array(patients, dtype=int)
    
    plotting = 'on'
    if plotting == 'on':
        markers = [ 'o', '<','o','o']
        colors = ['k','k','orange', 'red']
        facecols = [ 'black','black',  'white', 'white']  
        columns = [plt.plot([], [], markers[m], markerfacecolor=facecols[m], markeredgecolor=colors[m])[0] for m in range(2)]
        labels = [ 'cluster 0' , 'cluster 1', 'LACT/NAA class 0', 'LACT/NAA class 1']

        overlapcheck = 'off'
        if overlapcheck == 'on':
            markers = [ 'o','o', '--', '--', '--']
            colors = ['orange', 'red', 'gray', 'tomato', 'green']
            facecols = [ 'white', 'white','gray', 'tomato', 'greem']                  
            columns = [plt.plot([], [], markers[m], markerfacecolor=facecols[m], markeredgecolor=colors[m])[0] for m in range(0,5)]
            labels = [ 'mild HIE', 'moderate HIE', 'sess1 decision', 'sess2 decision', 'sess3 decision']

        axs.legend(handles=columns, labels=labels, loc='best')
        markers = ['o', '<', 's']
        centrings = ['center', 'left', 'right']
        labs = labs[~np.isnan(labs)]
        colors = pl.cm.jet(np.linspace(0,1,int(len(labs))))
        for m in range(len(labs)):
            if np.isnan(labs[m])==False:                    
                if labs[m] == 0 :
                    marker = markers[0]
                if labs[m] == 1:
                    marker = markers[1]
                if labs[m] == 2:
                    marker = markers[2]
                   
                axs.scatter(x[m], y[m], marker = marker, color = colors[m])
                plt.text(x=x[m], y=y[m], s=patients[m], ha=centrings[np.random.randint(0,2)])

        axs.set_ylabel("PCA 1")
        axs.set_xlabel("PCA 2")
        axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        plt.ylim(-0.1,1.1)
        plt.xlim(-0.1,1.1)
        plt.tight_layout()
    return labs, columns, labels

def pigletclustering(results,i,vv,idx,x,name1,name0,diagnosis,session, N, usingcontrols,axs):
    """
    results: empty array
    df: data
    i: number of type
    vv: number of comb
    idx: only for ref plotting not experiment
    """
    

    scoreagr = 0
    scoredis  = 0
    #featnames = ['AftInsAvg' , 'AUCInsult', 'BasAvg', "InsAvg","PeakLoc","PeakVal","RecAvg","RF1","RF2"]
    #featnames = ['minima', 'medians', 'maxima', 'auc']
    
    mycl = KMeans(n_clusters=N, algorithm='full', random_state=34)
    labs = np.array(mycl.fit(x).labels_)
  
    for k in range(len(labs)):
        results[i,vv,k] = labs[k]
    plotting = 'on'
    if plotting == 'on':
        
        if N == 3 and usingcontrols == 'no':
            raise Warning('double check the number of clusters requested')
        
        if usingcontrols == 'yes':
            markers = ['s', 'o', '<','o','o','o']
            colors = ['k','k','k','blue','orange', 'red']
            facecols = ['black', 'black','black', 'white', 'white', 'white']
                            
            columns = [plt.plot([], [], markers[m], markerfacecolor=facecols[m], markeredgecolor=colors[m])[0] for m in range(6)]
            labels = ['cluster 0', 'cluster 1', 'cluster 2','controls', 'mild HIE', 'severe HIE']
            
        if usingcontrols == 'no':
            markers = [ 'o', '<','o','o']
            colors = ['k','k','orange', 'red']
            facecols = [ 'black','black',  'white', 'white']
                            
            columns = [plt.plot([], [], markers[m], markerfacecolor=facecols[m], markeredgecolor=colors[m])[0] for m in range(4)]
            labels = ['cluster 0', 'cluster 1', 'mild HIE', 'severe HIE']
            #axs.legend(handles=columns, labels=labels, loc=[1.05, 0])
        colors = pl.cm.jet(np.linspace(0,1,len(labs)))
        x = np.array(x)
        for m in range(len(labs)):
            if labs[m] != np.nan:                    
                if labs[m] == 0 :
                    marker = markers[0]
                if labs[m] == 1:
                    marker = markers[1]
                if labs[m] == 2:
                    marker = markers[2]
                if diagnosis[m] == 0:
                    edgecolor = 'blue'
                if diagnosis[m] == 1:
                    edgecolor = 'orange'
                if diagnosis[m] == 2:
                    edgecolor = 'red'

                if (labs[m] == 0 and diagnosis[m] == 0) or (labs[m] == 1 and diagnosis[m] == 1):
                    scoredis += 1
                if (labs[m] == 0 and diagnosis[m] == 1) or (labs[m] == 1 and diagnosis[m] == 0):
                    scoreagr += 1
                                    
                axs.scatter(x[:,1][m], x[:,0][m], marker=marker, color=colors[m], edgecolor=edgecolor)
                    
        
        #plt.savefig(f"C:\\Users\\user\\mres\\data_analysis\\data\\pigletsforml\\comp_{featnames[i]}_z({name0})_z({name1}).png")
        #plt.savefig(f"C:\\Users\\rmapdbi\\mres\\data_analysis\\data\\pigletsforml\\{usingcontrols}_{N}_{featnames[i]}_z({name0})_vs_z({name1}).png")
        #plt.show()
        #plt.close()

        
    return results, scoreagr, scoredis, labs, columns, labels