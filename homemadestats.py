"""    def rungroup(self):
        interval = 500
        go = []; bo = []
        for n in range(len(self.allfiles)):
            if self.filechoice(n) == True and self.fetchoutcome(n) == 1:
                go.append(self.allfiles[n])
            if self.filechoice(n) == True and self.fetchoutcome(n) == 0:
                bo.append(self.allfiles[n])
        # repeat for every variable
        P = []; varsa = []
        for var in self.variables:
            datago = []; databo = []
            # gather group data
            for file in go:
                f = scanning("D:\\uclh_data")
                data,df = f.individual_data(file, interval, var)
                datago.append(list(data[0]))
            for file in bo:
                f = scanning("D:\\uclh_data")
                data,df = f.individual_data(file, interval, var)
                databo.append(list(data[0]))
            
            
            # statistical testing
            #datago = np.concatenate(datago,axis=1).flatten(); databo = np.concatenate(databo,axis=1).flatten()
            #s, p = stats.kruskal(datago, databo, nan_policy='omit')
            #P.append(p); varsa.append(var)
        #saving = list(zip(varsa, P))
        #name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\kruskal_groups_full.csv"
        #np.savetxt(name, saving, header=str(self.variables), fmt='%s', delimiter=';')
        print("calculations completed")


    def runcorr(self):
        #q = pwrcst("C:\\Users\\user\mres\\data_analysis\\data\\uclh_data")
        #q.calcuniversal(interval)
        interval = 500
        group = 'no'
        if group == 'yes':
            go = []; bo = []
            # retrieve patient outcome
            for n in range(len(self.allfiles)):
                if self.filechoice(n) == True and self.fetchoutcome(n) == 1:
                    go.append(self.allfiles[n])
                if self.filechoice(n) == True and self.fetchoutcome(n) == 0:
                    bo.append(self.allfiles[n])
            datago = []; databo = []
            R = np.eye(len(self.variables)); P = np.eye(len(self.variables)); M = np.eye(len(self.variables))
            for var in self.variables:

                datago = []; databo = []
                # gather group data
                for file in go:
                    f = scanning("D:\\uclh_data")
                    data,df = f.individual_data(file, interval, var)
                    datago.append(list(data[0]))
                for file in bo:
                    f = scanning("D:\\uclh_data")
                    data,df = f.individual_data(file, interval, var)
                    databo.append(list(data[0]))
                
                data1 =  np.concatenate(datago, axis=0)
                data2 = np.concatenate(databo, axis=0)
                
                data2 = np.pad(data2, int((len(data1)-len(data2))/2))
                data = [data1,data2]
                out = ['good','bad']
                df1 = pd.DataFrame(dict(list(zip(out,data))))
                print(df1)
                g = sns.pairplot(df1, corner=True, dropna='True', plot_kws=dict({"s":12}), diag_kws=dict({'color':'blue'}))
                plt.show()
                
                if len(data1) > 2:
                    r, p = stats.spearmanr(data1, data2, nan_policy='omit')
                    R[self.variables.index(var),self.variables.index(var)] = round(r,4)
                    P[self.variables.index(var),self.variables.index(var)] = round(p,4)

        #name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\spearmanr_byoutcome.csv"
        #np.savetxt(name, R, delimiter=';', fmt='%s', header=str(self.variables))
        #name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\spearmanp_byoutcome.csv"
        #np.savetxt(name, P, delimiter=';', fmt='%s', header=str(self.variables))

        if group == 'no':
            for n in range(len(self.allfiles)):  
                if self.filechoice(n) == True:
                    f = scanning("D:\\uclh_data")
                    print(self.variables)
                    data,df = f.individual_data(self.allfiles[n], interval, *self.variables)
                    #self.visualise(self.allfiles[n],df)
                    self.boxplotting(df, n)
                    
        print("calculations completed")
      
        return self.allfiles
    

    def runts(self):
        filearr = []
        # find pairs of datasets
        for n in range(len(self.allfiles)):
            if self.filechoice(n) == True:
                filearr.append(self.allfiles[n])
        combs = list(set(combinations(filearr,2)))
        combs = self.checkcombs(combs)
        TP = []
        # for each variable, calculate file by file  
        for k in range(len(self.variables)):
            TP.append([np.eye(len(filearr))])
        TP = np.array(TP, dtype='object').reshape(len(self.variables), len(filearr), len(filearr))
        TT = TP.copy()
        
        for m in range(len(combs)):  
            f = scanning("D:\\uclh_data")

            data1, df1 = f.individual_data(combs[m][0], 500,*self.variables)
            data2, df2 = f.individual_data(combs[m][1], 500,*self.variables)
           
            for vv in range(len(self.variables)):

                data1v = df1[self.variables[vv]]
                data2v = df2[self.variables[vv]]
                
                #t = stats.ks_2samp(data1v, data2v)
                #TP[vv, filearr.index(str(combs[m][1])),filearr.index(str(combs[m][0]))] = round(t[1],4)
                #TT[vv, filearr.index(str(combs[m][1])),filearr.index(str(combs[m][0]))] = round(t[0],4)
              
        for k in range(TP.shape[0]):
            name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\"+ str(self.variables[k]) +"_ttest_p_vals.csv"
            np.savetxt(name, TP[k,:,:], delimiter=';', fmt='%s', header = str(filearr))
            name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\"+ str(self.variables[k]) +"_ttest_t_vals.csv"
            np.savetxt(name, TT[k,:,:], delimiter=';', fmt='%s', header = str(filearr))

        return name

    def rs(self,group,*args): 
        # calculate pearson r as before, but data is group data 
            
        if group == 'no':    
            file,df,data,interval = args
            tops = list(set(combinations(self.variables,2)))
            filetops = []
            R = np.eye(len(self.variables)); P = np.eye(len(self.variables)); M = np.eye(len(self.variables))
            for m in range(len(tops)):
                data1 = np.array(df[tops[m][0]]); data2 = np.array(df[tops[m][1]])
                idx = np.where(~np.isnan(data1+data2))
                data1 = data1[idx]; data2 = data2[idx]
                #print(data1, data2)
                if len(data1) > 2:
                    r, p = stats.pearsonr(data1, data2)
                    error = np.cov(m=data1, y=data2)
                    R[self.variables.index(str(tops[m][1])),self.variables.index(str(tops[m][0]))] = round(r,4)
                    P[self.variables.index(str(tops[m][1])),self.variables.index(str(tops[m][0]))] = round(p,4)
            print(R)
            #name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\"+ str(file[0:25]) + "_" + str(interval)+"_r_vals.csv"
            #np.savetxt(name, R, delimiter=';', fmt='%s', header=str(self.variables))
            #name = "C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\corrs\\averagingtest\\"+ str(file[0:25]) + "_" + str(interval) +"_p_vals.csv"
            #np.savetxt(name, P, delimiter=';', fmt='%s', header=str(self.variables))

        return R,P"""