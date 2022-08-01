from sklearn.feature_selection import SelectKBest, chi2
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score,confusion_matrix, silhouette_score, accuracy_score, classification_report, mean_squared_error
from sklearn.decomposition import PCA
from sklearn.manifold import LocallyLinearEmbedding
import matplotlib
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pandas as pd
from scanfiles import scanning

def lle(x,y):
    acb,errb = randomfor(x,y)
    lle1 = LocallyLinearEmbedding(n_components=3)
    xx = lle1.fit_transform(x)
    aca, erra = randomfor(xx, y)
    return acb, errb, aca,erra

def chichi(df,xtrain,ytrain,xtest, ytest, supports):
        
    select_feature = SelectKBest(chi2, k=2).fit(xtrain, ytrain)
    selecttrain = select_feature.transform(xtrain)
    selecttest = select_feature.transform(xtest)
        
    clf_rf_2 = RandomForestClassifier(random_state=34).fit(selecttrain,ytrain)     
    ac_2 = accuracy_score(ytest,clf_rf_2.predict(selecttest))
    print('Accuracy of chi squared is: ',ac_2)
    rep = classification_report(ytest, clf_rf_2.predict(selecttest))
    print('classification report:', rep)
    cm_2 = confusion_matrix(ytest,clf_rf_2.predict(selecttest))
    sns.heatmap(cm_2,annot=True,fmt="d", cmap='viridis')
    plt.show()
    plt.close()
    indices = select_feature.get_support(indices=True)
    print(indices)
    signvars = []
    for i in range(len(indices)):
        signvars.append(supports[indices[i]])
    sns.jointplot(df.loc[:,signvars[0]], df.loc[:,signvars[1]], kind="reg", color="#ce1414")
    plt.tight_layout()

    plt.show()
    return ac_2
    
def mypca(x,y,n):
    print(x.shape)
    print("y:", y)
    acb,errb = randomfor(x,y)
    pca = PCA(n_components=n,svd_solver='full', random_state=101)
    xx = pca.fit_transform(x)
    df = pd.DataFrame(np.hstack((xx,y[:,np.newaxis])))
    print(df.head)
    aca,erra = randomfor(xx,y)
    fig, axs = plt.subplots(figsize=(14, 13))
    plt.clf()
    plt.plot(pca.explained_variance_ratio_, color='lime')
    plt.axis('tight')
    axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.ylabel("explained variance ratio")
    plt.xlabel("n")
        #plt.savefig(f"C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\testplots\\featselectvis\\pca_{specs[species]}.png", dpi=400)
       
    plt.show()
    return acb,errb,aca,erra
    
def randomfor(xs,ys):
    splitted = train_test_split(xs, ys, test_size=0.4, random_state=34)
    xtrain = splitted[0]; xtest = splitted[1]; ytrain = splitted[2]; ytest = splitted[3]
    xtrain = np.abs(xtrain); xtest = np.abs(xtest)
    #random forest classifier with n_estimators=10 (default)
    clf_rf = RandomForestClassifier(n_estimators=700, random_state=34).fit(xtrain,ytrain)
    predictionforest = clf_rf.predict(xtest)
    print(confusion_matrix(ytest,predictionforest))
    print(classification_report(ytest,predictionforest))
    ac = accuracy_score(ytest,predictionforest)
    err = mean_squared_error(ytest,predictionforest, squared='False')
    print('Accuracy of random forest is: ',ac)
    cm = confusion_matrix(ytest,predictionforest)
    sns.heatmap(cm,annot=True,fmt="d")
    #plt.savefig(f"C:\\Users\\user\\mres\\data_analysis\\data\\uclh_data\\randforest.png", dpi=400)

    plt.show()
    plt.close('all')
    return ac, err