import numpy as np
import matplotlib.pyplot as plt

def rotation(theta):
    R = np.eye(2)
    R[0,0] = np.cos(theta)
    R[1,1] = np.cos(theta)
    R[0,1] = np.sin(theta)
    R[1,0] = -np.sin(theta)
    return R


def myaxis(xlims,ylims):
    """
    xarray: array of x axis in graph
    """
  
    if (xlims[1]-xlims[0])>= (ylims[1]-ylims[0]):
        limm = np.max(np.abs(xlims))
    else:
        limm = np.max(np.abs(ylims))
    xarray = np.linspace(-limm, limm, 10)
    X = np.array(list(zip(xarray,np.zeros(xarray.shape))))
  
    R1 = rotation(np.pi/4)
    rpwr = np.dot(R1,X.T).T
    
    plt.plot(rpwr[:,0],rpwr[:,1], color='blue')
    scaling = (np.max(rpwr[:,0])+np.abs(np.min(rpwr[:,0])))/100
    plt.arrow(rpwr[0,0],rpwr[0,1], dx=-scaling, dy=scaling, color='blue', width=scaling/2, length_includes_head=False)
    plt.text(rpwr[0,0],rpwr[0,1], "rCST", color='blue')

    rcst = np.dot(rotation(np.pi/2),rpwr.T)
    disx = np.mean(rcst[0,:])-np.mean(rpwr[:,0])
    disy = np.mean(rcst[1,:])-np.mean(rpwr[:,1])

    plt.plot(rcst[0,:]-disx, rcst[1,:]-disy, color='deepskyblue')
    
    plt.arrow(rcst[0,0]-disx, rcst[1,0]-disy, dx=scaling, dy=scaling, color='deepskyblue', width=scaling/2, length_includes_head=False)
    plt.text(rcst[0,0]-disx, 0.9*rcst[1,0]-disy, "rPWR", color='deepskyblue')
   