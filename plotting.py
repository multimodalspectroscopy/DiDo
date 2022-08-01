import matplotlib.pyplot as plt
import numpy as np

def plotcomp(name,wavelengths,low,high,mydata,data):
    mymua = mydata[0]; mymus = mydata[1]
    mua = data[0]; mus = data[1]
    for m in (low,high):
        fig, axs = plt.subplots(2,1)
        axs[0].plot(wavelengths,mus,color='orange', label='solution')
        axs[0].fill_between(wavelengths,mus+np.std(mus), mus-np.std(mus), alpha=0.2, color='tomato')
        axs[0].plot(wavelengths, mymus, color='red', label='proposed')
        axs[0].legend()
        axs[0].set_ylabel('$\mu_s$ / ')
        axs[0].set_xlabel("wavelengths / nm")
        axs[1].plot(wavelengths,mymua,color='green', label='proposed')
        axs[1].plot(wavelengths,mua,color='lime', label='solution')
        axs[1].fill_between(wavelengths,mua+np.std(mua), mua-np.std(mua), alpha=0.2, color='lime')
        axs[1].legend()
        axs[1].set_ylabel('$\mu_a$ / ')
        axs[1].set_xlabel("wavelengths / nm")
        plt.tight_layout()   
        plt.savefig(f"code_development\\tuning\\muamus_{name}_{wavelengths[m]}", dpi=400)
        #plt.show()
        plt.close()
