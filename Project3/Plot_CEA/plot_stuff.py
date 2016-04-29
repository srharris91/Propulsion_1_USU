import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad

def get_data(filename):
#    header=np.genfromtxt(filename,dtype=str)[0,:]
    data=np.genfromtxt(filename,skip_header=1)
    return data
def get_data_dtypes(filename):
    data=np.genfromtxt(filename,names=True,dtype=None)
    return data


def plot_data(xlabel,x,ylabel,y,filename):
    fig=plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(x,y,'o')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.axis('equal')
    fig.savefig(filename,bbox_inches='tight') 

#user defined values
# number of pressure values
PS = 8

#x1=get_data_dtypes('MAE5540_Nitrous_HTBN.plt')

x1=get_data_dtypes('Part1_try1.plt')
#of  =x1['of']
#p   =x1['p']
#t   =x1['t']
#rho =x1['rho']
#h   =x1['h']
#g   =x1['g']
#m   =x1['m']
#mw  =x1['mw']
#cp  =x1['cp']
#gam =x1['gam']
#pip =x1['pip']
#mach=x1['mach']
#cf  =x1['cf']
#ivac=x1['ivac']
#isp =x1['isp']

fig=plt.figure(figsize=(28,42))
n=1
ax1=fig.add_subplot(np.size(x1.dtype.names),3,1)
for i in x1.dtype.names:
    for j in range(1,4):
        ax1=fig.add_subplot(np.size(x1.dtype.names),3,n)
        n=n+1
        for k in range(2,3*PS,3):
            ax1.plot(x1['of'][k::3*PS],x1[i][k::3*PS],'.-')
        if (i == 'isp'): ax1.set_xlabel(r'$O/F$')
        if (j == 1): ax1.set_ylabel(i)
fig.savefig('all.png',bbox_inches='tight')
