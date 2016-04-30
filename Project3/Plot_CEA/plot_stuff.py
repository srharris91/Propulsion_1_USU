import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad
import sys

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

# get argument list using sys module
sys.argv
filename1= str(sys.argv[1])


#user defined values
# number of pressure values
PS = 8
W=2

x1=get_data_dtypes(filename1)

# calculate C*
c_star = np.zeros(np.size(x1['gam']))
j=0
for i in x1['gam']:
    g=x1['gam'][j]
    t=x1['t'][j]
    mw=x1['mw'][j]
    c_star[j]=(np.sqrt(g*8.314*1000.)/(g*np.sqrt(2./(g+1))**((g+1.)/(g-1.))) )*np.sqrt(t/mw)
    if (x1['pip'][j]==1.7428 and x1['cf'][j]>=0.6613 and x1['cf'][j]<=0.6614 and x1['of'][j]==8.5 and x1['aeat'][j]==1.0 and x1['mach'][j]==1.0 and x1['pip'][j]!=1.0): print 'c_star = ',c_star[j],j
    j=j+1



#print c_star




# plot big plot
fig=plt.figure(figsize=(24,72))
n=1
ax1=fig.add_subplot(np.size(x1.dtype.names)+1,W,1)
for i in x1.dtype.names:
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+1,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(x1['of'][j+k::W*PS],x1[i][j+k::W*PS],'.-')
        #if (i == 'isp'): ax1.set_xlabel(r'$O/F$')
        if (j == 0): ax1.set_ylabel(i)



for j in range(0,W):
    ax1=fig.add_subplot(np.size(x1.dtype.names)+1,W,n)
    n=n+1
    for k in range(0,W*PS-2,W):
        ax1.plot(x1['of'][j+k::W*PS],c_star[j+k::W*PS],'.-')
    ax1.set_xlabel(r'$O/F$')
    if (j == 0): ax1.set_ylabel('c*')





fig.savefig('all.png',bbox_inches='tight')
