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
OFS= 10

#x=get_data('MAE5540_Nitrous_HTBN.plt')
#x1=get_data_dtypes('MAE5540_Nitrous_HTBN.plt')
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
#mach=x1['mach']
#ivac=x1['ivac']
#isp =x1['isp']

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

# fig=plt.figure(figsize=(28,42))
# 
# ax1=fig.add_subplot(14,3,1)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
# ax1.set_ylabel(r'$P$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,2)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
# ax1.set_ylabel(r'$P$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,3)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
# ax1.set_ylabel(r'$P$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,4)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
# ax1.set_ylabel(r'$T$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,5)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
# ax1.set_ylabel(r'$T$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,6)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
# ax1.set_ylabel(r'$T$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,7)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
# ax1.set_ylabel(r'$\rho$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,8)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
# ax1.set_ylabel(r'$\rho$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,9)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
# ax1.set_ylabel(r'$\rho$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,10)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
# ax1.set_ylabel(r'$h$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,11)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
# ax1.set_ylabel(r'$h$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,12)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
# ax1.set_ylabel(r'$h$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,13)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
# ax1.set_ylabel(r'$g$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,14)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
# ax1.set_ylabel(r'$g$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,15)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
# ax1.set_ylabel(r'$g$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,16)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
# ax1.set_ylabel(r'$m$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,17)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
# ax1.set_ylabel(r'$m$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,18)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
# ax1.set_ylabel(r'$m$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,19)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
# ax1.set_ylabel(r'$mw$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,20)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
# ax1.set_ylabel(r'$mw$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,21)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
# ax1.set_ylabel(r'$mw$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,22)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
# ax1.set_ylabel(r'$cp$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,23)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
# ax1.set_ylabel(r'$cp$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,24)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
# ax1.set_ylabel(r'$cp$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,25)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
# ax1.set_ylabel(r'$gam$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,26)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
# ax1.set_ylabel(r'$gam$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,27)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
# ax1.set_ylabel(r'$gam$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,28)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
# ax1.set_ylabel(r'$mach$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,29)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
# ax1.set_ylabel(r'$mach$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,30)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
# ax1.set_ylabel(r'$mach$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,31)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
# ax1.set_ylabel(r'$ivac$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,32)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
# ax1.set_ylabel(r'$ivac$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,33)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
# ax1.set_ylabel(r'$ivac$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,34)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
# ax1.set_ylabel(r'$isp$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,35)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
# ax1.set_ylabel(r'$isp$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,36)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
# ax1.set_ylabel(r'$isp$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,37)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],pip[i::3*PS],'.-')
# ax1.set_ylabel(r'$pip$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,38)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],pip[i::3*PS],'.-')
# ax1.set_ylabel(r'$pip$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,39)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],pip[i::3*PS],'.-')
# ax1.set_ylabel(r'$pip$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# ax1=fig.add_subplot(14,3,40)
# for i in range(0,3*PS,3):
#     ax1.plot(of[i::3*PS],cf[i::3*PS],'.-')
# ax1.set_ylabel(r'$cf$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,41)
# for i in range(1,3*PS,3):
#     ax1.plot(of[i::3*PS],cf[i::3*PS],'.-')
# ax1.set_ylabel(r'$cf$')
# ax1.set_xlabel(r'$O/F$')
# 
# ax1=fig.add_subplot(14,3,42)
# for i in range(2,3*PS,3):
#     ax1.plot(of[i::3*PS],cf[i::3*PS],'.-')
# ax1.set_ylabel(r'$cf$')
# ax1.set_xlabel(r'$O/F$')
# 
# 
# fig.savefig('all.png',bbox_inches='tight')





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
