import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad

def get_data(filename):
#    header=np.genfromtxt(filename,dtype=str)[0,:]
    data=np.genfromtxt(filename,skip_header=1)
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

x=get_data('MAE5540_Nitrous_HTBN.plt')
of  =x[:, 0]
p   =x[:, 1]
t   =x[:, 2]
rho =x[:, 3]
h   =x[:, 4]
g   =x[:, 5]
m   =x[:, 6]
mw  =x[:, 7]
cp  =x[:, 8]
gam =x[:, 9]
mach=x[:,10]
ivac=x[:,11]
isp =x[:,12]


fig=plt.figure(figsize=(20,42))

ax1=fig.add_subplot(12,3,1)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
ax1.set_ylabel(r'$P$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,2)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
ax1.set_ylabel(r'$P$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,3)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],p[i::3*PS],'.-')
ax1.set_ylabel(r'$P$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,4)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
ax1.set_ylabel(r'$T$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,5)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
ax1.set_ylabel(r'$T$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,6)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],t[i::3*PS],'.-')
ax1.set_ylabel(r'$T$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,7)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
ax1.set_ylabel(r'$\rho$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,8)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
ax1.set_ylabel(r'$\rho$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,9)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],rho[i::3*PS],'.-')
ax1.set_ylabel(r'$\rho$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,10)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
ax1.set_ylabel(r'$h$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,11)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
ax1.set_ylabel(r'$h$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,12)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],h[i::3*PS],'.-')
ax1.set_ylabel(r'$h$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,13)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
ax1.set_ylabel(r'$g$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,14)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
ax1.set_ylabel(r'$g$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,15)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],g[i::3*PS],'.-')
ax1.set_ylabel(r'$g$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,16)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
ax1.set_ylabel(r'$m$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,17)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
ax1.set_ylabel(r'$m$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,18)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],m[i::3*PS],'.-')
ax1.set_ylabel(r'$m$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,19)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
ax1.set_ylabel(r'$mw$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,20)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
ax1.set_ylabel(r'$mw$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,21)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],mw[i::3*PS],'.-')
ax1.set_ylabel(r'$mw$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,22)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
ax1.set_ylabel(r'$cp$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,23)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
ax1.set_ylabel(r'$cp$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,24)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],cp[i::3*PS],'.-')
ax1.set_ylabel(r'$cp$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,25)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
ax1.set_ylabel(r'$gam$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,26)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
ax1.set_ylabel(r'$gam$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,27)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],gam[i::3*PS],'.-')
ax1.set_ylabel(r'$gam$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,28)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
ax1.set_ylabel(r'$mach$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,29)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
ax1.set_ylabel(r'$mach$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,30)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],mach[i::3*PS],'.-')
ax1.set_ylabel(r'$mach$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,31)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
ax1.set_ylabel(r'$ivac$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,32)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
ax1.set_ylabel(r'$ivac$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,33)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],ivac[i::3*PS],'.-')
ax1.set_ylabel(r'$ivac$')
ax1.set_xlabel(r'$O/F$')


ax1=fig.add_subplot(12,3,34)
for i in range(0,3*PS,3):
    ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
ax1.set_ylabel(r'$isp$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,35)
for i in range(1,3*PS,3):
    ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
ax1.set_ylabel(r'$isp$')
ax1.set_xlabel(r'$O/F$')

ax1=fig.add_subplot(12,3,36)
for i in range(2,3*PS,3):
    ax1.plot(of[i::3*PS],isp[i::3*PS],'.-')
ax1.set_ylabel(r'$isp$')
ax1.set_xlabel(r'$O/F$')


fig.savefig('all.png',bbox_inches='tight')
