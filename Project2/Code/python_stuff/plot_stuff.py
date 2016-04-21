import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad

def get_data(filename):
#    header=np.genfromtxt(filename,dtype=str)[0,:]
    data=np.genfromtxt(filename,skip_header=0)
    return data

def plot_data(xlabel,x,ylabel,y,filename):
    fig=plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(x,y,'o')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.axis('equal')
    fig.savefig(filename,bbox_inches='tight') 

x=get_data('answer.txt')

#fig=plt.figure(figsize=(16,28))
#
#ax1=fig.add_subplot(6,2,1)
#ax1.semilogx(x[:,0],x[:,1],'.')
#ax1.set_ylabel(r'$\dot{P_0}$')
#
#ax1=fig.add_subplot(6,2,2)
#ax1.semilogx(x[:,0],x[:,2],'.')
#ax1.set_ylabel(r'$\dot{r}$')
#
#ax1=fig.add_subplot(6,2,3)
#ax1.semilogx(x[:,0],x[:,3]/1000,'.')
#ax1.set_ylabel(r'$P_0 (kPa)$')
#
#ax1=fig.add_subplot(6,2,4)
#ax1.semilogx(x[:,0],x[:,4],'.')
#ax1.set_ylabel(r'$r$')
#
#
#ax1=fig.add_subplot(6,2,5)
#ax1.semilogx(x[:,0],x[:,5],'.')
#ax1.set_ylabel(r'$Burn Area (m^2)$')
#
#
#ax1=fig.add_subplot(6,2,6)
#ax1.semilogx(x[:,0],x[:,6]/1.,'.')
#ax1.set_ylabel(r'$Prop Mass (kg)$')
#
#ax1=fig.add_subplot(6,2,7)
#ax1.semilogx(x[:,0],x[:,7]/1.,'.',label='MassDepletion')
#ax1.set_ylabel(r'$(kg/s)$')
#ax1.semilogx(x[:,0],x[:,8],'.',label='ChokingMassFlow')
#ax1.legend(loc='best',numpoints=1)
#
#ax1=fig.add_subplot(6,2,9)
#ax1.semilogx(x[:,0],x[:,9]/1.,'.')
#ax1.set_ylabel(r'$Thrust (N)$')
#
#ax1=fig.add_subplot(6,2,10)
#ax1.semilogx(x[:,0],x[:,10],'.')
#ax1.set_ylabel(r'$I_{sp}(s)$')
#
#ax1=fig.add_subplot(6,2,11)
#ax1.semilogx(x[:,0],x[:,11],'.')
#ax1.set_xlabel(r'$t$')
#ax1.set_ylabel(r'$Port (Ma)$')
#
#ax1=fig.add_subplot(6,2,12)
#ax1.semilogx(x[:,0],x[:,12],'.')
#ax1.set_xlabel(r'$t$')
#ax1.set_ylabel(r'$Nozzle (Ma)$')
#
#fig.savefig('answer.png',bbox_inches='tight')





fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,1],'.')
ax1.set_ylabel(r'$\dot{P_0}$')
ax1.set_xlabel(r'$t$')
fig.savefig('P0_dot.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,2],'.')
ax1.set_ylabel(r'$\dot{r}$')
ax1.set_xlabel(r'$t$')
fig.savefig('r_dot.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,3]/1000,'.')
ax1.set_ylabel(r'$P_0 (kPa)$')
ax1.set_xlabel(r'$t$')
fig.savefig('P0.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,4],'.')
ax1.set_ylabel(r'$r$')
ax1.set_xlabel(r'$t$')
fig.savefig('r.pdf',bbox_inches='tight')


fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,5],'.')
ax1.set_ylabel(r'$Burn Area (m^2)$')
ax1.set_xlabel(r'$t$')
fig.savefig('BurnArea.pdf',bbox_inches='tight')


fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,6]/1.,'.')
ax1.set_ylabel(r'$Prop Mass (kg)$')
ax1.set_xlabel(r'$t$')
fig.savefig('PropMass.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,7]/1.,'.',label='MassDepletion')
ax1.set_ylabel(r'$(kg/s)$')
ax1.set_xlabel(r'$t$')
ax1.semilogx(x[:,0],x[:,8],'.',label='ChokingMassFlow')
ax1.legend(loc='best',numpoints=1)
fig.savefig('MassFlows.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,9]/1.,'.')
ax1.set_ylabel(r'$Thrust (N)$')
ax1.set_xlabel(r'$t$')
fig.savefig('Thrust.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,10],'.')
ax1.set_ylabel(r'$I_{sp}(s)$')
ax1.set_xlabel(r'$t$')
fig.savefig('Isp.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,11],'.')
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$Port (Ma)$')
fig.savefig('PortMach.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.semilogx(x[:,0],x[:,12],'.')
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$Nozzle (Ma)$')
fig.savefig('NozzMach.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.plot(x[:,2],x[:,3]/1000.,'.')
ax1.set_xlabel(r'$\dot{r}$')
ax1.set_ylabel(r'$P_0 (kPa)$')
fig.savefig('Reg_v_P0.pdf',bbox_inches='tight')
