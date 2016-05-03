import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad
import sys
import glob
import os

def get_data(filename):
#    header=np.genfromtxt(filename,dtype=str)[0,:]
    data=np.genfromtxt(filename,skip_header=1)
    return data
def get_data_dtypes(filename):
    data=np.genfromtxt(filename,names=True,dtype=None)
    return data

def calcs(thrust): 
    # calculate C*,P0,t0,at,mdot,isp_opt
    for i in range(0,np.size(x1['gam'])):
        g=x1['gam'][i]
        x1['p'][i]=x1['p'][i]*100000. # convert bar to Pa
        x1['isp'][i]=x1['isp'][i]/9.806 # convert European Isp to American Isp
        x1['ivac'][i]=x1['ivac'][i]/9.806 # convert European Isp to American Isp
        p=x1['p'][i]
        t=x1['t'][i]
        M=x1['mach'][i]
        mw=x1['mw'][i]
        cf=x1['cf'][i]

        t0[i]=t*(1+((g-1.)/2.)*M**2)
        p0[i]=p*(1+((g-1.)/2.)*M**2)**(g/(g-1.))
        #print "P = ",p,p0[i]
        c_star[i]=(np.sqrt(g*8.3144598*1000.)/(g*np.sqrt(2./(g+1))**((g+1.)/(g-1.))) )*np.sqrt(t0[i]/mw)
        if cf!=0: at[i]=thrust/(p0[i]*cf)
        p0_at[i]=p0[i]*at[i]
        if cf!=0: mdot[i]=thrust/(c_star[i]*cf)

        if (x1['pip'][i]==1.7428 and x1['cf'][i]>=0.6613 and x1['cf'][i]<=0.6614 and x1['of'][i]==8.5 and x1['aeat'][i]==1.0 and x1['mach'][i]==1.0 and x1['pip'][i]!=1.0): print 'c_star = ',c_star[i],i

        if (cf!=0): isp_opt[i]=(
                p0_at[i]/(9.806*mdot[i])*(
                    g*np.sqrt(
                        (2./(g-1))*(2./(g+1))**((g+1)/(g-1))
                        )*np.sqrt(
                            1-(p/p0[i])**((g-1)/g))
                    +
                    x1['aeat'][i]*(p-101325.)/p0[i]
                    ))

def big():# plot big plot
    fig=plt.figure(figsize=(24,72))
    n=1
    ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,1)
    for i in x1.dtype.names:
        for j in range(0,W):
            ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
            n=n+1
            for k in range(0,W*PS-(W-1),W):
                ax1.plot(100-x1['f'][j+k::W*PS],x1[i][j+k::W*PS],'.-')
            #if (i == 'isp'): ax1.set_xlabel(r'$\%\ H_2O_2$')
            if (j == 0): ax1.set_ylabel(i)

    # Isp optimal
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],isp_opt[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('Optimal $I_{sp}$')

    # c*
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],c_star[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('c*')

    # t0
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],t0[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('t0')

    # p0
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],p0[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('p0')
        
    # a_throat
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],at[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('at')

    # P0*a_throat
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],p0_at[j+k::W*PS],'.-')
        #ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel('p0 at')

    # m_dot
    for j in range(0,W):
        ax1=fig.add_subplot(np.size(x1.dtype.names)+7,W,n)
        n=n+1
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x1['f'][j+k::W*PS],mdot[j+k::W*PS],'.-')
        ax1.set_xlabel(r'$\%\ H_2O_2$')
        if (j == 0): ax1.set_ylabel(r'$\dot{m}$')

    fig.savefig('all.pdf',bbox_inches='tight')

def Part1_plots(): # Part 1 specific plots desired i
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['t'][W-1+k::W*PS],'.-',label='$T$ Nozzle Exit' if k==0 else "")
        ax1.plot(100-x1['f'][W-1+k::W*PS],t0[W-1+k::W*PS],'s:',label='$T_0$ Stagnation Nozzle Exit' if k==0 else "")
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$T\ (K)$')
    ax1.legend(numpoints=1,loc='upper center',bbox_to_anchor=(0.5,1.3))
    fig.savefig('Part1_i.pdf',bbox_inches='tight')

    # Part 1 specific plots desired ii
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-2+k::W*PS],t0[W-2+k::W*PS],'.-',label='$T_0$ Stagnation Combuster' if k==0 else "")
        ax1.plot(100-x1['f'][W-1+k::W*PS],t0[W-1+k::W*PS],'s:',label='$T_0$ Stagnation Nozzle Exit' if k==0 else "")
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$T\ (K)$')
    ax1.legend(numpoints=1,loc='upper center',bbox_to_anchor=(0.5,1.3))
    fig.savefig('Part1_ii.pdf',bbox_inches='tight')

    # Part 1 specific plots desired iii
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['mach'][W-1+k::W*PS],'.-',label='$Ma$ Nozzle Exit' if k==0 else "")
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$Ma$')
    #ax1.legend(numpoints=1,loc='upper center',bbox_to_anchor=(0.5,1.2))
    fig.savefig('Part1_iii.pdf',bbox_inches='tight')

    # Part 1 specific plots desired iv
    fig=plt.figure(figsize=(3,8))
    ax1=fig.add_subplot(211)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],c_star[W-1+k::W*PS],'.-',label='$C*$ Nozzle Exit' if k==0 else "")
    ax1.set_ylabel(r'$C*$')
    #ax1.legend(numpoints=1,loc='upper center')

    ax1=fig.add_subplot(212)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['cf'][W-1+k::W*PS],'.-',label='$C_F$ Nozzle Exit' if k==0 else "")
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$C_F$')
    #ax1.legend(numpoints=1,loc='upper center')
    fig.savefig('Part1_iv.pdf',bbox_inches='tight')

    # Part 1 specific plots desired v
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['ivac'][W-1+k::W*PS],'.-',label='$I_{vac}$ Nozzle Exit' if k==0 else "")
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['isp'][W-1+k::W*PS],'s:',label='$I_{sp}$ Nozzle Exit' if k==0 else "")
        ax1.plot(100-x1['f'][W-1+k::W*PS],isp_opt[W-1+k::W*PS],'^:',label='Optimal $I_{sp}$ Nozzle Exit' if k==0 else "")
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$I_{sp}\ (s)$')
    ax1.legend(numpoints=1,loc='upper center',bbox_to_anchor=(0.5,1.425))
    fig.savefig('Part1_v.pdf',bbox_inches='tight')

    # Part 1 specific peroxide concentration required?
    print "Peroxide concentration requires a 5:1 Oxide to fuel ratio to vaporize all water"

    # Part 1 What at is required for the optimal operating chamber pressure at 90% H2O2 concentration?, corresponding Isp 
    #print 0,W*PS-(W-1),W
    #print np.size(at[2+0::W*PS]),W*PS
    #print "A* and optimal Isp at 90% H2O2 for 3114N thrust"
    for k in range(0,W*PS-(W-1),W):
        for i in range(W-1+k,np.size(at),W*PS):
            #print i,k,x1['f'][i],at[i]
            if (x1['f'][i]==10 and x1['p'][i]==100440.): 
                print "90% H2O2 and optimal pressure gives:"
                print "  at =       ",at[i]
                print "  isp=       ",x1['isp'][i]
                print "  isp_opt=   ",isp_opt[i]
                print "  mdot=      ",mdot[i]

    # Part 1 specific plots desired vi
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],mdot[W-1+k::W*PS],'.-')
    ax1.set_xlabel(r'$\%\ H_2O_2$')
    ax1.set_ylabel(r'$\dot{m}\ \frac{kg}{s}$')
    fig.savefig('Part1_vi.pdf',bbox_inches='tight')

def output_file():# output file
    f=open('output.txt','w')
    # headings
    f.write('# ')
    for t in x1.dtype.names:
        f.write(t+' ')
    f.write('isp_opt ')
    f.write('t0 ')
    f.write('p0 ')
    f.write('c_star ')
    f.write('at ')
    f.write('p0at ')
    f.write('m_dot \n')
    #for i in x1.dtype.names:
        #for j in np.size(x1[i]):
    for i in range(0,np.size(x1['gam'])):
        for t in x1.dtype.names:
            f.write(str(x1[t][i]))
            f.write(' ')
        f.write(str(isp_opt[i]))
        f.write(' ')
        f.write(str(t0[i]))
        f.write(' ')
        f.write(str(p0[i]))
        f.write(' ')
        f.write(str(c_star[i]))
        f.write(' ')
        f.write(str(at[i]))
        f.write(' ')
        f.write(str(p0_at[i]))
        f.write(' ')
        f.write(str(mdot[i]))
        f.write(' \n')
    f.close()

def plot_data(xlabel,x,ylabel,y,filename):
    fig=plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(x,y,'o')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.axis('equal')
    fig.savefig(filename,bbox_inches='tight') 

def example(): # example from section 7.1
    global PS
    global W
    PS=5
    W=2
    # get argument list using sys module
    global x1
    x1=get_data_dtypes(filename1)

    # subroutines for part 1
    # initialize values
    global c_star
    global p0    
    global p0_at
    global t0  
    global at 
    global mdot
    global isp_opt

    c_star  = np.zeros(np.size(x1['gam']))
    p0      = np.zeros(np.size(x1['gam']))
    p0_at   = np.zeros(np.size(x1['gam']))
    t0      = np.zeros(np.size(x1['gam']))
    at      = np.zeros(np.size(x1['gam']))
    mdot    = np.zeros(np.size(x1['gam']))
    isp_opt = np.zeros(np.size(x1['gam']))

    calcs(110.)
    big()
    #Part1_plots()
    #output_file()


def Part1(): #  Part 1
    global PS
    global W
    PS=8
    W=3
    # get argument list using sys module
    global x1
    x1=get_data_dtypes(filename1)

    # subroutines for part 1
    # initialize values
    global c_star
    global p0    
    global p0_at
    global t0  
    global at 
    global mdot
    global isp_opt

    c_star  = np.zeros(np.size(x1['gam']))
    p0      = np.zeros(np.size(x1['gam']))
    p0_at   = np.zeros(np.size(x1['gam']))
    t0      = np.zeros(np.size(x1['gam']))
    at      = np.zeros(np.size(x1['gam']))
    mdot    = np.zeros(np.size(x1['gam']))
    isp_opt = np.zeros(np.size(x1['gam']))

    calcs(3114.)
    big()
    Part1_plots()
    output_file()

def Part2_iter(): #  Part 2
    global PS
    global W
    PS=8
    W=3
    # get argument list using sys module
    global x1
    x1=get_data_dtypes(filename1)

    # subroutines for part 1
    # initialize values
    global c_star
    global p0    
    global p0_at
    global t0  
    global at 
    global mdot
    global isp_opt

    c_star  = np.zeros(np.size(x1['gam']))
    p0      = np.zeros(np.size(x1['gam']))
    p0_at   = np.zeros(np.size(x1['gam']))
    t0      = np.zeros(np.size(x1['gam']))
    at      = np.zeros(np.size(x1['gam']))
    mdot    = np.zeros(np.size(x1['gam']))
    isp_opt = np.zeros(np.size(x1['gam']))

    calcs(3114.)
    # Part 2 specific plots for iteratin
    fig=plt.figure(figsize=(3,4))
    ax1=fig.add_subplot(111)
    for k in range(0,W*PS-(W-1),W):
        ax1.plot(100-x1['f'][W-1+k::W*PS],x1['p'][W-1+k::W*PS],'.-',label=k)
    ax1.set_xlabel(r'$\%\ ABS$')
    ax1.set_ylabel(r'$P\ (Pa)$')
    ax1.legend(numpoints=1,loc='best')
    fig.savefig('Part2_p.pdf',bbox_inches='tight')
    output_file()

def Part2(): #  Part 2, combined files
    # get all txt files
    W=3
    PS=8
    s= np.sort(glob.glob("./Part2/*/*.txt"))
    max_values=np.zeros(np.size(s))
    loc_of=np.zeros(np.size(s))
    loc_isp_opt=np.zeros(np.size(s))
    #loc_isp=np.zeros(np.size(s))
    fig = plt.figure(figsize=(3,3))
    fig2 = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax2=fig2.add_subplot(111)
    for f in range(0,np.size(s)):
        filename=s[f]
        x_data=get_data_dtypes(filename)
        max_values[f] = np.max(x_data['c_star'][2::3])
        for i in range(0,np.size(x_data['c_star'])):
            if (x_data['c_star'][i] == max_values[f] and x_data['isp'][i] !=0):
                loc_of[f]=x_data['of'][i]
                loc_isp_opt[f]=x_data['isp_opt'][i]
        for k in range(0,W*PS-(W-1),W):
            ax1.plot(100-x_data['f'][W-1+k::W*PS],x_data['isp_opt'][W-1+k::W*PS],'.-',label=filename if k==0 else '')
            ax2.plot(100-x_data['f'][W-1+k::W*PS],x_data['c_star'][W-1+k::W*PS],'.-',label=filename if k==0 else '')
    ax1.set_xlabel('O/F')
    ax2.set_xlabel('O/F')
    ax1.set_ylabel('$I_{sp}$')
    ax2.set_ylabel('$C*$')
    ax1.legend(numpoints=1,loc='best')
    ax2.legend(numpoints=1,loc='best')
    fig.savefig('Part2_Isps.pdf',bbox_inches='tight')
    fig2.savefig('Part2_C_stars.pdf',bbox_inches='tight')
    # c*
    fig = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(np.array([80,85,90,95,99]),max_values,'o-')
    ax1.set_xlabel('Peroxide Mass fraction')
    ax1.set_ylabel('C*')
    fig.savefig('Part2_C_star.pdf',bbox_inches='tight')

    # O/F
    fig = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(np.array([80,85,90,95,99]),loc_of,'o-')
    ax1.set_xlabel('Peroxide Mass fraction')
    ax1.set_ylabel('O/F')
    fig.savefig('Part2_OF.pdf',bbox_inches='tight')
    
    # Isp
    fig = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(np.array([80,85,90,95,99]),loc_isp_opt[:],'o-')
    #ax1.plot(np.array([80,85,90,95]),loc_isp[:-1:],'o-')
    ax1.set_xlabel('Peroxide Mass fraction')
    ax1.set_ylabel('Optimal $I_{sp}$')
    fig.savefig('Part2_OptimalIsp.pdf',bbox_inches='tight')

def Part3(): #  Part 2, combined files
    # get all txt files
    s= np.sort(glob.glob("./Part2/*/*.txt"))
    max_values=np.zeros(np.size(s))
    loc_of=np.zeros(np.size(s))
    loc_isp_opt=np.zeros(np.size(s))
    loc_isp=np.zeros(np.size(s))
    for f in range(0,np.size(s)):
        filename=s[f]
        x_data=get_data_dtypes(filename)
        max_values[f] = np.max(x_data['c_star'][2::3])
        for i in range(0,np.size(x_data['c_star'])):
            if (x_data['c_star'][i] == max_values[f] and x_data['isp'][i]!=0):
                loc_of[f]=x_data['of'][i]
                loc_isp_opt[f]=x_data['isp_opt'][i]
                loc_isp[f]=x_data['isp'][i]

    # ratio of Isp/Isp (hybrid/monopropollent)
    fig = plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(np.array([80,85,90,95,99]),(loc_isp_opt[:])/137.69999,'o-')
    ax1.set_xlabel('Peroxide Mass fraction')
    ax1.set_ylabel(r'$\frac{t_2}{t_1}$')
    fig.savefig('Part3.pdf',bbox_inches='tight')




def main():
    #user defined values
    #sys.argv
    global filename1
    filename1 = str(sys.argv[1])

    # example stuff
    #example()

    # Part 1 stuff
    #Part1()

    # part 2 stuff
    #Part2_iter()
    Part2()

    # part 3 stuff
    #Part3()

if __name__ == '__main__':
    main()
