import numpy as np
import matplotlib.pyplot as plt
import time

def follow(thefile):
    while True:
        line = np.genfromtxt(thefile)
        ax.cla()
        ax.plot(line[:,1],line[:,2],'b-')
        ax.draw()
        time.sleep(2.1)

def plot_cont(thefile,thefile2):
    # plt.axis([0, 1, 0, 1])
    plt.ion()
    plt.show()

    while True:
        line = np.genfromtxt(thefile)
        line2 = np.genfromtxt(thefile2)
        plt.cla()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.plot(line2[:,1],line2[:,2],'ro',mew='0',label='Coast')
        plt.plot(line[:,1],line[:,2],'bo',mew='0',label='Continuous Thrust')
        plt.legend(loc='center left',bbox_to_anchor=(1,0.815),numpoints=1)
        plt.draw()
        time.sleep(1.1)
    
def plot_cont2(thefile,thefile2):
    plt.figure(figsize=(3,3))

    while True:
        line = np.genfromtxt(thefile)
        line2 = np.genfromtxt(thefile2)
        plt.cla()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.plot(line2[:,1],line2[:,2],'r.',mew=0,label='Coast')
        plt.plot(line[:,1],line[:,2],'b.',mew=0,label='Thrust')
        plt.plot(0,0,'k.',mew=0,label='Earth')
        plt.legend(loc='upper left',bbox_to_anchor=(1,0.815),numpoints=1)
        #plt.legend(loc='best',numpoints=1)
        plt.savefig('file.png',bbox_inches='tight')
        #plt.savefig('file.pdf',bbox_inches='tight')
        time.sleep(3.1)

def plot_once(thefile,thefile2):
    plt.figure(figsize=(3,3))
    plt.axis('equal')
    line = np.genfromtxt(thefile)
    line2 = np.genfromtxt(thefile2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(line2[:,1],line2[:,2],'r.',mew=0,label='Coast')
    plt.plot(line[:,1],line[:,2],'b.',mew=0,label='Thrust')
    circle1 = plt.Circle((0,0),6371393.,color='c',label='Earth')
    plt.gca().add_artist(circle1)
    plt.legend(loc='upper left',bbox_to_anchor=(1,0.815),numpoints=1)
    #plt.legend(loc='best',numpoints=1)
    #plt.savefig('file.png',bbox_inches='tight')
    plt.savefig('file.pdf',bbox_inches='tight')
def follow_2(thefile):
    # thefile.seek(0,2)
    while True:
        line = thefile.readline()
        if not line:
            time.sleep(2.1)
            continue
        yield np.matrix(line)



continuous = "../cpp_stuff/continuous.txt"
coast = "../cpp_stuff/coast.txt"
#   logfile = open(continuous,'r')
#plot_cont2(continuous,coast)
plot_once(continuous,coast)
#   plt.ion()
#   plt.show()
#   plt.xlabel('x')
#   plt.ylabel('y')
#   # initial=np.genfromtxt(filename)
#   # plt.plot(initial[:,1],initial[:,2],'bo')
#   loglines = follow_2(logfile)
#   #loglines = follow(filename)
#   for line in loglines:
#       # print line
#       #data = np.genfromtxt(line.strip())
#       plt.plot(line[0,1],line[0,2],'b.')
#       plt.draw()
#       #time.sleep(1.1)
