import numpy as np
import matplotlib.pyplot as plt
import time

def follow(thefile):
    while True:
        line = np.genfromtxt(thefile)
        ax.plot(line[:,0],line[:,1])
        ax.draw()
        time.sleep(0.1)

def plot_cont(thefile):
    # plt.axis([0, 1, 0, 1])
    plt.ion()
    plt.show()
    plt.xlabel('x')
    plt.ylabel('y')

    while True:
        line = np.genfromtxt(logfile)
        plt.plot(line[:,1],line[:,2],'bo')
        plt.draw()
        time.sleep(1.)
    

logfile = "../cpp_stuff/output_file.txt"
plot_cont(logfile)
#plt.figure()
#ax = fig.add_subplot(111)
#plt.ion()
#plt.figure()
#plt.show()
#time.sleep(10)
#while True:
#line = np.genfromtxt(logfile)
#plt.plot(line[:,0],line[:,1],'bo')
#plt.draw()
#plt.draw()
#time.sleep(10.1)
#loglines = follow(logfile)
#for line in loglines:
#    print line



#import time
#import numpy as np
#import matplotlib.pyplot as plt

# plt.axis([0, 1, 0, 1])
# plt.ion()
# plt.show()
# 
# while True:
    # line = np.genfromtxt(logfile)
    # plt.plot(line[:,0],line[:,1],'bo')
    # plt.draw()
    # time.sleep(0.05)
# for i in range(1000):
    # plt.scatter(1000*x, 1000*y)
    # plt.draw()
 #   time.sleep(0.05)
