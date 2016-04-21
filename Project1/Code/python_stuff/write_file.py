import numpy as np
import matplotlib.pyplot as plt
import time

for i in range(10000):
    logfile=open("file.txt",'a')
    # logfile.write("%i %i\n" % (i,2*i))
    logfile.write("%f %f\n" % (np.random.random(),np.random.random()))
    time.sleep(1.05) 
    logfile.close()
