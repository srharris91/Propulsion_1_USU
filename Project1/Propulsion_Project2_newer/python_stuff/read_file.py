import numpy as np
import matplotlib.pyplot as plt
import time

def follow(thefile):
    thefile.seek(0,2) 
    while True:
        line = thefile.readline()
        # print "line says ",line
        if not line:
            time.sleep(2.1)
            continue
        yield line


logfile = open("file.txt",'r')
loglines = follow(logfile)
for line in loglines:
    print line
