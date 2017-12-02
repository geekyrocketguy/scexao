#!/usr/bin/env python

#Takes as input a filename and two frame numbers. Displays a CDS image.

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pdb
import sys
import time
   
args=sys.argv

if len(args) != 4: 
    print
    print "Usage: python cds.py myfilename.fits firstframe lastframe"
    print "Displays a cds frame of the image."
    print
    #return
else:
    filename=str(args[1])
    f1 = int(args[2])
    f2 = int(args[3])


#def cds(filename, f1, f2):
img = pyfits.getdata(filename)
cds = img[f1] - img[f2]
if np.median(cds) < 0: cds *= -1

print "stddev is:", np.std(cds, ddof=1)

mymin=np.sort(cds.flatten())[0.01*np.size(cds)]
mymax=np.sort(cds.flatten())[0.9999*np.size(cds)]

plt.imshow(cds, interpolation='none', vmin=mymin, vmax=mymax)
#plt.imshow(cds, interpolation='none', norm=LogNorm(), vmin=2000, vmax=mymax)
plt.title(filename + ' CDS ' + str(f2) + " - " + str(f1))
plt.colorbar()
plt.show()

stddevs=np.array([])
for i in range(10, np.shape(img)[0], 2):
    stddevs=np.append(stddevs, np.std(img[i+1]-img[i], ddof=1) )

print "Average RFI is:", np.median(stddevs)
time.sleep(5)
