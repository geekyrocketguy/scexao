#!/usr/bin/env python

#Takes as input a filename and two frame numbers. Displays a CDS image.

#import pyfits
from astropy.io import fits as pyfits
import numpy as np
import sys

args=sys.argv

if len(args) != 4: 
    print()
    print( "Usage: python trimfits.py myfilename.fits firstframe lastframe")
    print('0=first frame')
    print( "Crops fits cube to user-specified frames.")
    print()
    #return
else:
    filename=str(args[1])
    f1 = int(args[2])
    f2 = int(args[3])

img = pyfits.getdata(filename)
img = img[f1:f2]

print("New image size is "+str(np.shape(img)))
pyfits.writeto(filename[:-5]+"_trimmed.fits", img, clobber='True')
