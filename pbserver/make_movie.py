#Takes a fits file and saves a video from it

import pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
import pdb
#from scipy.stats import mode
import os
#import time

dir = 'images/'
filename = 'saphira_14:46:51.541667336_cleanedaligned.fits'
print "Reading image..."
img = pyfits.getdata(dir+filename)[:600]#2500:3500]
print "Image read. Size is ", np.shape(img)

#make image logarithmic
#loc = np.where(img > 0)
#img[loc] = np.log(img[loc])
#loc = np.where(img < 0)
#img[loc] = -1.*np.log(-1*img[loc])
img += 10
loc = np.where(img < 1)
img[loc] = 1

i=0
while np.std(img[i], ddof=0) ==0: #don't use a blank image for min/max
    i+= 1
mymin = np.sort(img[i].flatten())[np.size(img[i])*.01]
mymax = np.sort(img[i].flatten())[np.size(img[i])*.999]
counter = 0
for z in range(np.shape(img)[0]):
    if z % ((np.shape(img)[0])/10) == 0:
        print str(int(round(100. * z / np.shape(img)[0]))).replace(' ',''), \
            '% done.'
        
    if np.std(img[z]) != 0: #if it's not a blank image
        counter += 1
        fig = plt.figure()
        ax = plt.Axes(fig, [0.04, -0.05, .95, 1.1])
        fig.add_axes(ax)
        plt.imshow(img[z], interpolation='none', norm=LogNorm(), vmin=5, vmax=mymax)
        #plt.show()
        #pdb.set_trace()
        #plt.colorbar()
        plt.savefig('imagestovid/img%03d.png' % (counter,))
        plt.clf()


os.system('avconv -r 30 -i imagestovid/img%03d.png -b:v 1000k -c:v libx264 speckles_ao188.mp4') #30 fps. %03d means 3 digit numerical extension of filenames
