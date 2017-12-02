#!/usr/bin/env python

import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt

os.system('shmim2fits pbimagediff recentimage.fits')
img = pyfits.getdata('recentimage.fits')
img = img + (img < -1*2**12)*2**16 #this is necessary to recast type

plt.figure(num=1, figsize=(11, 4), dpi=100) 
plt.subplot(121)
plt.imshow(img, interpolation='none')
plt.colorbar()

plt.subplot(122)
plt.imshow(img, interpolation='none', vmin=-100, vmax=500)
plt.colorbar()

plt.show()


