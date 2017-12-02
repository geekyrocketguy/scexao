import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.interpolation
#from scipy.ndimage.filters import gaussian_filter
import pdb

#Corrects for tip/tilt in no-ao images. Was forked from tiptilt3.
#writes shifted pixels as -7777

def align(file):
    img = pyfits.getdata(file)
    
#### FIND AVERAGE position
    #Calculate approximate center of PSF
    im_coadd = (np.sum(img,0)/len(img)).astype('float')

    im_coadd -= np.median(im_coadd)
    im_coadd[im_coadd<0] = 0

    x0 = np.sum(np.sum(im_coadd, 0)*range(np.shape(im_coadd)[1])) / \
        np.sum(im_coadd)
    y0 = np.sum(np.sum(im_coadd, 1)*range(np.shape(im_coadd)[0])) / \
        np.sum(im_coadd)

#### FIND INDIVIDUAL POSITIONS
    for i in range(len(img)):
        if i%(round(len(img)/100.))==0: #update status
            print str(int(round(float(i)/len(img)*100.)))+"% complete."
        im = np.copy(img[i])
        im -= np.median(im)
        im[im<0] = 0

        x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
        y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)

        img[i] = scipy.ndimage.interpolation.shift(img[i], [y0-y1, x0-x1], \
                                                   mode='mirror')

##### SAVE IMAGE
    if 1:
        newfilename = file[:file.find('.fits')] + '_aligned.fits'
        print "Saving image as "+newfilename
        pyfits.writeto(newfilename, img, clobber='true') #save file
        print "Saved."
