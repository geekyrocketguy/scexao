import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.interpolation
from scipy.ndimage.filters import gaussian_filter
from skimage.feature import register_translation
import pdb
import os.path

#Corrects for tip/tilt using a cross-correlation. 

def main(file, recalculate=False):
    img = pyfits.getdata(file)#[0:1000]
    before = np.sum(img, 0) #sum of image. For before/after alignment comparison.

    if recalculate or not os.path.isfile(file.replace('.fits', '_offsets.csv') ):
        print "Calculating image shifts."

        shifts = np.zeros((2, len(img)-1))

        for i in range(len(img)-1):
            shifts[:,i] = register_translation(img[i], img[i+1], 1000)[0]

            if i%(len(img)/100)==0:
                print int(round(float(i) / len(img) * 100.)), "% analyzed."

        np.savetxt(file.replace('.fits', '_offsets.csv'), shifts, delimiter=",")
    else:
        print "Restoring image shifts from file."
        shifts = np.loadtxt(file.replace('.fits', '_offsets.csv'), delimiter=',')

    if len(img) != np.size(shifts)/2+1:
        print "Image size does not match shifts file!"
        return

    x_shift_arr_filtered = gaussian_filter(shifts[1,:], 1)
    y_shift_arr_filtered = gaussian_filter(shifts[0,:], 1)
    x_shift_arr_cum = np.cumsum(x_shift_arr_filtered)
    y_shift_arr_cum = np.cumsum(y_shift_arr_filtered)

    coeffs_x = np.polyfit(range(len(x_shift_arr_cum)), x_shift_arr_cum, 3)
    coeffs_y = np.polyfit(range(len(y_shift_arr_cum)), y_shift_arr_cum, 3)
    fit_x = np.poly1d(coeffs_x)
    fit_y = np.poly1d(coeffs_y)
    driftfit_x = fit_x(range(len(x_shift_arr_cum)))
    driftfit_y = fit_y(range(len(y_shift_arr_cum)))
        
    plt.figure(1, figsize=(10,10), dpi=100)
    plt.subplot(211)
    #plt.plot(shifts[0,:]-0.03, label='Unfiltered Y')
    #plt.plot(shifts[1,:]+0.03, label='Unfiltered X')
    #plt.plot(x_shift_arr_filtered + 0.03, label='Filtered X')
    #plt.plot(y_shift_arr_filtered - 0.03, label='Unfiltered Y')
    #plt.title("Shift from one frame to next")
    #plt.legend()

    plt.plot(x_shift_arr_cum, label='X shift')
    plt.plot(y_shift_arr_cum, label='Y shift')
    plt.plot(driftfit_x, label='X fit')
    plt.plot(driftfit_y, label='Y fit')
    plt.title('Cumulative Pixel Shift')
    plt.legend()

    x_shift_arr_cum -= driftfit_x
    y_shift_arr_cum -= driftfit_y

    plt.subplot(212)
    plt.plot(x_shift_arr_cum, label='X shift')
    plt.plot(y_shift_arr_cum, label='Y shift')
    plt.title('Drift-corrected Cumulative Pixel Shift')
    plt.legend()
    plt.show()

    for i in np.arange(len(x_shift_arr_filtered))+1:
        #check shifting
        if 0:
            if i < len(img)-1:
                plt.figure(1, figsize=(12, 5), dpi=100)
                plt.subplot(121)
                mymin = np.sort((img[i-1]-img[i]).flatten())[0.001*np.size(img[i])]
                mymax = np.sort((img[i-1]-img[i]).flatten())[0.999*np.size(img[i])]
                plt.imshow(img[i-1] - img[i], interpolation='none', vmin=mymin, vmax=mymax)
                plt.title(np.std(img[i-1] - img[i], ddof=1))
                
                plt.subplot(122)
                shifted = scipy.ndimage.interpolation.shift(img[i], \
                                                            [1*y_shift_arr_filtered[i-1],\
                                                             1*x_shift_arr_filtered[i-1]],\
                                                            mode='wrap')
                plt.imshow(img[i-1] - shifted, interpolation='none', vmin=mymin, vmax=mymax)
                plt.title(np.std(img[i-1] - shifted, ddof=1))
                plt.show()


        img[i] = scipy.ndimage.interpolation.shift(img[i], [1*y_shift_arr_cum[i-1],\
                                                            1*x_shift_arr_cum[i-1]],\
                                                   mode='wrap')

        if i%(len(img)/100)==0:
            print int(round(float(i) / len(img) * 100.)), "% aligned."

            
    #Check if images were shifted correctly
    if 1:
        print "Calculating image shifts."
        shifts = np.zeros((2, len(img)-1))
        
        for i in range(len(img)-1):
            shifts[:,i] = register_translation(img[i], img[i+1], 100)[0]
            
            if i%(len(img)/100)==0:
                print int(round(float(i) / len(img) * 100.)), "% analyzed."

        plt.plot(shifts[0,:]-0.03, label='Unfiltered Y')
        plt.plot(shifts[1,:]+0.03, label='Unfiltered X')
        plt.legend()
        plt.show()

    plt.figure(1, figsize=(12, 5), dpi=100)
    plt.subplot(121)
    plt.imshow(before, interpolation='none')
    plt.colorbar()
    plt.title("Before")
    plt.subplot(122)
    plt.imshow(np.sum(img, 0), interpolation='none')
    plt.colorbar()
    plt.title("After")
    plt.show()

    if 1: #Save image?
        newfilename = file[:file.find('.fits')] + '_aligned.fits'
        print "Saving image as "+newfilename
        pyfits.writeto(newfilename, img, clobber='true') #save file
        print "Saved."
