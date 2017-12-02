import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.interpolation
from scipy.ndimage.filters import gaussian_filter
import pdb

#Corrects for tip/tilt using oscillating astrometric speckles, as were used
# for the 5/31 observing run. 
# v1 is on the local computer and uses the PSF core.
# v2 is on scexao and uses crappily modulated speckles

def align(file):
    #dir = '/media/data/20170313/saphira/processed/'
    #filename = 'pbimage_12:10:16.115625254_p.fits'
    img = pyfits.getdata(file)[:-1]
    
    n_brightest = 17 #number of pixels to include in centroid calculation
    brightness_threshold = 1000 #when do you stop looking for an astrometric speckle?
    #1600 for 12:16:42

    x_shift_arr = np.zeros(len(img))
    y_shift_arr = np.zeros(len(img))

#### FIND AVERAGE LOCATION OF EACH ASTROMETRIC SPECKLE
    #Calculate approximate center of PSF
    im_coadd = (np.sum(img,0)/len(img)).astype('float')

    x0 = np.sum(np.sum(im_coadd, 0)*range(np.shape(im_coadd)[1])) / \
        np.sum(im_coadd)
    y0 = np.sum(np.sum(im_coadd, 1)*range(np.shape(im_coadd)[0])) / \
        np.sum(im_coadd)

    #Get ALL PSFs sorta close to this
    if 0:
        print "coarse aligning"
        x_mask_shift = np.zeros(len(img))
        y_mask_shift = np.zeros(len(img))
        for z in range(len(img)):
            x_mask_shift[z] = np.sum(np.sum(img[z], 0)*range(np.shape(img[z])[1])) / \
                              np.sum(img[z]) - x0
            y_mask_shift[z] = np.sum(np.sum(img[z], 1)*range(np.shape(img[z])[0])) / \
                              np.sum(img[z]) - y0
            
            #img[z] = np.roll(img[z], int(round(x0-x1)), axis=1)
            #img[z] = np.roll(img[z], int(round(y0-y1)), axis=0)

        if 0:
            plt.figure(1, figsize=(10, 5), dpi=100) 
            plt.subplot(121)
            plt.imshow(im_coadd, interpolation='none', vmin=1, vmax=np.max(im_coadd)*1.1)
            plt.title('before')
            plt.subplot(122)
            plt.imshow((np.sum(img,0)/len(img)).astype('float'), interpolation='none', vmin=1,
                       vmax=np.max(im_coadd)*1.1)
            plt.title('after')
            plt.show()

    #Q2 Q1
    #Q3 Q4

    #Q1
    im = np.copy(im_coadd)
    im[y0 - 25 : , :] = 0
    im[ : , : x0 + 25 ] = 0
    q1_x0, q1_y0 = centroid(im, n_brightest+10)

    #Q2
    im = np.copy(im_coadd)
    im[y0 - 25 : , :] = 0
    im[ : , x0 - 25 : ] = 0
    q2_x0, q2_y0 = centroid(im, n_brightest+10)

    #Q3
    im = np.copy(im_coadd)
    im[ : y0 + 25 , :] = 0
    im[ : , x0 - 25 : ] = 0
    q3_x0, q3_y0 = centroid(im, n_brightest+10)

    #Q4
    im = np.copy(im_coadd)
    im[ : y0 + 25 , :] = 0
    im[ : , : x0 + 25 ] = 0
    q4_x0, q4_y0 = centroid(im, n_brightest+10)


    print "x avg", np.mean((q1_x0, q2_x0, q3_x0, q4_x0))
    print "y avg", np.mean((q1_y0, q2_y0, q3_y0, q4_y0))

    if 1:
        plt.imshow(im_coadd, interpolation='none', vmin=0, vmax=5000)
        plt.plot([x0] , [y0], 'wx')
        plt.plot([q1_x0] , [q1_y0], 'kx')
        plt.plot([q2_x0] , [q2_y0], 'kx')
        plt.plot([q3_x0] , [q3_y0], 'kx')
        plt.plot([q4_x0] , [q4_y0], 'kx')
        plt.show()


### ALIGN IMAGES

    #make masks around each astrometric speckle
    mask_q1_avg = np.zeros(np.shape(im))
    mask_q2_avg = np.zeros(np.shape(im))
    mask_q3_avg = np.zeros(np.shape(im))
    mask_q4_avg = np.zeros(np.shape(im))

    for x in range(np.shape(im)[1]):
        for y in range(np.shape(im)[0]):
            if (x-q1_x0)**2 + (y-q1_y0)**2 < 11**2: mask_q1_avg[y,x] = 1
            if (x-q2_x0)**2 + (y-q2_y0)**2 < 11**2: mask_q2_avg[y,x] = 1
            if (x-q3_x0)**2 + (y-q3_y0)**2 < 11**2: mask_q3_avg[y,x] = 1
            if (x-q4_x0)**2 + (y-q4_y0)**2 < 11**2: mask_q4_avg[y,x] = 1


    #Recenter PSFs to reduce tip/tilt
    for i in range(0, np.shape(img)[0]):
        if i % 100 ==0:
            print int(np.round(np.float(i)/np.shape(img)[0]*100)), "% done."

        if 0:
            plt.imshow(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none')
            plt.show()

        if i == 0:
            diff = img[0] - img[1]
        elif (i>0) & (i<len(img)-1):
            diff = img[i+1] - img[i-1]
        else: #if i==len(img)-1
            diff = img[i] - img[i-1]
        #n_good_quads = 0

        #shift masks to compensate for crazy tip/tilt, but only if they exist
        if 'x_mask_shift' in vars():
            mask_q1 = np.roll(mask_q1_avg, int(round(x_mask_shift[i])), axis=1)
            mask_q1 = np.roll(mask_q1, int(round(y_mask_shift[i])), axis=0)
            mask_q2 = np.roll(mask_q2_avg, int(round(x_mask_shift[i])), axis=1)
            mask_q2 = np.roll(mask_q2, int(round(y_mask_shift[i])), axis=0)
            mask_q3 = np.roll(mask_q3_avg, int(round(x_mask_shift[i])), axis=1)
            mask_q3 = np.roll(mask_q3, int(round(y_mask_shift[i])), axis=0)
            mask_q4 = np.roll(mask_q4_avg, int(round(x_mask_shift[i])), axis=1)
            mask_q4 = np.roll(mask_q4, int(round(y_mask_shift[i])), axis=0)
        else:
            mask_q1 = mask_q1_avg
            mask_q2 = mask_q2_avg
            mask_q3 = mask_q3_avg
            mask_q4 = mask_q4_avg

        if 0:#i%100==0:
            plt.imshow(img[i] + 1000*(mask_q1+mask_q2+mask_q3+mask_q4), interpolation='none',
                       vmin=0, vmax=1e4)
            plt.plot([y0], [x0], 'wx')
            plt.plot([y0+y_mask_shift[i]], [x0+x_mask_shift[i]], 'kx')
            plt.plot([q1_x0] , [q1_y0], 'wx')
            plt.plot([q2_x0] , [q2_y0], 'wx')
            plt.plot([q3_x0] , [q3_y0], 'wx')
            plt.plot([q4_x0] , [q4_y0], 'wx')
            plt.colorbar()
            plt.show()

        x_arr, y_arr, weight = get_centroids(diff, mask_q1, mask_q2, mask_q3, mask_q4, \
                                             brightness_threshold, n_brightest, q1_x0, \
                                             q1_y0, q2_x0, q2_y0, q3_x0, q3_y0, q4_x0, q4_y0)
        #print "1. len x_arr, y_arr", len(x_arr), len(y_arr)
        #IF NOT ENOUGH POINTS FOUND          
        #if  0: #(len(x_arr) < 3) & (i>0) & (i<len(img)-1):
        if  (i>0) & (i<len(img)-1):
            diff = img[i] - img[i-1]
            x_arr2, y_arr2, weight2 = get_centroids(diff, mask_q1, mask_q2, mask_q3, mask_q4, \
                                                    brightness_threshold, n_brightest, q1_x0, \
                                                    q1_y0, q2_x0, q2_y0, q3_x0, q3_y0, q4_x0, q4_y0)

            x_arr = np.append(x_arr, x_arr2)
            y_arr = np.append(y_arr, y_arr2)
            weight = np.append(weight, weight2)
            #print "2. len x_arr, y_arr", len(x_arr), len(y_arr)
            
            diff = img[i] - img[i+1]
            x_arr2, y_arr2, weight2 = get_centroids(diff, mask_q1, mask_q2, mask_q3, mask_q4, \
                                                    brightness_threshold, n_brightest, q1_x0, \
                                                    q1_y0, q2_x0, q2_y0, q3_x0, q3_y0, q4_x0, q4_y0)

            x_arr = np.append(x_arr, x_arr2)
            y_arr = np.append(y_arr, y_arr2)
            weight = np.append(weight, weight2)
            #print "3. len x_arr, y_arr", len(x_arr), len(y_arr)

        #if i==110: pdb.set_trace()
            
        #compute weighted average

        if len(x_arr)==0:
            print np.max(diff * (mask_q1 + mask_q2 + mask_q3 + mask_q4))
            print np.min(diff * (mask_q1 + mask_q2 + mask_q3 + mask_q4))
            diff = img[i+1] - img[i-1]
            print np.max(diff * (mask_q1 + mask_q2 + mask_q3 + mask_q4))
            print np.min(diff * (mask_q1 + mask_q2 + mask_q3 + mask_q4))
            pdb.set_trace()

        x_shift = np.sum(x_arr * weight) / np.sum(weight)
        y_shift = np.sum(y_arr * weight) / np.sum(weight)

        x_shift_arr[i] = x_shift
        y_shift_arr[i] = y_shift

        #check shifting
        if 0:
            plt.figure(1, figsize=(15, 5), dpi=100) 

            plt.subplot(131)
            plt.imshow(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none')
            plt.plot([q1_x0, q2_x0, q3_x0, q4_x0], [q1_y0, q2_y0, q3_y0, q4_y0], 'wx')
            plt.title('Reference Positions')
        
            plt.subplot(132)
            plt.imshow(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none')
            plt.plot(np.median(x_arr) + np.array([q1_x0, q2_x0, q3_x0, q4_x0]), \
                     np.median(y_arr) + np.array([q1_y0, q2_y0, q3_y0, q4_y0]), 'kx')
            plt.title('Detected Positions')
                        
            crop = scipy.ndimage.interpolation.shift(im, [-1*x_shift, -1*y_shift], mode='wrap')
            plt.subplot(133)
            plt.imshow(crop, interpolation='none', vmin=0, vmax=4000)
            plt.plot([q1_x0, q2_x0, q3_x0, q4_x0], [q1_y0, q2_y0, q3_y0, q4_y0], 'wx')
            plt.title('Shifted to Reference Positions')

            plt.show()

        #print "xarr", x_arr
        #print "yarr", y_arr
        #print "weight", weight
        #print
        
        #im = np.copy(img[i])
        #img[i,:,:] = scipy.ndimage.interpolation.shift(im,
        #                                    [-1*y_shift, -1*x_shift], mode='wrap')

    x_shift_arr_filtered = gaussian_filter(x_shift_arr, 1)
    y_shift_arr_filtered = gaussian_filter(y_shift_arr, 1)
    for i in range(len(img)):
        img[i] = scipy.ndimage.interpolation.shift(img[i], [-1*y_shift_arr_filtered[i],\
                                                            -1*x_shift_arr_filtered[i]],\
                                                   mode='wrap')

    plt.figure(1, figsize=(10, 10), dpi=100) 
    plt.subplot(211)
    plt.plot(x_shift_arr)
    plt.plot(x_shift_arr_filtered)
    plt.title('X displacement')
    plt.subplot(212)
    plt.plot(y_shift_arr)
    plt.plot(y_shift_arr_filtered)
    plt.title('Y displacement')
    plt.show()
    #pdb.set_trace()
        
    plt.figure(1, figsize=(10, 5), dpi=100) 
    plt.subplot(121)
    plt.imshow(im_coadd, interpolation='none', vmin=1)#, vmax=np.max(im_coadd)*1.1)
    plt.colorbar()
    plt.title('before')
    plt.subplot(122)
    plt.imshow((np.sum(img,0)/len(img)).astype('float'), interpolation='none', vmin=1)#,
               #vmax=np.max(im_coadd)*1.1)
    plt.colorbar()
    plt.title('after')
    plt.show()


    if 1: #Save image?
        newfilename = file[:file.find('.fits')] + '_aligned_tt3.fits'
        print "Saving image as "+newfilename
        pyfits.writeto(newfilename, img, clobber='true') #save file
        print "Saved."


def centroid(im, n_brightest):
    #takes as arguments the 2D image and number of pixels to include in the calculation
    
    threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
    im -= threshold
    im[im<0] = 0
    #plt.imshow(im, interpolation='none')
    #plt.show()
    x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
    y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)

    return [x1, y1]


def get_centroids(diff, mask_q1, mask_q2, mask_q3, mask_q4, brightness_threshold, \
                  n_brightest, q1_x0, q1_y0, q2_x0, q2_y0, q3_x0, q3_y0, q4_x0, q4_y0):
    xarr = []
    yarr = []
    wt = []

    #sorta compensate for RFI
    for y in range(np.shape(diff)[0]):
        for x in range(np.shape(diff)[1]/32):
            #if np.max((mask_q1+mask_q2+mask_q3+mask_q4)[y, x*32:(x+1)*32])==1:
            diff[y, x*32:(x+1)*32] -= np.median(diff[y, x*32:(x+1)*32])
            #np.sort(diff[y, x*32:(x+1)*32])[5] #median is too high
            
    
    if 0:
        #plt.imshow(diff*(mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none', \
        #           vmin=-1*brightness_threshold, vmax=brightness_threshold)
        plt.imshow(diff, interpolation='none', \
                   vmin=-500, vmax=500)
        plt.colorbar()
        plt.show()
    
    #Q1
    im = diff*mask_q1
    if abs(np.min(im)) > np.max(im): im *= -1 #invert if negative
    #print np.max(im)
    if np.max(im) > brightness_threshold:
        x1, y1 = centroid(im, n_brightest)
        xarr = np.append(xarr, x1-q1_x0)
        yarr = np.append(yarr, y1-q1_y0)
        wt = np.append(wt, np.sum(im))

        
    #Q2
    im = diff*mask_q2
    if abs(np.min(im)) > np.max(im): im *= -1 #invert if negative
    #print np.max(im)
    if np.max(im) > brightness_threshold:
        x1, y1 = centroid(im, n_brightest)
        xarr = np.append(xarr, x1-q2_x0)
        yarr = np.append(yarr, y1-q2_y0)
        wt = np.append(wt, np.sum(im))

    #Q3
    im = diff*mask_q3
    if abs(np.min(im)) > np.max(im): im *= -1 #invert if negative
    #print np.max(im)
    if np.max(im) > brightness_threshold:
        x1, y1 = centroid(im, n_brightest)
        xarr = np.append(xarr, x1-q3_x0)
        yarr = np.append(yarr, y1-q3_y0)
        wt = np.append(wt, np.sum(im))

    #Q4
    im = diff*mask_q4
    if abs(np.min(im)) > np.max(im): im *= -1 #invert if negative
    #print np.max(im)
    if np.max(im) > brightness_threshold:
        x1, y1 = centroid(im, n_brightest)
        xarr = np.append(xarr, x1-q4_x0)
        yarr = np.append(yarr, y1-q4_y0)
        wt = np.append(wt, np.sum(im))

    return xarr, yarr, wt
