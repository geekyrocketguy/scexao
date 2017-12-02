import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage.interpolation
import pdb

def align():
    dir = '/media/data/20170313/saphira/processed/'
    filename = 'pbimage_12:10:16.115625254_p.fits'
    img = pyfits.getdata(dir+filename)[0:1000]

#### FIND AVERAGE LOCATION OF EACH ASTROMETRIC SPECKLE
    
    #Calculate approximate center of PSF
    im_coadd = (np.sum(img,0)/len(img)).astype('float')

    x0 = np.sum(np.sum(im_coadd, 0)*range(np.shape(im_coadd)[1])) / \
        np.sum(im_coadd)
    y0 = np.sum(np.sum(im_coadd, 1)*range(np.shape(im_coadd)[0])) / \
        np.sum(im_coadd)

    n_brightest = 15 #number of pixels to include in centroid calculation

    #Q2 Q1
    #Q3 Q4

    #Q1
    im = np.copy(im_coadd)
    im[y0 - 25 : , :] = 0
    im[ : , : x0 + 25 ] = 0
    threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
    im -= threshold
    im[im<0] = 0
    q1_x0 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
    q1_y0 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)


    #Q2
    im = np.copy(im_coadd)
    im[y0 - 25 : , :] = 0
    im[ : , x0 - 25 : ] = 0
    threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
    im -= threshold
    im[im<0] = 0
    q2_x0 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
    q2_y0 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)


    #Q3
    im = np.copy(im_coadd)
    im[ : y0 + 25 , :] = 0
    im[ : , x0 - 25 : ] = 0
    threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
    im -= threshold
    im[im<0] = 0
    q3_x0 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
    q3_y0 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)


    #Q4
    im = np.copy(im_coadd)
    im[ : y0 + 25 , :] = 0
    im[ : , : x0 + 25 ] = 0
    threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
    im -= threshold
    im[im<0] = 0
    q4_x0 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
    q4_y0 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)

    if 0:
        plt.imshow(im_coadd, interpolation='none')
        plt.plot([x0] , [y0], 'wx')
        plt.plot([q1_x0] , [q1_y0], 'wx')
        plt.plot([q2_x0] , [q2_y0], 'wx')
        plt.plot([q3_x0] , [q3_y0], 'wx')
        plt.plot([q4_x0] , [q4_y0], 'wx')
        plt.show()


### ALIGN IMAGES

    #make masks around each astrometric speckle
    mask_q1 = np.zeros(np.shape(im))
    mask_q2 = np.zeros(np.shape(im))
    mask_q3 = np.zeros(np.shape(im))
    mask_q4 = np.zeros(np.shape(im))

    for x in range(np.shape(im)[1]):
        for y in range(np.shape(im)[0]):
            if (x-q1_x0)**2 + (y-q1_y0)**2 < 10**2: mask_q1[y,x] = 1
            if (x-q2_x0)**2 + (y-q2_y0)**2 < 10**2: mask_q2[y,x] = 1
            if (x-q3_x0)**2 + (y-q3_y0)**2 < 10**2: mask_q3[y,x] = 1
            if (x-q4_x0)**2 + (y-q4_y0)**2 < 10**2: mask_q4[y,x] = 1


    #Recenter PSFs to reduce tip/tilt
    brightness_threshold = 7000 #when do you stop looking for an astrometric speckle?
    for i in range(0, np.shape(img)[0]):
        if i % 100 ==0:
            print int(np.round(np.float(i)/np.shape(img)[0]*100)), "% done."

        if 0:
            plt.imshow(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none')
            plt.show()

        #Q1
        im = np.copy(img[i])
        im *= mask_q1
        if np.max(im) > brightness_threshold:
            threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
            im -= threshold
            im[im<0] = 0
            #plt.imshow(im, interpolation='none')
            #plt.show()
            q1_x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
            q1_y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)
        else:
            q1_x1 = -1
            q1_y1 = -1

        #Q2
        im = np.copy(img[i])
        im *= mask_q2
        if np.max(im) > brightness_threshold:
            threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
            im -= threshold
            im[im<0] = 0
            q2_x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
            q2_y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)
        else:
            q2_x1 = -1
            q2_y1 = -1

        #Q3
        im = np.copy(img[i])
        im *= mask_q3
        if np.max(im) > brightness_threshold:
            threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
            im -= threshold
            im[im<0] = 0
            q3_x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
            q3_y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)        
        else:
            q3_x1 = -1
            q3_y1 = -1

        #Q4
        im = np.copy(img[i])
        im *= mask_q4
        if np.max(im) > brightness_threshold:
            threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
            im -= threshold
            im[im<0] = 0
            q4_x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
            q4_y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)
        else:
            q4_x1 = -1
            q4_y1 = -1

        #calculate shift amount
        x_arr = np.array([])
        y_arr = np.array([])

        if q1_x1 != -1:
            x_arr = np.append(x_arr, q1_x1-q1_x0)
            y_arr = np.append(y_arr, q1_y1-q1_y0)
        
        if q2_x1 != -1:
            x_arr = np.append(x_arr, q2_x1-q2_x0)
            y_arr = np.append(y_arr, q2_y1-q2_y0)
        
        if q3_x1 != -1:
            x_arr = np.append(x_arr, q3_x1-q3_x0)
            y_arr = np.append(y_arr, q3_y1-q3_y0)
        
        if q4_x1 != -1:
            x_arr = np.append(x_arr, q4_x1-q4_x0)
            y_arr = np.append(y_arr, q4_y1-q4_y0)
        
        if len(y_arr) == 0:
            print "You giving me empty images?!"
            print "max speckle:", np.max(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4))
            plt.imshow(img[i], vmin=0, vmax=brightness_threshold, interpolation='none')
            plt.colorbar()
            plt.show()
            pdb.set_trace()


        #check shifting
        if 0:
            plt.figure(1, figsize=(15, 5), dpi=100) 

            plt.subplot(131)
            plt.imshow(im_coadd, interpolation='none')
            plt.plot([q1_x0, q2_x0, q3_x0, q4_x0], [q1_y0, q2_y0, q3_y0, q4_y0], 'wx')
        
            plt.subplot(132)
            plt.imshow(img[i] * (mask_q1 + mask_q2 + mask_q3 + mask_q4), interpolation='none')
            plt.plot(x_arr + [q1_x0, q2_x0, q3_x0, q4_x0], y_arr + [q1_y0, q2_y0, q3_y0, q4_y0], 'kx')
                        
            crop = scipy.ndimage.interpolation.shift(
                img[i], [-1*np.median(y_arr), -1*np.median(x_arr)], mode='wrap')
            plt.subplot(133)
            plt.imshow(crop, interpolation='none')
            plt.plot([q1_x0, q2_x0, q3_x0, q4_x0], [q1_y0, q2_y0, q3_y0, q4_y0], 'wx')

            plt.show()

        im = np.copy(img[i])
        img[i,:,:] = scipy.ndimage.interpolation.shift(im, [-1*np.median(y_arr), -1*np.median(x_arr)], mode='wrap')

    plt.figure(1, figsize=(10, 5), dpi=100) 
    plt.subplot(121)
    plt.imshow(im_coadd, interpolation='none')
    plt.subplot(122)
    plt.imshow((np.sum(img,0)/len(img)).astype('float'), interpolation='none')
    plt.show()

    
    if 1: #Save image?
        newfilename = filename[:filename.find('.fits')] + 'aligned.fits'
        print "Saving image as "+'important_images/'+newfilename
        pyfits.writeto('important_images/'+newfilename, img, clobber='true') #save file
        print "Saved."
