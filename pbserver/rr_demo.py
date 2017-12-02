import pyfits
import matplotlib.pyplot as plt
import numpy as np

def demo():
    data_img = pyfits.getdata('/media/data/20170113/saphira/pbimage_22:05:02.538868029.fits')

    dark_img = pyfits.getdata('/media/data/20170113/saphira/pbimage_22:04:54.693597787.fits')

    median_img = np.median(dark_img[:, 128:, :], 0)

    stddev_rrr = np.array([])
    stddev_rr = np.array([])

    for i in range(np.shape(data_img)[0]):
        stddev_rrr = np.append(stddev_rrr, 
                               np.std(data_img[i,128:,:] - data_img[i,:128,:], ddof=1) )

        stddev_rr = np.append(stddev_rr, 
                               np.std(data_img[i,128:,:] - median_img, ddof=1) )


    print "median RFI for RRR image:", np.median(stddev_rrr)
    print "median RFI for RR image:", np.median(stddev_rr)

    for i in range(0,1000, 100):
        rrr_img = data_img[i,128:,:] - data_img[i,:128,:]
        rr_img  = data_img[i,128:,:] - median_img


        plt.figure(num=1, figsize=(10, 5), dpi=100) 

        plt.subplot(121)
        plt.imshow(rrr_img, interpolation='none', vmin=-1000, vmax=1000)
        #plt.colorbar()
        plt.title("RRR subtraction, stddev="+str(np.std(rrr_img, ddof=1))[:5])

        plt.subplot(122)
        plt.imshow(rr_img, interpolation='none', vmin=-1000, vmax=1000)
        #plt.colorbar()
        plt.title("RR subtraction, stddev="+str(np.std(rr_img, ddof=1))[:5])

        plt.show()
