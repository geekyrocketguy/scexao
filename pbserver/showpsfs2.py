import matplotlib.pyplot as plt
import numpy as np
import pyfits
import pdb
from matplotlib.colors import LogNorm

def main():
    dir = '/media/data/20170531/saphira/processed/'
    p = 'pbimage_12:17:12.892986502_p.fits'
    aligned = 'pbimage_12:17:12.892986502_p_aligned.fits'

    img_p = pyfits.getdata(dir+p)
    img_aligned = pyfits.getdata(dir+aligned)

    #0.0107''/pixel plate scale
    #lambda/D in units of pixels
    lambda_D_pix = 1.63e-6 / 8.2 * 206265. / .0107 
    lambda_D_min = 4 #inner radius of annulus
    lambda_D_max = 10 #inner radius of annulus
    theta = np.arange(0, 2.*np.pi, 0.01)
    x = np.cos(theta)
    y = np.sin(theta)
    
    x0 = 64.44
    y0 = 62.54

    plt.imshow(np.mean(img_p, axis=0), interpolation='none', norm=LogNorm(), vmin=400)
    plt.colorbar()#shrink=0.8)
    plt.title('Unaligned Mean Image')
    plt.xlim(0, np.shape(img_p)[2]-1)
    plt.ylim(0, np.shape(img_p)[1]-1)
    plt.show()

    plt.imshow(np.mean(img_aligned,axis=0), interpolation='none', norm=LogNorm(), vmin=400)
    plt.plot(lambda_D_min*lambda_D_pix*x + x0, lambda_D_min*lambda_D_pix*y + y0, 'r-')
    plt.plot(lambda_D_max*lambda_D_pix*x + x0, lambda_D_max*lambda_D_pix*y + y0, 'r-')
    plt.colorbar()#shrink=0.8)
    plt.title('Aligned Mean Image')
    plt.xlim(0, np.shape(img_aligned)[2]-1)
    plt.ylim(0, np.shape(img_aligned)[1]-1)
    plt.show()

    plt.imshow(img_aligned[1000], interpolation='none', norm=LogNorm(), vmin=400)
    plt.plot(lambda_D_min*lambda_D_pix*x + x0, lambda_D_min*lambda_D_pix*y + y0, 'r-')
    plt.plot(lambda_D_max*lambda_D_pix*x + x0, lambda_D_max*lambda_D_pix*y + y0, 'r-')
    plt.colorbar()#shrink=0.8)
    plt.title('Aligned Individual Image')
    plt.xlim(0, np.shape(img_aligned)[2]-1)
    plt.ylim(0, np.shape(img_aligned)[1]-1)
    plt.show()

