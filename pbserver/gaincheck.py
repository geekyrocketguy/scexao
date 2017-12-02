import numpy as np
import pyfits
import matplotlib.pyplot as plt
import pdb

def gain():
    im_names=['gain_0.890-5-0.fits', 
              'gain_0.939-6-0.fits',
              'gain_0.967-3-0.fits',
              'gain_1.022-7-0.fits']

    v = np.array([0.89, 0.939, 0.967, 1.022])

    increments=np.arange(10)*32

    for j in range(32):
        adus=np.zeros(4)
        for i in range(4):
            img = pyfits.getdata(im_names[i])
            adus[i] = np.median(img[: , :, increments+j])
            
        plt.plot(v, adus, 'ro')
        coeffs=np.polyfit(v, adus, 1)
        p = np.poly1d(coeffs)
        plt.plot(v, p(v), 'b-')
        #plt.title(str(j)+' ' + str(coeffs))
        #plt.show()

        print j, np.round(coeffs)+242e3

    plt.title('Voltage Gain')
    plt.xlabel('Volts')
    plt.ylabel('ADU')
    plt.show()
