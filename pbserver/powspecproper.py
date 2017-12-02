#started on 2017-04-24
#plots power spectrum of RFI, taking into account bonus clock cycles
#This is the latest version

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb
#from scipy.stats import mode
#import os
#import time
from scipy import signal

welch = True

#Declare some global variables
dir = '../pbserver.main.23APR17/DATA/'
filename = 'noise_utr-1-0.fits'#rfi2-1-0.fits'

t_pix = 1e-6 #seconds per pixel
t_endrow = 540e-9 #seconds at the end of each row
t_endframe = 310e-9 #seconds at the end of each frame
t_clock = 10e-9 #duration of clock cycle
#end of row clock at end of frame?

print "Reading image..."
img = pyfits.getdata(dir+filename)[210:]
#img = np.random.normal(loc=1e4, scale=1000, size=(10, 256, 320))
#img = np.random.uniform(low=100, high=200, size=(10, 256, 320))
img_avg = np.median(img, 0)
for z in range(np.shape(img)[0]):
    img[z] -= img_avg
print "Image read and median subtracted. Size is ", np.shape(img)
print "Median, stddev:", np.median(img), np.std(img, ddof=1)

def main():

    if 0:
        for i in np.arange(0, len(img), 20):
            plt.imshow(img[i] , interpolation='none', vmin=-90, vmax=90)
            plt.colorbar()
            #plt.title("Example median-subtracted frame")
            plt.title(i)
            plt.show()

    row_clocks = np.shape(img)[2]/32*t_pix/t_clock + t_endrow/t_clock
    pixel_arr = np.zeros(((row_clocks * np.shape(img)[1] + t_endframe/t_clock) * \
                          np.shape(img)[0] ) )#+500#.astype(int)

    for z in range(np.shape(img)[0]):
        if z % (np.shape(img)[0] /100.) < np.shape(img)[0]/100.: #update progress
            print str(int(round(float(z) / np.shape(img)[0] * 100.)))+ "% complete."
        for y in range(np.shape(img)[1]):
            for x in range(np.shape(img)[2]/32):
                pixel_arr[z*(row_clocks*np.shape(img)[1] + t_endframe/t_clock) + \
                          y*row_clocks +     x*t_pix/t_clock : \
                          z*(row_clocks*np.shape(img)[1] + t_endframe/t_clock) + \
                          y*row_clocks + (x+1)*t_pix/t_clock] \
                    = np.median(img[z,y,x*32:(x+1)*32])


    #fig = plt.figure()
    #plt.loglog(f, Pxx_den)

    if welch:
        print "Computing Welch FFT"
        bin = 7
        xaxis_w, powspec_w = signal.welch(pixel_arr, t_clock**-1, nperseg=len(pixel_arr)/bin)#, 
        #scaling='spectrum')
        #pdb.set_trace()
    #else:
    #bin = 1 #only applies to Welsh filter

    print "Computing FFT."
    #pad array with 0s to make fft faster
    #print "before length", len(pixel_arr)
    mylen = 2
    while mylen < len(pixel_arr):
        mylen *= 2
    pad = mylen-len(pixel_arr)
    pixel_arr = np.append(np.zeros(pad/2), pixel_arr)
    pixel_arr = np.append(pixel_arr, np.zeros(mylen - len(pixel_arr)))
    #print "after length", len(pixel_arr)
    
    fft = (np.fft.fft(pixel_arr)) /len(pixel_arr)
    xaxis = np.fft.fftfreq(pixel_arr.size, d=t_clock)
    powspec = (abs(fft))**2

    cumsum = np.cumsum(powspec)

    print "Plotting"
    plt.figure(1, figsize=(12, 8), dpi=100)
    plt.subplot(211)
    plt.plot(xaxis, powspec/2, 'b', label='Unfiltered power spectrum')
    plt.plot(xaxis_w, powspec_w, 'r', label="Welch filtered power spectrum") #welch filtered version
    plt.legend(loc=1, fontsize='small', framealpha=0.5)
    plt.xlim((np.shape(img)[0]*(row_clocks*np.shape(img)[1] + \
                                t_endframe/t_clock)*t_clock)**(-1), \
             0.5*t_pix**(-1))
    plt.ylim((1e-6, 1e1))
    plt.yscale('log')
    plt.xscale('log')
    plt.title('Power Spectrum '+filename + ' ' +str(np.shape(img)))
    plt.ylabel('Power [arb. units]')
    plt.xlabel('Frequency [Hz]')

    plt.subplot(212)
    plt.plot(xaxis, cumsum)
    plt.xlim((np.shape(img)[0]*(row_clocks*np.shape(img)[1] + \
                                t_endframe/t_clock)*t_clock)**(-1) , \
             0.5*t_pix**(-1))
    plt.title('Cumulative sum of power')
    plt.ylabel('Cumulative Sum [arb. units]')
    plt.xlabel('Frequency [Hz]')
    plt.yscale('log')
    plt.xscale('log')

    plt.tight_layout() #prevent labels from going off the edge
    plt.show()
