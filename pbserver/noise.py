#computes noise for different readout modes

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb

def utr(file, f1, f2):
    if len(locals()) != 3:
        print('Syntax: noise.utr(file, frame_start, frame_end)')
        return
    
    print("Reading image...")
    img = pyfits.getdata(file).astype(float)
    print("Image read. Computing average frame.")
    
    img = img[f1:f2] #crop image

    bias = np.median(img, 0)
    print("Done.")
    
    noise_1frame = np.array([])
    noise_2frame = np.array([])

    for i in range(len(img)):
        #print("Calculating noise in image " + str(i) + " of " + str(len(img)) + ".")
        noise_1frame = np.append(noise_1frame, np.std(img[i]-bias, ddof=1))

        if (i%2==0) & (i>0):
            noise_2frame = np.append(noise_2frame, np.std(img[i]-img[i-1], ddof=1))

            #plt.imshow(img[i]-bias, interpolation='none', vmin=-50, vmax=50)
            #plt.colorbar()
            #plt.show()

    print("Noise (frame[i] - median frame): "+ str(np.median(noise_1frame)))
    print("Noise (frame[i] - frame[i-1]): "+ str(np.median(noise_2frame)))

    print()


def linearity(file):
    if len(locals()) != 1:
        print("Syntax: noise.linearity('file')")
        return
    
    print("Reading image...")
    img = pyfits.getdata(file).astype(float)
    print("Image read. Computing average frame.")
    
    bias = np.median(img, 0)
    print("Done.")
    
    meds = np.zeros(len(img))
    for i in range(len(img)):
        meds[i] = np.median(img[i]-bias)

    plt.plot(meds, 'o')
    plt.xlabel("Frame Number")
    plt.ylabel("ADUs")
    plt.show()


def rrr(file, f1, f2):
    if len(locals()) != 3:
        print("Syntax: noise.rrr('file', f1, f2)")
        return
    
    print("Reading image...")
    img = pyfits.getdata(file).astype(float)
    
    img = img[f1:f2] #crop image

    stddevs = np.zeros(len(img)-1)
    ydim = np.shape(img)[1]

    for i in range(len(img)-1):
        stddevs[i] = np.std(img[i, ydim/2: , 96:] - img[i+1, :ydim/2 , 96:], ddof=1)

    print("Noise: "+str(np.median(stddevs)))


def rr(file, f1, f2, f3, f4):
    if len(locals()) != 5:
        print("Syntax: noise.rrr('file', f1, f2, f3, f4)")
        print("f1, f2 denote bias frames. f3, f4 denote data frames.")
        return
    
    print("Reading image...")
    img = pyfits.getdata(file).astype(float)
    
    bias = np.median(img[f1:f2], 0)
    img = img[f3:f4] #crop image

    stddevs = np.zeros(len(img))
 
    for i in range(len(img)):
        stddevs[i] = np.std(img[i] - bias, ddof=1)

    print(stddevs)
    print("Noise: "+str(np.median(stddevs)))
    
