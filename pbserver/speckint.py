#For speckle interferometry


import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pdb

bin = 1 #integer number of frames that should be averaged before computing FFT

def makefft(filelist):
    #GET LIST OF FILES TO PROCESS
    file = open(filelist, 'r')
    params = np.array(file.readlines())

    dir = params[np.squeeze(np.where(params == 'dir\n'))+1].replace('\n', '')
    dataloc = np.squeeze(np.where(params == 'data\n'))

    datafiles = np.array([])

    for i in range(len(params)):
        line = params[i]
        if (line[0] != '#') & (line[0] != '\n'): #if not a comment or blank line
            line = line.replace('\n', '') #stupid new line character
            if (i>dataloc):
                datafiles = np.append(datafiles, line)
    file.close()

    
    for file in datafiles:
        img = pyfits.getdata(dir+file)
        frame = np.zeros((512,512)) #pad with 0s for better fft
        xdim = np.shape(img)[2]
        ydim = np.shape(img)[1]

        for i in range(len(file)/bin):

            #if bin != 1:
            lilimage = np.mean(img[i*bin : (i+1)*bin , : , :], 0)
            #else: 
            #    lilimage = img[i , : , :]

            lilimage -= np.median(lilimage) #background subtract
            frame[256-ydim/2 : 256+ydim/2 , 
                  256-xdim/2 : 256+xdim/2 ] = lilimage #embed image in 0s

            try:
                fftcube = np.append(fftcube, np.fft.fftshift(np.fft.fft2(frame)), axis=0)
            except:
                fftcube = np.fft.fftshift(np.fft.fft2(frame))

    fftcube = np.abs(fftcube**2) #make real
    fftcube = np.median(fftcube, 0)

    #save file
    pyfits.writeto(dir + datafiles[0][:16]+'_fft.fits', clobber='T')
