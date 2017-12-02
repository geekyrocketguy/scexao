import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb

#Computes the vibrations in images using a centroid on a non-saturated PSF.

framerate = 5848 #frames per second
print "Frame rate is assumed to be", framerate, "fps."

def main(files):

    if type(files) == str:
        files = [files]

    n_brightest = 17 #number of pixels to include in centroid calculation

    for i in range(len(files)):

        print "Working on file", i+1, "of", np.size(files), "."
        img = pyfits.getdata(files[i])

        if i==0:
            shifts = np.zeros((2, len(img)))
        else:
            shifts = np.append(shifts, np.zeros((2, len(img))), axis=1)

        for z in range(len(img)):

            if z%1000 == 0:
                print int(round(float(z)/len(img)/np.size(files)*100.)), "% done."

            im = np.copy(img[z])

            threshold = np.sort(im.flatten())[np.size(im)-n_brightest] #29th brightest pixel
            im -= threshold
            im[im<0] = 0

            x1 = np.sum(np.sum(im, 0)*range(np.shape(im)[1])) / np.sum(im)
            y1 = np.sum(np.sum(im, 1)*range(np.shape(im)[0])) / np.sum(im)

            shifts[0, shifts.size/2 - len(img) + z] = x1
            shifts[1, shifts.size/2 - len(img) + z] = y1

    shifts[0,:] -= np.median(shifts[0,:]) #subtract average X value
    shifts[1,:] -= np.median(shifts[1,:]) #subtract average Y value

    print "Saving shifts."
    np.savetxt(files[0].replace('.fits', '_shifts.csv'), shifts, delimiter=",")

    time_axis = np.arange(shifts.size / 2.) / framerate

    plt.figure(1, figsize=(10,10), dpi=100)
    plt.subplot(211)
    plt.plot(time_axis, shifts[0,:])
    plt.title('X Shift')
    plt.xlabel('time (s)')
    plt.ylabel('Pixels')

    plt.subplot(212)
    plt.plot(time_axis, shifts[1,:])
    plt.title('Y Shift')
    plt.xlabel('time (s)')
    plt.ylabel('Pixels')

    plt.show()

    
def compare_psds(off_file, on_file):
    off = np.loadtxt(off_file, delimiter=',')    
    on = np.loadtxt(on_file, delimiter=',')

    x_xaxis_off = np.fft.fftfreq(off[0,:].size, d=1./framerate)
    x_fft_off = (np.fft.fft(off[0,:])) /len(off[0,:])
    x_powspec_off = (abs(x_fft_off))**2

    x_xaxis_on = np.fft.fftfreq(on[0,:].size, d=1./framerate)
    x_fft_on = (np.fft.fft(on[0,:])) /len(on[0,:])
    x_powspec_on = (abs(x_fft_on))**2

    y_xaxis_off = np.fft.fftfreq(off[1,:].size, d=1./framerate)
    y_fft_off = (np.fft.fft(off[1,:])) /len(off[1,:])
    y_powspec_off = (abs(y_fft_off))**2

    y_xaxis_on = np.fft.fftfreq(on[1,:].size, d=1./framerate)
    y_fft_on = (np.fft.fft(on[1,:])) /len(on[1,:])
    y_powspec_on = (abs(y_fft_on))**2

    plt.figure(1, figsize=(10,10), dpi=100)

    plt.subplot(221)
    plt.plot(x_xaxis_on, x_powspec_on, label="SC On")
    plt.plot(x_xaxis_off, x_powspec_off, label="SC Off")
    plt.yscale('log')
    plt.xscale('log')
    plt.title("X axis")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.xlim((1,3e3))
    plt.legend()

    plt.subplot(222)
    plt.plot(y_xaxis_on, y_powspec_on, label="SC On")
    plt.plot(y_xaxis_off, y_powspec_off, label="SC Off")
    plt.yscale('log')
    plt.xscale('log')
    plt.title('Y axis')
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.xlim((1,3e3))
    plt.legend()

    plt.subplot(223)
    plt.plot(x_xaxis_on, np.cumsum(x_powspec_on) - np.mean(np.cumsum(x_powspec_on[:100])), \
             label="SC On")
    plt.plot(x_xaxis_off, np.cumsum(x_powspec_off) - np.mean(np.cumsum(x_powspec_off[:100])), \
             label="SC Off")
    #plt.yscale('log')
    plt.xscale('log')
    plt.title("X axis")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.xlim((1,3e3))
    plt.legend()

    plt.subplot(224)
    plt.plot(y_xaxis_on, np.cumsum(y_powspec_on) - np.mean(np.cumsum(y_powspec_on[:100])), \
             label="SC On")
    plt.plot(y_xaxis_off, np.cumsum(y_powspec_off) - np.mean(np.cumsum(y_powspec_off[:100])), \
             label="SC Off")
    #plt.yscale('log')
    plt.xscale('log')
    plt.title('Y axis')
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.xlim((1,3e3))
    plt.legend()

    plt.show()
