"""Strehl ratio calculator"""
#This software was written by Nem Jovanovic on the 18th of Feb 2013
#It can take a data frame and dark frames and calcualte the strehl ratio of the image as compared to an ideal PSF.
#Tips for using the software:
#1. You need to change the working directory "wder" to the location of the frames
#2. You need to change the name of the birght and drak files in the first two sections as required.
#3. If you dont have dark frames, then you can comment out the dark frame section and the "-dark" in finalNT on line 59.
#4. If there are cosmic rays that have not calibrated out the Gaussian fit should still work as long as there is not more then one pixel.
#If there are more then one pixel of high signal noise then you need to preselect a region of 200x200 pixels around the PSF manually.
#These values should be entered into section "Finding area of interest"
#5. You can adjust the threshold of the selection of the region of the data above the noise. This can be done in section "Thresholding noisy data"
#6. If you change wavelength you will need to alter the radius of the pupil to get the right plate scale. This is the rad variable in section "Computing ideal PSF" 

import numpy as np
import pyfits as py
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import scipy as sp
interactive(True)
import scipy as sp
from scipy import optimize
from scipy.optimize import curve_fit
from numpy.fft import fftshift as shift
from numpy.fft import ifftshift as ishift
from numpy.fft import fft2 
from scipy import misc 
import cmath as cm
from scipy import stats
import sys
import os
import pdb
import hot_pixels as hp
import PSF_shifter          #imports a program used to shift the PSF to the center of a pixel
import PSF_generator as PSFgen
import cmath as cm
from scipy import stats
import math as mt
from matplotlib.colors import LogNorm

home = os.getenv('HOME')
sys.path.append(home+'/src/lib/python/')
from image_registration import chi2_shifts as CC

rcParams['image.origin'] = 'lower'                  #Flips data about the horizontal axis 
rcParams['image.interpolation'] = 'nearest'         # does not interpolate pixels in imshow



"""Finding centre and height of data"""
"""2D Gaussian curve fit"""
def Gaus(Parameters, v1, v2):
    return Parameters[0]*np.exp(-((v1-Parameters[1])**2/(2*Parameters[3]**2))-((v2-Parameters[2])**2/(2*Parameters[4]**2)))+Parameters[5]

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = np.sum(data)
    X, Y = np.indices(np.shape(data)) 
    cx = (np.sum(Y*data)/total)
    cy = (np.sum(X*data)/total) 
    row = data[round(cx), :]
    width_x = (np.sqrt(np.sum(abs((np.arange(row.size)-cx)**2*row))/np.sum(row)))
    col = data[:, round(cy)]
    width_y = (np.sqrt(np.sum(abs((np.arange(col.size)-cy)**2*col))/np.sum(col)))
    DCoffset=np.median(data[(dist1>irad)*(dist1<orad)])
    height = np.max(data)-DCoffset
    return height, cx, cy, width_x, width_y, DCoffset

"""Fit Gaussian function"""
def fitGaussian(data, v1, v2):
    Parameters0=moments(data)
    errorfunction=lambda p, v1, v2: np.ravel(Gaus(p, v1, v2)-data)
    p, success = sp.optimize.leastsq(errorfunction, Parameters0, args=(v1, v2),maxfev=10000)
    return p

def sqrt(im):
    scale_min = np.min(im)
    scale_max = np.max(im)
    im.clip(scale_min, scale_max)
    im = im - scale_min
    indices = np.where(im < 0)
    im[indices] = 0.0
    im = np.sqrt(im)
    im = im / mt.sqrt(scale_max - scale_min)
    return im

def center(im):
    Linemean1 = np.mean(im,axis=0)
    Linemean2 = np.mean(im,axis=1)
    ysize, xsize=np.shape(im)
    x = np.linspace(0,xsize-1,xsize)
    y = np.linspace(0,ysize-1,ysize)
    xcoor = np.argmax(Linemean1)                
    ycoor = np.argmax(Linemean2)
    return ycoor, xcoor



###--------------------Importing ideal PSF----------------------###
PSFref  = py.getdata('./Support_files/New_PSF_SCExAO_for_chuck_in_Hband.fits')#Imported
ny,nx   = np.shape(PSFref)
ycoor1  = np.argmax(PSFref)/nx                        #finds y coordiante of peak flux
xcoor1  = np.argmax(PSFref)%nx                        #finds x coordiante of peak flux
PSFmini = PSFref[ycoor1-100:ycoor1+100,xcoor1-100:xcoor1+100]


###--------------Encircled flux for ideal PSF------------------###
Y1, X1   = np.indices(np.shape(PSFmini))
dist1    = np.hypot(Y1-100, X1-100)   
rad      = np.linspace(0,35,36)
rads     = np.size(rad)
EncFlux1 = np.zeros_like(rad)
pixnum1  = np.zeros_like(rad) 
Window   = np.zeros_like(PSFmini)
radial   = np.zeros_like(rad)
for i in range(0,rads):
    Window[dist1<=(rad[i])*2] = 1   
    EncFlux1[i] = np.sum(Window*PSFmini)
    pixnum1[i] = np.sum(Window)
    Window = np.zeros_like(PSFmini)

coefscomp = np.lib.polyfit(pixnum1[22:],EncFlux1[22:],1)
fit_y_comp = np.lib.polyval(coefscomp,pixnum1[22:])


###---------------PSF ratio of peak to total count------------###
irad=80
orad=90


#ParametersComp = fitGaussian(PSFmini, X1, Y1)
PeakComp = np.max(PSFmini)#ParametersComp[0]#
FluxComp = coefscomp[1]
RatioComp = PeakComp/FluxComp     


###----------------importing data-----------------------###
wder1='/media/data/20170531/ircam1log/'
data=py.getdata(wder1+'ircam1_13:00:05.297715089.fits')    #Change file name as required 
zsize,ysize,xsize=np.shape(data)                            #Size of frame


###---------------------Importing darks----------------------###
wder2='/media/data/20170531/ircam1log/darks/150043/'
dark=py.getdata(wder2+'dark_0000100.fits')
#dark=np.mean(im,axis=0)                                    #Taking the mean of the darks
zsize, ysize, xsize = np.shape(data)


###------------------Dark and Hot pixel removal---------------###
finalb=data-dark  #Dark subtraction

for i in range(zsize):
    useless,finalb[i]=hp.find_outlier_pixels(finalb[i])
#final[0]=np.mean(final,axis=0)


###-------------------------Binning---------------------------###

binsiz=np.array([1000])# ([1000,500,200,100,50,20,10,5,4,3,2,1])
b=np.size(binsiz)
AveS=np.zeros(b)

for j in range(b):
    final=np.zeros([zsize/binsiz[j],256,320])
    for i in range(0, zsize/binsiz[j]):
        final[i]=np.mean(finalb[i*binsiz[j]:i*binsiz[j]+binsiz[j]], axis=0)

    maxpts = zsize/binsiz[j]
    Strehl = np.zeros((maxpts))

    finalmini=np.zeros((200,200))
    finalmininorm=np.zeros_like((200,200))
    data_s=np.zeros_like((200,200))
    data_sn=np.zeros_like((200,200))


    ###--------------------Computing ideal PSF------------------###
        #This code generates the Subaru pupil and then calculates the Fourier transform
        #to get the PSF in the image plane. It is possible to also add a wavefront to
        #the pupil map. 
        #The code can fit a Gaussian to the data to check the width and centre position.     
    '''
    #PSF=PSFgen.PSF_gen(0.315)                       #Intensity profile of the PSF/calcualted
    PSFref=py.getdata('Reference_PSF_for_Strehl.fits') #Imported
    ny,nx=np.shape(PSFref)
    ycoor1=np.argmax(PSFref)/nx                        #finds y coordiante of peak flux
    xcoor1=np.argmax(PSFref)%nx                        #finds x coordiante of peak flux
    PixPeakComp=np.max(Int)                        
    #Intnorm=np.divide(Int,PixPeakComp)                 #Normalizes the PSF
    PSFmini=PSFref[ycoor1-100:ycoor1+100,xcoor1-100:xcoor1+100]
    '''

    ###---------------Encircled flux for calculated PSF----------###
    '''    
    Y1, X1 =  np.indices(np.shape(PSFmini))
    dist1 = np.hypot(Y1-100, X1-100)   
    rad = np.linspace(0,35,36)
    rads = np.size(rad)
    EncFlux1 = np.zeros_like(rad)
    pixnum1 = np.zeros_like(rad) 
    Window = np.zeros_like(PSFmini)
    radial = np.zeros_like(rad)
    for i in range(0,rads):
        Window[dist1<=(rad[i])*2] = 1
        EncFlux1[i] = np.sum(Window*PSFmini)
        pixnum1[i] = np.sum(Window)
        Window = np.zeros_like(PSFmini)
        #if i==0:    #calcualting radial profile of PSF
        #    radial[i]=(EncFlux1[i])/(pixnum1[i])
        #else:
        #    radial[i]=(EncFlux1[i]-EncFlux1[i-1])/(pixnum1[i]-pixnum1[i-1])

    #figure()
    #plot(radial)
    
    coefscomp = np.lib.polyfit(pixnum1[22:],EncFlux1[22:],1)
    fit_y_comp = np.lib.polyval(coefscomp,pixnum1[22:])
    #### Plot encircled flux vs # of encircled pixels ####
    #figure()
    #plt.scatter(pixnum1,EncFlux1)
    #plt.plot(pixnum1[22:],fit_y_comp, 'r')
    '''

    ###------------------PSF ratio of peak to total count---------------###
    '''    
    ParametersComp = fitGaussian(PSFmini, X1, Y1)
    PeakComp = ParametersComp[0]
    FluxComp = coefscomp[1]
    RatioComp = PeakComp/FluxComp     #Calculates the Peak to the total encircled flux for the
    '''

    for k in range(0, maxpts):
                             
        ###------------------Centre of gravity calculator-------------------###
        ycoor2, xcoor2 = center(final[k])
        
        finalmini=final[k,ycoor2-100:ycoor2+100,xcoor2-100:xcoor2+100]
        finalmini=finalmini-np.median(finalmini)

        PixPeakData=np.max(finalmini)
        finalmininorm=finalmini/float(PixPeakData)

        Y2, X2 =  np.indices(np.shape(finalmininorm))
        dist2 = np.hypot(Y2-100, X2-100)    #Creates a distance function from the centre
        #figure()
        #imshow(finalNTmininorm)

        ParametersData=fitGaussian(finalmininorm, X2, Y2)
        #print ParametersData

        deltax = 100-ParametersData[1]     #finds distance from the center of the box in x 
        deltay = ParametersData[2]-100     #finds distance from the center of the box in y
        
        #figure()
        #imshow(finalmininorm)
        data_s = PSF_shifter.PSF_shift(finalmininorm, deltax, deltay) 
        data_s = abs(data_s)
        #figure()
        #imshow(data_s)
        data_sn = data_s/np.max(data_s)
        #ParamData_s=fitGaussian(data_sn, X2, Y2)
        #print ParamData_s

        """Computing ideal PSF"""
        #This code generates the Subaru pupil and then calculates the Fourier transform
        #to get the PSF in the image plane. It is possible to also add a wavefront to
        #the pupil map. 
        #The code can fit a Gaussian to the data to check the width and centre position. 
        
        #PSF=PSFgen.PSF_gen(0.315)
        #Int=PSF[0]                                      #Intensity profile of the PSF
        #ny,nx=np.shape(PSF[0])
        #ycoor2=np.argmax(Int)/nx                        #finds y coordiante of peak flux
        #xcoor2=np.argmax(Int)%nx                        #finds x coordiante of peak flux
        #PixPeakComp=np.max(Int)                        
        #Intnorm=np.divide(Int,PixPeakComp)              #Normalizes the PSF
        #Intmininorm=Intnorm[ycoor2-100:ycoor2+100,xcoor2-100:xcoor2+100]
        
        """Normalise/resize"""
        #figure()
        #imshow(data_sn[70:130,70:130])
        #plt.savefig("Open_loop_PSFaveraged_Measured_PSF.pdf", dpi=300, facecolor='gray')
        #figure()
        #imshow(Intmininorm[70:130,70:130])
        #plt.savefig("Perfect_PSF_Computed_PSF.pdf", dpi=300, facecolor='gray')

        #diff=Intmininorm-data_sn
        #total = np.sum(diff[94:106,94:106]**2)
        #print total
        #fig = figure()
        #ax = fig.add_subplot(111)
        #cax = ax.imshow(diff)
        #cbar = fig.colorbar(cax, shrink = 0.8)
        #cbar.ax.set_yticklabels(['0.0','0.25','0.5','0.75','1.0'])
        #cbar.set_clim(-0.1, 0.3)
        #plt.savefig("Diff_PSF.pdf", dpi=300, facecolor='gray')

        """Encircled energy"""
        #dist3 = np.hypot(Y1-100, X1-100)  
        #Y2, X2 =  np.indices(np.shape(Intmininorm))              
        #FluxComp=np.sum(Intmininorm)#[dist4<iradDCO])             
        #Calcualtes the Flux within the predefined aperture for the computed frame

        """Photometric background extraction"""
        rad = np.linspace(0,35,36)
        rads = np.size(rad)
        EncFlux2 = np.zeros_like(rad)
        pixnum2 = np.zeros_like(rad)
        #radial = np.zeros_like(rad) 
        Window = np.zeros_like(data_sn)

        for i in range(0,rads):
            Window[dist2<=(rad[i])*2] = 1
            EncFlux2[i] = np.sum(Window*data_sn)
            pixnum2[i] = np.sum(Window)
            Window = np.zeros_like(data_sn)
            #if i==0:               #for creating a radial profile
            #    radial[i]=(EncFlux2[i])/(pixnum2[i])
            #else:
            #    radial[i]=(EncFlux2[i]-EncFlux2[i-1])/(pixnum2[i]-pixnum2[i-1])
        #figure()
        #plt.scatter(pixnum,EncFlux)
        #slope, intercept, r_value, p_value, std_err = stats.linregress(pixnum[22:],EncFlux[22:])
        coefs = np.lib.polyfit(pixnum2[22:],EncFlux2[22:],1)
        fit_y = np.lib.polyval(coefs,pixnum2[22:])
        #plt.plot(pixnum[22:],fit_y, 'g')

        #EncFlux2 = np.zeros_like(rad)
        #pixnum2 = np.zeros_like(rad) 

        #for i in range(0,rads):
        #    Window[dist3<=(rad[i]*2)] = 1
        #    EncFlux2[i] = np.sum(Window*Intmininorm)
        #    pixnum2[i] = np.sum(Window)
        #    Window = np.zeros_like(Intmininorm)

        #plt.scatter(pixnum2,EncFlux2)

        #coefscomp = np.lib.polyfit(pixnum2[22:],EncFlux2[22:],1)
        #fit_y_comp = np.lib.polyval(coefscomp,pixnum2[22:])
        #plt.plot(pixnum2[22:],fit_y_comp, 'r')
        #xlabel('Number of pixels')
        #xlim(-100,20000)
        #ylabel('Encircled count')
        #ylim(0)
        #plt.savefig("Open_loop_PSFaveraged_Encircled_flux_vs_number_of_pixels.pdf", dpi=300, facecolor='gray')


        """Strehl Ratio Calcualtion"""
        PeakData = np.max(data_sn)#ParamData_s[0]#
        #ParametersComp = fitGaussian(Intmininorm, X2, Y2)
        #print ParametersComp
        #PeakComp = ParametersComp[0]
        FluxData = coefs[1]
        #FluxComp = coefscomp[1]
        RatioData = PeakData/FluxData     
        #Calculates the Peak to the total encircled flux for the data
        #RatioComp = PeakComp/FluxComp     
        #Calculates the Peak to the total encircled flux for the computed frame
        Strehl[k] = RatioData/RatioComp   
        #Calcualtes the Strehl ratio by dividing the two ratios above
        AveS[j] = np.mean(Strehl)

    ###------------------------Writing to fits-----------------------------###
    #hdu = py.PrimaryHDU(radial)
    #hdu.writeto('PSF_radial_profile.fits')  #   

    #figure()
    #semilogy(rad[0:20], radial[0:20],label='Laboratory PSF')
    #plt.xlabel('Radius (pixels)')
    #plt.ylabel('Azimuthally averaged flux')
    #plt.legend()

    #print Strehl')
    #x=np.linspace(0,maxpts+1,maxpts)
    #figure()
    #plt.scatter(x, Strehl)
    #plt.xlim(0,maxpts)
    #plt.ylim(-0.05,1.05)
    #plt.xlabel('Data point number')
    #plt.ylabel('Strehl ratio')
    #plt.savefig("Strehl_ratio_vs_data_pt_no_UT15_41_32.pdf", dpi=300, facecolor='gray')

    #AveS[j] = np.mean(Strehl)
    #print AveS#np.max(Strehl)#AveS
    #print Strehl[0]

###------------------------Plotting data------------------------###
'''
figure()
semilogx(binsiz, AveS)
plt.xlabel('Bin size (# of frames)')
plt.ylabel('Average Strehl along cube')
plt.xlim(10**0,10**3)
plt.ylim(0.0,1.0)
plt.scatter(binsiz, AveS)
#plt.savefig("Strehl_ratio_vs_binsize_HD85672_set_07_58_46.pdf", dpi=300, facecolor='gray')
#plt.savefig("Strehl_ratio_vs_binsize_HD85672_set_07_58_46.png", dpi=300, facecolor='gray')
'''

figure()
a=np.linspace(0,maxpts-1,maxpts)
plt.scatter(a,Strehl)
plt.xlabel('Data pt number')
plt.ylabel('Instantaneous Strehl')
plt.xlim(-10,1010)
plt.ylim(0,1.0)
#plt.savefig("Strehl_ratio_vs_data_points_HD85672_set_07_58_46.pdf", dpi=300, facecolor='gray')
#plt.savefig("Strehl_ratio_vs_data_points_HD85672_set_07_58_46.png", dpi=300, facecolor='gray')

print AveS
print np.std(Strehl)

figure()
clipped=np.clip(finalmininorm,2e-4,1)
imshow(clipped,cmap='CMRmap',norm=LogNorm(vmin=1e-4,vmax=1))
cbar = plt.colorbar()
plt.show()
#plt.savefig("New_PSF_lab_turb_0225_lc03_02_13_53_log.png", dpi=300, bbox_inches='tight')
figure()
imshow(finalmininorm,cmap='CMRmap',vmin=0, vmax=1)
cbar = plt.colorbar()
plt.show()
#plt.savefig("New_PSF_lab_turb_0225_lc03_02_13_53_lin.png", dpi=300, bbox_inches='tight')

#First run with PyWFS
#Low order and  PyWFS on 28th October - ircam1_09:22:30.998536713.fits
#Low orders removed on 28th October - ircam1_09:24:23.207711953.fits
#AO188 only October 28th - ircam1_09:26:27.501809485.fits
#Darks for the data bias0000100.fits

#Second run with PyWFS
#PyWFS closed on 28th October second time - ircam1_10:49:22.465200051.fits
#AO188 only on 28th October second time - ircam1_10:51:08.687436533.fits

#29th of October - good PSF on HD15115
#PyWFS on 29th - ircam1_10:21:32.951022133.fits

#30th PSF calibration in the lab - ircam1_03:41:41.113198311.fits

#30th PSF on HR8799 looks good - ircam1_08:01:24.351939354.fits

#30th before FPWS - ircam1_04:40:38.067233310.fits 
#30th after closing FPWS - ircam1_04:41:25.438950075.fits
#500us int. dark - dark_0000500.fits Altair

#19th of March 2016 - good set with PyWFS closed after 3 RMs. ircam1_12:57:27.209495394.fits

#20th really good Strehl at the end of the night

#No modulation sci_13:51:39.475833_20000.0us.fits
#With Turbulance before image 815 sci_13:46:56.008727_20000.0us.fits
#With AO after image 815 sci_13:46:56.008727_20000.0us.fits
#sci_13:47:17.443385_20000.0us.fits for Strehl ratio measurement from November
