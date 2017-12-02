import pyfits
import numpy as np
import matplotlib.pyplot as plt
import os
import pdb
from shutil import copy2

#Define txt file that contains the list of images to process
#filelist = 'reduction_olivid.txt'

def main(filelist):

    #GET LIST OF FILES TO PROCESS
    file = open(filelist, 'r')
    params = np.array(file.readlines())

    dir = params[np.squeeze(np.where(params == 'dir\n'))+1].replace('\n', '')
    dataloc = np.squeeze(np.where(params == 'data\n'))
    darkloc = np.squeeze(np.where(params == 'dark\n'))

    darkfiles = np.array([])
    datafiles = np.array([])

    for i in range(len(params)):
        line = params[i]
        if (line[0] != '#') & (line[0] != '\n'): #if not a comment or blank line
            line = line.replace('\n', '') #stupid new line character
            if (i>dataloc) & (i<darkloc):
                datafiles = np.append(datafiles, line)
            elif (i>darkloc):
                darkfiles = np.append(darkfiles, line)
    file.close()


    #PROCESS FILES
    dark = get_dark(dir, darkfiles, recreate='T')

    subtract(dir, datafiles, dark)

  

def get_dark(dir, darkfiles, recreate='F'):
    filename = 'saphcamdarks/'+dir[12:20]+'_'+darkfiles[0][:-15]+'_darkmaster.fits'
    if (os.path.isfile(filename)) & (recreate=='F'): #if the file already exists
        print("Restoring dark frame.")
        return pyfits.getdata(filename)
    else: #create file
        print("Creating dark frame.")
        for i in range(len(darkfiles)):
            print(" Reading image "+ str(i+1) + " of " + str(len(darkfiles))+'.')
            params = darkfiles[i].split()
            img = pyfits.getdata(dir+params[0])
            if len(params) == 3: 
                x = np.arange(len(img))
                img = img[np.where((x < int(params[1])) | (x > int(params[2])))]
            elif len(params) != 1:
                print "Wrong number of parameters entered in file."
                pdb.set_trace()

            try:
                darkcube = np.append(darkcube, img, axis=0)
            except:
                darkcube = img #first iteration

        #Stupid Overflow correction
        if np.min(darkcube) < 0:
            print(" Correcting overflow...")
            darkcube = darkcube + (darkcube < 0)*2**16
            print(" Done.")

        print(" Darkcube shape:" + str(np.shape(darkcube)))
        print(" Medianing dark frame.")
        darkmaster = np.median(darkcube, axis=0)
        print(" Done.")

        if 1:
            plt.imshow(darkmaster, interpolation='none')
            plt.colorbar()
            plt.title(filename)
            plt.show()

        if not os.path.isdir('saphcamdarks'):  #Check if directory exists
            os.makedirs('saphcamdarks') #make directory
            print(" Created directory "+'saphcamdarks')

        print(" Writing fits file: " + filename)
        pyfits.writeto(filename, darkmaster, clobber='T')
        return darkmaster


def subtract(dir, datafiles, dark):
#Subtracts dark from data frames and saves the results

    for i in range(len(datafiles)):
        print('')
        print('Processing file '+str(i+1)+" of "+str(len(datafiles))+'.')
        params = datafiles[i].split()
        img = pyfits.getdata(dir+params[0])

        #Stupid Overflow correction
        if np.min(img) < 0:
            print(" Correcting overflow...")
            img = img + (img < 0)*2**16
            print(" Done.")
    
        if not os.path.isdir(dir+'processed'):  #Check if directory exists
            os.makedirs(dir+'processed') #make directory
            print(" Created directory "+dir+'processed')


        print(" Performing dark frame subtractions.")
        if len(params) == 1: #process all of file
            for j in range(np.shape(img)[0]):
                img[j,:,:] -= dark
                copy2(dir+params[0][:-5]+'.txt', dir+'processed/'+params[0][:-5]+'_p.txt')

        elif len(params) == 3: #file has bad frames that should be skipped
            file = open(dir+params[0][:-5]+'.txt', 'r')
            timestamps = np.array(file.readlines())
            file.close()

            #remove bad frames and associated timestamps
            x = np.arange(len(img))
            x = x[np.where((x < int(params[1])) | (x > int(params[2])))]
            img = img[x]
            timestamps = timestamps[x]

            #subtract dark frame
            for j in range(np.shape(img)[0]):
                img[j,:,:] -= dark
    
            #save timestamps
            file = open(dir+'processed/'+params[0][:-5]+'_p.txt', 'w')
            file.write(timestamps)
            file.close()

        else: #user screwed up
            print(" Wrong number of parameters entered in file.")
            print(" params:")
            pdb.set_trace()

        #save image
        filename = dir+'processed/'+params[0][:-5]+'_p.fits'
        print(" Writing fits file: " + filename)
        pyfits.writeto(filename, img, clobber='T')
        
