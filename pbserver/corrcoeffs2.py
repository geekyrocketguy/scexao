import matplotlib.pyplot as plt
import numpy as np
import pyfits
import pdb
import glob
from matplotlib.colors import LogNorm
#import pickle
import random

#like corrcoeffs, but doesn't read all the images into a cube. This uses less ram.

def main(restore=False, dataset='dataset'):

    #check if user specified which dataset to use
    dataset = dataset.lower()
    if (dataset != 'internal') & (dataset != 'onsky'):
        print "Please specify dataset='internal' or dataset='onsky'"
        return
        
    if dataset=='internal': 
        suffix='_internal.npy'
        plttitle="Pearson's Correlation Coefficients for Internal Source"
        datapath='/media/data/20180413/saphira/processed/pbimage*fits'
    else: 
        suffix='_onsky.npy'
        plttitle="Pearson's Correlation Coefficients for On-Sky Data"
        datapath='/media/data/20170911/saphira/processed/pbimage*fits'

    if restore != True:
        print "Calculating data, not restoring it."

        filelist = sorted(glob.glob(datapath))

        print "Identifying annulus."

        n_files_avg = 5 #number of files to average when finding the psf core

        for i in np.linspace(0, len(filelist)-1, n_files_avg):
            if i==0:
                cube = pyfits.getdata(filelist[int(i)])
            else:
                cube += pyfits.getdata(filelist[int(round(i))])

        cube /= n_files_avg
        avgcube = np.mean(cube, 0)
        im = np.copy(avgcube)

        n_brightest = 12
        threshold = np.sort(avgcube.flatten())[np.size(avgcube)-n_brightest] #29th brightest pixel
        avgcube -= threshold
        avgcube[avgcube<0] = 0

        x0 = np.sum(np.sum(avgcube, 0)*range(np.shape(avgcube)[1])) / \
            np.sum(avgcube)
        y0 = np.sum(np.sum(avgcube, 1)*range(np.shape(avgcube)[0])) / \
            np.sum(avgcube)

        #0.0107 is arcseconds per pixel
        lambda_D_pix = 1.63e-6 / 8.2 * 206265. / .0107 #lambda/D in pixels
        rmin = 3 * lambda_D_pix
        rmax = 10 * lambda_D_pix

        #make annulus
        annulus = np.zeros(np.shape(im))
        for y in range(np.shape(annulus)[0]):
            for x in range(np.shape(annulus)[1]):
                r2 = (x-x0)**2 + (y-y0)**2
                if (r2 > rmin**2) & (r2 < rmax**2):
                    annulus[y,x] = 1

        #check region of interest
        if 0:
            plt.imshow(im, interpolation='none', vmin=50, vmax=5e4, norm=LogNorm())
            plt.plot([x0], [y0], 'kx', markersize=15)

            theta = np.linspace(0, 2*np.pi, 100)
            plt.plot(x0 + rmin * np.cos(theta), y0 + rmin * np.sin(theta))
            plt.plot(x0 + rmax * np.cos(theta), y0 + rmax * np.sin(theta))

            plt.colorbar()
            plt.tight_layout()
            plt.show()

            plt.imshow(annulus, interpolation='none')
            plt.show()

        print "Done."
        print "Reading in data."
        loc = np.where(annulus == 1)
        size_selection = np.shape(loc)[1] #number of pixels in annulus

        separations = np.append(range(1, 199), range(200, 1680*5, 8))
        separations = np.append(separations, range(1680*6, 1680*360, 1680))
        #save separations
        np.save('corr_coeff_separations'+suffix, separations)

        #np.append is really slow, so let's create a huge cube and then
        # delete the extra space at the end.
        #Assume no more than 1e4 frames per cube
        img2 = np.int32(np.zeros((len(filelist)*1e4 , size_selection)))
        stds = np.zeros(len(filelist)*1e4)

        z_index=0 #how many frames have I processed?
        for i in range(len(filelist)):
            img = pyfits.getdata(filelist[i])

            print ' '+ str(int(round(float(i)/len(filelist)*100.))) + '% complete'

            for z in range(len(img)):
                img2[z_index] = img[z][loc]
                img2[z_index] -= np.mean(img2[z_index])
                stds[z_index] = np.std(img2[z_index])

                z_index += 1

        #delete excess space
        img2 = img2[:z_index]
        stds = stds[:z_index]
        print "img2 type:", type(img2[30,30])
        framenos = range(0, len(img2)-100, 1)
        random.shuffle(framenos) #so it samples randomly in time

        print "Analyzing data now."
        matrix = np.float32(np.zeros((len(framenos), len(separations)))) #i, sep
        matrix -= 99 #makes it easier to identify unfilled areas

        for no in range(len(framenos)):
            i = framenos[no]
    
            if no%10 == 0:
                print ' '+ str(round(float(no)/len(framenos)*100.*1000.)/1000.)+ "% done."

            for s in range(len(separations)):
                sep = separations[s]

                if i+sep >= len(img2): continue #skip to next iteration

                #clear definition
                #rho = np.sum((sel1 - np.mean(sel1)) * (sel2 - np.mean(sel2)) ) / \
                #             (np.std(sel1) * np.std(sel2) * np.size(sel1))

                #efficient definition
                rho = np.sum(img2[i]  * img2[i+sep]) / \
                             (stds[i] * stds[i+sep] * size_selection)

                matrix[no, s] = rho
    
            #save matrix
            if (no % 10000 == 0) & (no > 0):
                print "Saving."
                print "matrix type:", type(matrix[30,30])
                #default is 64-bit. Let's cut the size in half.
                np.save('corr_coeff_matrix'+suffix, matrix)
                print "Saved."
                print

    else: #restore data from text file
        print "Restoring data from file."

        matrix = np.load('corr_coeff_matrix'+suffix)
        separations = np.load('corr_coeff_separations'+suffix)
        print "Done."

    if 0:
        plt.imshow(matrix, interpolation='none')
        plt.colorbar()
        plt.xlabel('separation')
        plt.ylabel('i')
        plt.show()

    print "Removing empty data, flattening matrix."
    flattened = np.zeros(len(separations))
    for i in range(len(flattened)):
        flattened[i] = np.mean((matrix[:,i])[np.where(matrix[:,i] != -99)])
    print "Done."

    plt.plot(separations/1680., flattened)
    plt.title(plttitle)
    plt.xlabel('Seconds')
    plt.ylabel("Pearson's Correlation Coefficient")
    plt.show()

