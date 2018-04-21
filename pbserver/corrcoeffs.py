import matplotlib.pyplot as plt
import numpy as np
import pyfits
import pdb
import glob
from matplotlib.colors import LogNorm
#import pickle
import random
#OLD, NO LONGER SUPPORTED
#calculates Pearson's correlation coefficients for on-sky SAPHIRA data.
# corr_coeffs_2.py uses less memory and is currently written for internal
# data (but this could easily be changed).

def main(restore=False):
    if restore != True:
        print "Calculating data."

        filelist = sorted(glob.glob('/media/data/20170911/saphira/processed/pbimage*fits'))
        
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

        loc = np.where(annulus == 1)

        #n = 0
        img = np.zeros((10e3*len(filelist), 128, 128))
        #img = np.zeros((10e3, 128, 128))

        for i in range(int(len(img)/1e4)):
            print "Populating image cube, " + \
                str(int(round(float(i)/(len(img)/1e4)*100))) + "% done."
            img[1e4*i : 1e4*(i+1)] = pyfits.getdata(filelist[i])

        framenos = range(0, len(img)-100, 1)
        random.shuffle(framenos) #so it samples randomly in time
    
        separations = np.append(range(1, 168), range(200, 1680*5, 8))
        separations = np.append(separations, range(1680*6, 604800, 1680))

        #save separations
        np.save('corr_coeff_separations_onsky.npy', separations)
        #file1 = open('corr_coeff_separations.txt','w') 
        #pickle.dump(separations, file1)
        #file1.close()

        #move the selected parts of img into im2 for quicker processing
        size_selection = np.shape(loc)[1] #number of pixels in annulus
        img2 = np.zeros((len(img) , size_selection))
        stds = np.zeros(len(img))
        for i in range(len(img2)):
            if i%1000==0:
                print "Reducing amount of data, " + \
                    str(int(round(float(i)/len(img2)*100.))) + '% done.'

            img2[i] = img[i][loc]
            img2[i] -= np.mean(img2[i])
            stds[i] = np.std(img2[i])
        del img #save ram
            
        matrix = np.zeros((len(framenos), len(separations))) #i, sep

        for no in range(len(framenos)):
            i = framenos[no]
            #sel1 = img[i][loc] #faster to put it up here
    
            if no%10 == 0:
                #print str(int(round(float(i)/len(img)*100.)))+ "% done."
                print str(round(float(no)/len(framenos)*100.*1000.)/1000.)+ "% done."

            for s in range(len(separations)):
                #if s%100==0: print float(s)/len(separations)*100.
    
                sep = separations[s]

                if i+sep >= len(img2): continue #skip to next iteration

                #exclude period in which saphira was blanked off
                if (i > (1e4*8+6000)) & (i < (1e4*10+900)):
                    continue
                if (i+sep > (1e4*8+6000)) & (i+sep < (1e4*10+900)):
                    continue

                #sel2 = img[i+sep][loc]

                #clear definition
                #rho = np.sum((sel1 - np.mean(sel1)) * (sel2 - np.mean(sel2)) ) / \
                #             (np.std(sel1) * np.std(sel2) * np.size(sel1))

                #efficient definition
                rho = np.sum(img2[i]  * img2[i+sep]) / \
                             (stds[i] * stds[i+sep] * size_selection)

                matrix[no, s] = rho
    
            #save matrix
            if no % 10000 == 0:
                print "saving"
                np.save('corr_coeff_matrix_onsky.npy', matrix)
                print "saved"
                print
            #print "saving matrix"
            #file = open('corr_coeff_matrix.txt','w') 
            #pickle.dump(matrix, file)
            #file.close()
            #print done
            #print

    else: #restore data from text file
        print "Restoring data from file."

        matrix = np.load('corr_coeff_matrix_onsky.npy')
        separations = np.load('corr_coeff_separations_onsky.npy')
        print "done."
        #file = open('corr_coeff_matrix.txt','r') 
        #matrix = pickle.load(file)
        #file.close()

        #file1 = open('corr_coeff_separations.txt','r') 
        #separations = pickle.load(file1)
        #file1.close()


    #matrix[10,0] = 2
    if 0:
        plt.imshow(matrix, interpolation='none')
        plt.colorbar()
        plt.xlabel('separation')
        plt.ylabel('i')
        plt.show()

    #flattened = np.mean(matrix, 0)

    flattened = np.zeros(len(separations))
    for i in range(len(flattened)):
        flattened[i] = np.mean((matrix[:,i])[np.where(matrix[:,i] != 0)])

    plt.plot(separations/1680., flattened)
    plt.xlabel('seconds')
    plt.ylabel("Pearson's Correlation Coefficient")
    plt.show()

