import pyfits
import matplotlib.pyplot as plt
import numpy as np

def check(filename):

    stddevs = np.array([])
    for i in range(10, 30):

        filenum = str(i)
        while len(filenum) < 4:
            filenum = '0'+filenum
        filenum2 = str(i+1)
        while len(filenum2) < 4:
            filenum2 = '0'+filenum2
        

        cds = pyfits.getdata(filename + filenum + '.fits') - \
              pyfits.getdata(filename + filenum2 + '.fits')
        stddevs = np.append(stddevs, np.std(cds, ddof=1))

        if i%5==0:
            plt.imshow(cds, interpolation='none', vmin=0, vmax=300)
            plt.colorbar()
            plt.show()

    print "Average standard deviation of CDS is", np.median(stddevs)
