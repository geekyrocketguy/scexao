#Displays/saves images showing why RRR has a post-reset settling issue

import matplotlib.pyplot as plt
import numpy as np
import pyfits
import pdb

def main():
    file = '/home/scexao/pizzabox/SCExAO/pbserver.main.64.data/seansdatafiles/noise_rrr-1-0.fits'

    img = pyfits.getdata(file)

    imafter = img[400, 256:, :] 
    imbefore = img[399, :256, :]
    diff = imafter - imbefore

    plt.imshow(imbefore, vmin=18e3, vmax=25e3, interpolation='none')
    plt.colorbar(shrink=0.9)
    plt.title('Raw Image Following Reset', fontsize=20)
    plt.tight_layout()
    #plt.show()
    plt.savefig('../figures/rrrfail_after.png', bbox_inches='tight')
    plt.close()

    plt.imshow(imafter, vmin=18e3, vmax=25e3, interpolation='none')
    plt.colorbar(shrink=0.9)
    plt.title('Raw Image Preceding Reset', fontsize=20)
    plt.tight_layout()
    #plt.show()
    plt.savefig('../figures/rrrfail_before.png', bbox_inches='tight')
    plt.close()

    plt.imshow(diff, vmin=0, vmax=400, interpolation='none')
    plt.colorbar(shrink=0.9)
    plt.title('Difference Image', fontsize=20)
    plt.tight_layout()
    #plt.show()
    plt.savefig('../figures/rrrfail_diff.png', bbox_inches='tight')
    plt.close()

