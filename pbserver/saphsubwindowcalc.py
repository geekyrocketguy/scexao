import numpy as np

def calc(dim):
    for i in np.arange(320)+1.:
        if (i*dim*10/64./64.) % 1. == 0.:
            print i
