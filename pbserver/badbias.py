import numpy as np
from astropy.io import fits as pf
from PIL import Image
import time
import os
import sys
home = os.getenv('HOME')
sys.path.append(home+'/src/lib/python/')
from scexao_shm import shm

cam = shm("/tmp/pbimage.im.shm")
ndark=1000

while(True):
    temp = np.squeeze(cam.get_data(True, True, timeout = 1.).astype(float)) #dimensions 1,y,x
    darkcube = np.zeros(np.append(ndark, np.shape(temp)))

    for idark in range(ndark):
        darkcube[idark, :, :] = np.squeeze(cam.get_data(True, True, timeout = 1.).astype(float))
        
    dark = np.median(darkcube, axis=0)
    newXImSize = np.shape(dark)[1]
    newYImSize = np.shape(dark)[0]
    fname = home+"/pizzabox/SCExAO/pbserver.main/saphcamdarks/"+\
            time.strftime("%Y%m%d%H%M%S")+'_'+str(newXImSize)+'_'+str(newYImSize)+'.fits'
    pf.writeto(fname, dark, clobber=True)

    print "Wrote", fname
    time.sleep(600)
