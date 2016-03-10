from __future__ import print_function, division, absolute_import
from well_reshape import filenamesToDict
import os
import numpy as np
import skimage.io
from matplotlib import pyplot as plt

xpix, ypix, nTimepoints, pixType = 128, 128, 524, "uint16"
testdir = "/Volumes/HDD/Huygens_SYNC/_SYNC/CollectiveMigrationAnalysis/Examplemovies/160304_well13_128x128"
filenames = os.listdir(testdir)

wellDict = filenamesToDict(filenames)

for k in wellDict.keys():
    ncols = wellDict[k]['ncols']
    nrows = wellDict[k]['nrows']
    outWidth = xpix*ncols
    outHeight = ypix*nrows
    arr = np.empty((nTimepoints, outWidth, outHeight), dtype=pixType)

    for c in range(ncols):
        for r in range(nrows):
            startX = (ncols-c-1)*xpix
            startY = r*ypix
            inArray = skimage.io.imread(testdir+"/"+wellDict[k]['positions'][(r,c)])
            arr[:, startY:(startY+ypix), startX:(startX+xpix)] = inArray

skimage.io.imsave(testdir+"/"+str(k)+".tif",arr)