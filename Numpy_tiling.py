from __future__ import print_function, division, absolute_import
from well_reshape import filenamesToDict
import os
import numpy as np
import skimage.io
import time
from matplotlib import pyplot as plt

xpix, ypix, nTimepoints, pixType = 512, 512, 524, "uint16"
#testdir = "/Volumes/HDD/Huygens_SYNC/_SYNC/CollectiveMigrationAnalysis/Examplemovies/160304_well13_128x128"
indir=r"O:\Jens\160304 peprocess"
outdir=r"O:\Jens\160304_processed"
filenames = os.listdir(indir)

wellDict = filenamesToDict(filenames)
wellDict1 = {wellDict.keys()[0]: wellDict[wellDict.keys()[0]]}

def stitchWells(wellDict, inputDir, outputDir):
    for well in wellDict.keys():
        t0=time.time()
        print("Starting on wellID:", well)
        ncols = wellDict[well]['ncols']
        nrows = wellDict[well]['nrows']
        outWidth = xpix*ncols
        outHeight = ypix*nrows
        outArray = np.empty((nTimepoints, outWidth, outHeight), dtype=pixType)

        for c in range(ncols):
            for r in range(nrows):
                startX = (ncols-c-1)*xpix
                startY = r*ypix
                loadme = os.path.join(inputDir, wellDict[well]['positions'][(r,c)])
                print("Working on...", loadme)
                inArray = skimage.io.imread(loadme)
                outArray[:, startY:(startY+ypix), startX:(startX+xpix)] = inArray

        saveme=os.path.join(outputDir, str(well)+"_stitched.tif")
        skimage.io.imsave(saveme, outArray)
        print("Done with wellID: ", well, "in ", round(time.time()-t0,2), " s")



stitchWells(wellDict, indir, outdir)