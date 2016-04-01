from __future__ import print_function, division, absolute_import
from well_reshape import filenamesToDict
import os
import numpy as np
import skimage.io
import time
from matplotlib import pyplot as plt

#xpix, ypix, nTimepoints, pixType = 256, 256, 406, "uint8"
#testdir = "/Volumes/HDD/Huygens_SYNC/_SYNC/CollectiveMigrationAnalysis/Examplemovies/160304_well13_128x128"
indir=r"O:\SYNC with anlysis\160321 scratch_factors"
outdir=r"O:\SYNC with anlysis\160321 scratch_factors_out"
filenames = os.listdir(indir)

wellDict = filenamesToDict(indir, filenames)
wellDict1 = {wellDict.keys()[0]: wellDict[wellDict.keys()[0]]}

#flatfield = skimage.io.imread(r"C:\Users\Franken_Scope\Desktop\160126 Flatfield correction\20X 0.5NA\MED_20x 0.5NA_DICII_4x4 bin.tif")
#background = skimage.io.imread(r"C:\Users\Franken_Scope\Desktop\160126 Flatfield correction\Dark references\160205_dark_frame_4x4_bin.tif")
#flatfield = skimage.io.imread(r"O:\Jens\160304_processed_ffcorr\MED_20x 0.5NA_DICII_4x4 bin.tif")
#flatfield_image = np.divide(np.asarray(flatfield, dtype='float'), np.amax(flatfield))

def flatfileld_correction(frame, flatfield_image):
    """

    Args:
        frame:
        backgroun:
        flatfield_image: background subtracted already

    Returns:

    """

    return np.asarray(np.multiply(frame, flatfield_image), dtype='uint16')

def normalize_frame(frame):
    return np.divide(np.asarray(frame, dtype='float'), np.amax(frame,axis=0))

def stitchWells(wellDict, inputDir, outputDir):
    for well in wellDict.keys():
        t0=time.time()
        print("Starting on wellID:", well)
        ncols = wellDict[well]['ncols']
        nrows = wellDict[well]['nrows']
        nChans, nTimepoints = wellDict[well]['nChannels'], wellDict[well]['timepoints']
        pixType = wellDict[well]['pixelType']
        xpix, ypix = wellDict[well]['xpix'], wellDict[well]['ypix']
        outWidth = xpix*ncols
        outHeight = ypix*nrows
        outArray = np.empty((nTimepoints, nChans, outHeight, outWidth), dtype=pixType)
        print("well array shape:", outArray.shape)
        for r in range(nrows):
            for c in range(ncols):
                startX = (ncols-c-1)*xpix
                startY = r*ypix
                loadme = os.path.join(inputDir, wellDict[well]['positions'][(r,c)])
                print("Working on...", loadme)
                inArray = skimage.io.imread(loadme)
                #inArray = flatfileld_correction(inArray, flatfield_image)
                outArray[:,:, startY:(startY+ypix), startX:(startX+xpix)] = inArray

        saveme=os.path.join(outputDir, str(well)+"_stitched.tif")
        skimage.io.imsave(saveme, outArray, plugin='tifffile', metadata={'axes': 'TCYX'})
        print("Done with wellID: ", well, "in ", round(time.time()-t0,2), " s")

def getFlatfieldWells(wellDict, inputDir, outputDir):
    for well in wellDict.keys():
        t0=time.time()
        print("Starting on wellID:", well)
        ncols = wellDict[well]['ncols']
        nrows = wellDict[well]['nrows']
        #outArray = np.ones((nTimepoints, ypix, xpix), dtype='float')
        #print("well array shape:", outArray.shape)
        for r in range(nrows):
            for c in range(ncols):
                loadme = os.path.join(inputDir, wellDict[well]['positions'][(r,c)])
                print("Working on...", loadme)
                inArray = skimage.io.imread(loadme)
                outArray = normalize_frame(inArray)


        saveme=os.path.join(outputDir, str(well)+"_nomalized.tif")
        skimage.io.imsave(saveme, outArray, )
        print("Done with wellID: ", well, "in ", round(time.time()-t0,2), " s")


stitchWells(wellDict, indir, outdir)