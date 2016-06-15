from __future__ import print_function, division, absolute_import
from frankenScope_well_reshape import filenamesToDict
import os
import numpy as np
from skimage.transform import rescale, resize
try :
    import Stitchit.Stitchit.tiffile_mod as tiffile

except ImportError:
    import Stitchit.tiffile_mod as tiffile

import time
from matplotlib import pyplot as plt


#indir=r"O:\Jens\160408_HaCaT-H2B_12well_Calcium_T0_1h_2"
#indir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/2C_2x2gridx2_2T_1"
#outdir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/"
#outdir=r"O:\temp"

#indir=r"O:\Jens\160610_HaCaT-H2B_Insert_serum_EGF_4"
#outdir=r"O:\tempout"



wellNames = {2:"20_uM_CBD",
             4:"50_uM_CBD",
             5:"100_uM_CBD",
             6:"150_uM_CBD",
             7:"4_mM_EGTA",
             8:"2_mM_EGTA",
             10:"1_mM_EGTA",
             11:"0.5_mM_EGTA",
             12:"serum",
             14:"no_serum",
             15:"20ug_EGF_EpiLife",
             }



#wellDict = filenamesToDict(indir, wellNames)
#wellDict1 = {wellDict.keys()[0]: wellDict[wellDict.keys()[0]]}

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

def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean', 'max', 'min']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [int(l) for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray

def stitchWells(wellDict, inputDir, outputDir, resizeTo=None):

    tStart=time.time()

    for well in wellDict.keys():
        t0=time.time()
        print("Starting on wellID:", well)
        ncols = wellDict[well]['ncols']
        nrows = wellDict[well]['nrows']
        nChans, nTimepoints = wellDict[well]['nChannels'], wellDict[well]['nTimepoints']
        nSlices = wellDict[well]['nSlices']
        frame_interval, time_unit  = wellDict[well]['frame_interval'], wellDict[well]['timeunit']
        pixType = wellDict[well]['pixelType']
        xpix, ypix, pixel_resolution = wellDict[well]['xpix'], wellDict[well]['ypix'], wellDict[well]['pixel_resolution']
        outWidth = xpix*ncols
        outHeight = ypix*nrows

        outArray = np.empty((nTimepoints, nSlices, nChans, outHeight, outWidth), dtype="uint16")

        print("well array shape:", outArray.shape)

        for r in range(nrows):
            for c in range(ncols):
                startX = (ncols-c-1)*xpix
                startY = r*ypix
                loadme = os.path.join(inputDir, wellDict[well]['positions'][(r,c)][0])
                print("Working on: ", str(loadme))

                with tiffile.TiffFile(loadme) as tif:
                    inArray = tif.asarray()
                    print("Array loaded, shape: ", inArray.shape, "Elapsed time: ", time.time() - t0)

                try:
                    outArray[:,:,:, startY:(startY+ypix), startX:(startX+xpix)] = inArray
                except:
                    inArray = np.reshape(inArray, (nTimepoints, nSlices, nChans, xpix, ypix))
                    print("Array reshaped to: ", inArray.shape, "Elapsed time: ", time.time() - t0)
                    outArray[:,:,:, startY:(startY + ypix), startX:(startX + xpix)] = inArray

                print("Array appended to OutArray. Elapsed time: ", time.time() - t0)

        saveme=os.path.join(outputDir, str(well)+"_stitched.tif")

        if resizeTo != None:
            print("Resizing outArray...")
            t=time.time()
            outArray = bin_ndarray(outArray, ((nTimepoints, nSlices, nChans,
                                               outHeight*resizeTo,
                                               outWidth*resizeTo))).astype("uint16")
            print("Done!", time.time()-t, "s")
        bigTiffFlag = outArray.size * outArray.dtype.itemsize > 2000 * 2 ** 20
        print("bigTillFlag set to:", bigTiffFlag)

        metadata = {"zStack" : bool(1-nSlices),
                    "unit":"um",
                    "tunit": time_unit,
                    "finterval" : frame_interval
                    }

        save_data = {"bigtiff" : bigTiffFlag,
                    "imagej" : True,
                    "resolution" : (pixel_resolution, pixel_resolution),
                    "metadata" : metadata
                    }
        tiffile.imsave(saveme, outArray, **save_data)
        print("Done with wellID: ", well, "in ", round(time.time()-t0,2), " s")
    print("All done, it took: ", round(time.time() - tStart, 2), " s")

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

#stitchWells(wellDict, indir, outdir, 0.25)