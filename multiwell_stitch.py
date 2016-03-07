from __future__ import print_function, division
import os.path
import javabridge as jb
import bioformats as bf
from matplotlib import pyplot as plt
import tifffile as tiff
import re
import numpy as np
import time
from skimage import io

#Python-bioformats includes a buggy bioformats .jar, so it is replaced
bioformats_package_path = ["C:\\bioformats\\bioformats_package.jar"]
jb.start_vm(class_path=bf.JARS[:-1]+bioformats_package_path, max_heap_size="2G", run_headless=True)


"""
List files in Dir
160304_HaCaTwt_12well_T0_1h_8min_1_MMStack_24-Pos_000_000.ome
regex = r"\d+\D+\d{3}[_]\d{3}"  #findall -> 24-Pos_000_000
for file in files
"""
directory = "O:\Jens\\160304_HaCaTwt_12well_T0_1h_8min_1"

filenames = os.listdir(directory)

f=open("filenames.txt",'w')
f.writelines(filenames)
f.close()

def get_imageDimensions(filename):
    """
    Returns a tuple of the basic image dimensions from an ome-file.

    Args:
        filename: (str)

    Returns: (nTimepoints, nChannels, nXpixels, nYpixels, pixelType)

    """
    metadata = bf.OMEXML(bf.get_omexml_metadata(filename))
    nTimepoints = metadata.image().Pixels.get_SizeT()
    nChannels = metadata.image().Pixels.get_SizeC()
    nXpixels = metadata.image().Pixels.get_SizeX()
    nYpixels = metadata.image().Pixels.get_SizeY()
    pixelType = metadata.image().Pixels.get_PixelType()

    return nTimepoints, nChannels, nXpixels, nYpixels, pixelType

def OMEtoArray(filename):
    """
    Args:
        file: (str) Path to OME-TIFF file

    Returns: numpy array of shape (time, X, Y)

    """
    #Get basic dimensions of file
    nTimepoints, nChannels, nXpixels, nYpixels, pixelType = get_imageDimensions(filename)

    #Create empty array to fill with data from file
    outStack = np.empty((nTimepoints, nXpixels, nYpixels), dtype=pixelType)

    #Read frames in to outStack
    with bf.ImageReader(filename) as reader:
        for t in range(nTimepoints):
            outStack[t] = reader.read(t=t, rescale=False)
        reader.close()
    return outStack


loadme1 = os.path.join(directory, filenames[0])
loadme2 = os.path.join(directory, filenames[1])

stack1 = OMEtoArray(loadme1)
stack2 = OMEtoArray(loadme2)

stack = np.concatenate((stack2,stack1),axis=2)
plt.imshow(stack[100], cmap='viridis')
plt.show()

jb.kill_vm()

