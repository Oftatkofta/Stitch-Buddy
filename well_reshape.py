from __future__ import print_function, division, absolute_import
import re
import skimage.io
import os
import cPickle as pickle

#with open("filenames.txt",'r') as f:
#    filenames = f.read()

#filenames = filenames.split("f1")

indir=r"O:\SYNC with anlysis\160321 scratch_factors"
outdir=r"O:\SYNC with anlysis\160321 scratch_factors_out"
filenames = os.listdir(indir)


def filenamesToDict(indir, filenames, useBioformats=False):
    """
    Transforms a list of filenames into a dictionary.

    Filenames must conform to the following general pattern:
    xxxxx_wellID-xxx_threeDigitRowNumber_threeDigitColNumber.xxx

    :param
    filenames: list of filenames from a multiwell mosaik experiment
    useBioformats: (bool) Should the files be read using bioformats?

    property_dict = 'nrows': (int) No. of rows in well
                    'ncols': (int) No. of columns in well
                    'nChans': (int) No. of channels in image, defaults to 1
                    'xpix': (int) No. pixels in X-dimension of images
                    'ypix': (int) No. pixels in Y-dimension of images
                    'dT': (float) No. of time units between frames
                    'timeUnit': (str) time unit
                    'pixres': (float) size of pixels in resolution units
                    'resUnit': (str) spatial resolution unit
                    'pixelType':(str) pixeldepth of image
                    'positions': position_dict gives file at position (row,col)
                    'files':(list) names of files in well

    position_dict = (row, col):filename



    :return: Dictionary with wellID:property_dict
    """
    if useBioformats:
        import javabridge as jb
        import bioformats as bf
        jb.start_vm(class_path=bf.JARS, run_headless=True)

    wellDict = {}

    #Regex used to flexibly identify wells, rows, and columns from filenames

    #Matches any number of digits preceded by "_" and followed by "-"
    well_regex = re.compile("(?<=_)\d+(?=-)")

    #Matches three digits preceded by "_" and followed by "."
    column_regex = re.compile("(?<=_)\d{3}(?=\.)")

    #Matches three digits preceded by "_" and followed by "_"
    row_regex = re.compile("(?<=_)\d{3}(?=_)")
    first_file = skimage.io.imread(os.path.join(indir,filenames[0]))
    if first_file.ndim >3:
        nTimepoints, nChannels, ypix, xpix = first_file.shape
    else:
        nTimepoints, ypix, xpix = first_file.shape
        nChannels = 1
    pixelType = first_file.dtype
    for f in filenames:
        #Extract positioning information from filename with regex
        wellID = int(well_regex.search(f).group())
        rowNumber = int(row_regex.search(f).group())
        columnNumber = int(column_regex.search(f).group())

        #If there is no key for wellID in wellDict -> create a dict of properties
        if wellDict.get(wellID) == None:
            wellDict[wellID] = {'nrows':1,
                                'ncols':1,
                                'nChannels':int(nChannels),
                                'xpix':int(xpix),
                                'ypix':int(ypix),
                                'timepoints':int(nTimepoints),
                                'dT':0,
                                'timeunit':None,
                                'pixres':0,
                                'pixelType':str(pixelType),
                                'positions':{},
                                'files':[]
                                }

        #Populate Properties
        wellDict[wellID]['nrows'] = max(rowNumber+1, wellDict[wellID]['nrows'])
        wellDict[wellID]['ncols'] = max(columnNumber+1, wellDict[wellID]['ncols'])

        #List of filenames for the well
        wellDict[wellID]['files'].append(f)

        #Dict with (row, column):filename
        wellDict[wellID]['positions'][(rowNumber, columnNumber)]=f

    return wellDict


wellDict = filenamesToDict(indir, filenames)
for well in wellDict.keys():
    print(well, wellDict[well])