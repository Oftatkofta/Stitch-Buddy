from __future__ import print_function, division, absolute_import
import re
import skimage.io
import os
import cPickle as pickle

#with open("filenames.txt",'r') as f:
#    filenames = f.read()

#filenames = filenames.split("f1")

indir=indir=r"O:\temp"
outdir=r"O:\tempout"

#Ingore non-.tif files in indir
filenames = [fname for fname in os.listdir(indir) if ".tif" in fname]


def filenamesToDict(indir, filenames):
    """
    Transforms a list of filenames into a dictionary.

    Filenames must conform to the following general pattern:
    ..._wellID-xxx_threeDigitRowNumber_threeDigitColNumber...

    :param
    filenames: list of filenames from a multiwell mosaik experiment

    property_dict = 'nrows': (int) No. of rows in well
                    'ncols': (int) No. of columns in well
                    'nChans': (int) No. of channels in image
                    'xpix': (int) No. pixels in X-dimension of images
                    'ypix': (int) No. pixels in Y-dimension of images
                    'dT': (float) No. of time units between frames
                    'timeUnit': (str) time unit
                    'pixres': (float) size of pixels in resolution units
                    'resUnit': (str) spatial resolution unit
                    'pixelType':(str) pixeldepth of image
                    'positions': position_dict gives file at position (row,col)
                    'files':(list) names of files in well
                    'isConcat':(bool) If the sequence is split in to multiple files

    position_dict = {(row, col):filename(s)}



    :return: Dictionary with wellID:property_dict
    """

    wellDict = {}

    isConcat = False

    #Regex used to flexibly identify wells, rows, columns, and split files from filenames

    #Matches any number of digits preceded by "_" and followed by "-"
    well_regex = re.compile("(?<=_)\d+(?=-)")

    #Matches three digits preceded by three digits and "_"
    column_regex = re.compile("(?<=\d{3}_)\d{3}")

    #Matches three digits preceded by three digits and "_"
    row_regex = re.compile("(?<=_)\d{3}(?=_)")

    #Matches a digit preceded by three digits and "_", followed by ".ome"
    concat_regex = re.compile("(?<=\d{3}_\d{3}_)\d(?=\.ome)")

    first_file = skimage.io.imread(os.path.join(indir,filenames[0]))

    print(first_file.shape)

    if first_file.ndim == 4:
        nTimepoints, nChannels, ypix, xpix = first_file.shape

    if first_file.ndim == 5:
        test, nTimepoints, nChannels, ypix, xpix = first_file.shape

    else:
        try:
            nTimepoints, nChannels, ypix, xpix = first_file.shape
        except:
            nTimepoints, ypix, xpix = first_file.shape
            nChannels = 1

    pixelType = first_file.dtype

    for f in filenames:
        #Extract positioning information from filename with regex
        wellID = int(well_regex.search(f).group())
        rowNumber = int(row_regex.search(f).group())
        columnNumber = int(column_regex.search(f).group())
        concat = int(column_regex.search(f).group())
        if concat != None:
            isConcat = True

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
                                'files':[],
                                'isConcat': isConcat
                                }

        #Populate Properties
        wellDict[wellID]['nrows'] = max(rowNumber+1, wellDict[wellID]['nrows'])
        wellDict[wellID]['ncols'] = max(columnNumber+1, wellDict[wellID]['ncols'])

        #List of filenames for the well
        wellDict[wellID]['files'].append(f)

        #Dict with (row, column):(list) filename(s)
        if wellDict[wellID]['positions'].get((rowNumber, columnNumber)) == None:
            wellDict[wellID]['positions'][(rowNumber, columnNumber)]=[f]
        else:
            wellDict[wellID]['positions'][(rowNumber, columnNumber)].append(f)
            wellDict[wellID]['positions'][(rowNumber, columnNumber)].sort()
    return wellDict


wellDict = filenamesToDict(indir, filenames)
for well in wellDict.keys():
    print(well, wellDict[well]['positions'])