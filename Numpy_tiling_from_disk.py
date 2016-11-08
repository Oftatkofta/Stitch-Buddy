from __future__ import print_function, division, absolute_import
import re
import os
import numpy as np
import tiffile as tiffile
import time
from fractions import Fraction

def filenamesToDict(indir, wellNameDict=None):

    """
    Transforms the .tif-files in a directory into a dictionary.

    Filenames must conform to the following general pattern:
    ..._wellID-xxx_threeDigitRowNumber_threeDigitColNumber...

    :param
    indir: path to directory with files from a multiwell mosaic experiment
    wellNameDict: Dictionary of names to give the wells from the wellID-number

    property_dict = 'nrows': (int) No. of rows in well
                    'ncols': (int) No. of columns in well
                    'nChans': (int) No. of channels in image
                    'xpix': (int) No. pixels in X-dimension of images
                    'ypix': (int) No. pixels in Y-dimension of images
                    'frame_interval': (float) No. of time units between frames
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
    # Ingore non-.tif files in indir
    filenames = [fname for fname in os.listdir(indir) if ".tif" in fname]

    wellDict = {}

    isConcat = False

    #Regex used to flexibly identify wells, rows, columns, and split files from filenames

    #Matches any number of digits preceded by "MMStack_" and followed by "-"
    well_regex = re.compile("(?<=MMStack_)\d+(?=-)")

    #Matches three digits preceded by three digits and "_"
    column_regex = re.compile("(?<=\d{3}_)\d{3}")

    #Matches three digits preceded by three digits and "_"
    row_regex = re.compile("(?<=_)\d{3}(?=_)")

    #Matches a digit preceded by three digits and "_", followed by ".ome"
    concat_regex = re.compile("(?<=\d{3}_\d{3}_)\d(?=\.ome)")

    first_file = os.path.join(indir, filenames[0])

    with tiffile.TiffFile(first_file) as tif:
        meta_data = tif.micromanager_metadata
        frames = len(tif.pages)
        page=tif[0]
        pixres=page.tags['x_resolution'].value #Assumes equal x/y resolution
        resoloution_unit = page.tags['resolution_unit'].value
        pixelDepth = page.tags['bits_per_sample'].value



    frame_interval = meta_data['summary']['WaitInterval']
    ypix, xpix = meta_data['summary']['Height'], meta_data['summary']['Width']
    nChannels =  meta_data['summary']['Channels']
    nSlices =  meta_data['summary']['Slices']
    nTimepoints = frames/(nChannels*nSlices)
    resoloution_unit =  {1: 'none', 2: 'inch', 3: 'centimeter'}[resoloution_unit]

    for f in filenames:
        #Extract positioning information from filename with regex
        wellID = int(well_regex.search(f).group())


        if wellNameDict != None:
            wellID = wellNameDict[wellID]


        rowNumber = int(row_regex.search(f).group())
        columnNumber = int(column_regex.search(f).group())
        concat = concat_regex.search(f)
        if concat != None:
            isConcat = True

        #If there is no key for wellID in wellDict -> create a dict of properties
        if wellDict.get(wellID) == None:
            wellDict[wellID] = {'nrows':1,
                                'ncols':1,
                                'nChannels':int(nChannels),
                                'nSlices': int(nSlices),
                                'xpix':int(xpix),
                                'ypix':int(ypix),
                                'nTimepoints':int(nTimepoints),
                                'frame_interval':frame_interval,
                                'timeunit':'ms',
                                'pixel_resolution':pixres, #resolution stored as rational in tif tag
                                'resoloution_unit': resoloution_unit,
                                'pixelDepth':pixelDepth,
                                'positions':{},
                                'files':[],
                                'isConcat':isConcat,
                                'OME-XML':None
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
            wellDict[wellID]['isConcat'] = isConcat

    return wellDict

def bin_ndarray(ndarray, new_shape, operation='mean'):
	"""
	Bins an ndarray in all axes based on the target shape, by summing,
		averaging, or returning min/max pixel value

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
		nChans = wellDict[well]['nChannels']
		nTimepoints = wellDict[well]['nTimepoints']
		nSlices = wellDict[well]['nSlices']
		frame_interval = wellDict[well]['frame_interval']
		time_unit  = wellDict[well]['timeunit']
		pixelDepthDict = {8: "uint8", 16:"uint16", 32:"float32"}
		pixType = pixelDepthDict[wellDict[well]['pixelDepth']]
		xpix = wellDict[well]['xpix']
		ypix = wellDict[well]['ypix']
		pixel_resolution = wellDict[well]['pixel_resolution']

		outWidth = xpix*ncols
		outHeight = ypix*nrows

		if resizeTo != None:
			old_resolution = pixel_resolution[1]/float(pixel_resolution[0])
			print("old resolution; {}, rational: {}".format(old_resolution,
			                                                pixel_resolution))
			new_resolution = old_resolution*resizeTo
			rational_new_resolution = Fraction(new_resolution).limit_denominator()

			pixel_resolution = (rational_new_resolution.denominator,
								rational_new_resolution.numerator)
			print("new resolution; {}, rational: {}".format(new_resolution,
			                                                pixel_resolution))
			outArray = np.empty(
				(nTimepoints, nSlices, nChans, outHeight*resizeTo, outWidth*resizeTo),
				dtype=pixType)

		else:
			outArray = np.empty(
				(nTimepoints, nSlices, nChans, outHeight, outWidth),
				dtype=pixType)

		full_frame_buffer = np.empty((1, nSlices, nChans, outHeight, outWidth),
				dtype=pixType)

		print("output well array shape: {}, full_frame_array shape; {}".format(outArray.shape, full_frame_buffer.shape) )


		for frame in range(nTimepoints):
			t1 = time.time()
			print("Working on frame: {}".format(frame))
			for row in range(nrows):
				for col in range(ncols):
					#get name of file at this row, col position
					loadme = os.path.join(inputDir,
										  wellDict[well]['positions'][(row, col)][0]
										  )
					#X/Y coordinates for insertion of subframe
					startX = (ncols-col-1)*xpix
					startY = row*ypix

					#is_ome is set to False so that all .asarray calls will return the same shape
					with tiffile.TiffFile(loadme, is_ome = False) as tif:
						#only get current frame with channels and slices as an array
						if nSlices == 1:
							slice_to_load = slice(frame * nChans,
							                      (frame + 1) * nChans)
						else:
							slice_to_load = slice(frame*(nChans+nSlices),
							                      (frame+1)*(nChans+nSlices))
						print(slice_to_load)
						inArray = tif.asarray(key=slice_to_load)
						print("Frame: {}, Row: {}, Col: {} loaded from file: {} array shape: {}, loadtime: {} s".format(frame, row, col, tif.filename, inArray.shape, round(time.time() - t1)) )

						try:
							full_frame_buffer[:,:,:, startY:(startY + ypix), startX:(startX + xpix)] = inArray
						except:
							inArray = np.reshape(inArray, (1, nSlices, nChans, xpix, ypix))
							print("Input array reshaped to: {}".format(inArray.shape))
							full_frame_buffer[:,:,:, startY:(startY + ypix), startX:(startX + xpix)] = inArray

						print("Input array appended to fullframebuffer. Elapsed time for well: {}".format(round(time.time() - t0)) )

			if resizeTo != None:

				resized_frame_buffer = bin_ndarray(full_frame_buffer, (1, nSlices, nChans, outHeight*resizeTo, outWidth*resizeTo))
				print("full_fame_buffer resized to {}".format(resized_frame_buffer.shape))

				outArray[frame, :, :, :, :] = resized_frame_buffer
			else:
				outArray[frame,:,:,:,:] = full_frame_buffer

			print("OutArray appended with frame: {}, frame took {} s to stitch.".format(frame, round(time.time() - t1 )))

		saveme=os.path.join(outputDir, str(well)+"_stitched.tif")
		bigTiffFlag = outArray.size * outArray.dtype.itemsize > 2000 * 2 ** 20
		print("bigTillFlag set to:", bigTiffFlag, "Saving output...(may take a while)")

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

		print("Done with wellID: ", well, "in ", round(time.time() - t0, 2),
			  " s")
	print("All done, it took ", round(time.time() - tStart, 2), " s in total!")


