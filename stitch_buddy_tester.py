from __future__ import print_function
from Numpy_tiling_from_disk import *





indir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/4wells_2x2-mosaik_2-channel_4-frames_test_1"
indir2 = r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/2Z_2C_2x2gridx2_2T_2"

outdir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/testout"


wellNames = {1:"Ass",
			 2:"hat",
			 3:"100_uM_CBD",
			 4:"150_uM_CBD",
			 }
wellDict1 = filenamesToDict(indir, wellNames)
wellDict2 = filenamesToDict(indir2)

stitchWells(wellDict1, indir, outdir, 0.5)
stitchWells(wellDict2, indir2, outdir, 0.5)

