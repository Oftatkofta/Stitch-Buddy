from __future__ import print_function
from stitch_buddy import *





indir1 = r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/4wells_2x2-mosaik_2-channel_4-frames_test_1/"
#indir2 = r"Q:\Jens\161021_HaCaT_YFP-PML_26hStarve_Serum_16MINinterval_t0_2h_1"
#indir3 =r"Q:\Jens\061018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h_1"

outdir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/testout/"


wellNames1 = {1:"160808_asshat",
              2: "bajszz",
              3: "-interval_t0_110min",
              4: "sccoopp"}
#wellNames2 = {4:"161021_HaCaT_YFP-PML_26hStarve_Serum_16MINinterval_t0_2h"}
#wellNames3 = {1:"161018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h"}

wellDict1 = filenamesToDict(indir1, wellNames1)
#wellDict2 = filenamesToDict(indir2, wellNames2)
#wellDict3 = filenamesToDict(indir3, wellNames3)
print(wellDict1["sccoopp"]["OME-XML"])
#stitchWellsInRAM(wellDict1, indir1, outdir, 0.5)

