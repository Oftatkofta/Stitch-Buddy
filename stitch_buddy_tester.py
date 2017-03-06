from __future__ import print_function
from stitch_buddy import *
from skimage import io
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import numpy as np

indir1 = "Q:\\Jens\\09022017 High vs Low_1\\wells_3_4_5_6\\3\\convert\\"
indir2 = "Q:\\Jens\\161215_HaCaT_LAMP_sorted_8h_post_plate_1\\"
indir3 =r"Q:\Jens\061018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h_1"

outdir=r"Q:\Jens\Stich budy\09022017 High vs Low"


testfile1 = "09022017 High vs Low_1_MMStack_3-Pos_000_000.ome.tif"
testfile2 = "161215_HaCaT_LAMP_sorted_15h_post_plate_1_MMStack_1-Pos_000_000.ome.tif"


loadme_working = os.path.join(indir2, testfile2)

loadme_broken = os.path.join(indir1, testfile1)

print(loadme_working)
print(loadme_broken)


wellNames1 = {3:"3_high-LTG_ColIV_MDMM_10-RDI"}
#, 4:"4_high-LTG_ColIV_MDMM_10-RDI",
 #             5: "5_low-LTG_ColIV_MDMM_10-RDI", 6:"6_low-LTG_ColIV_MDMM_10-RDI"}
#wellNames2 = {4:"161021_HaCaT_YFP-PML_26hStarve_Serum_16MINinterval_t0_2h"}
#wellNames3 = {1:"161018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h"}

wellDict_working = filenamesToDict(indir2)
wellDict_broken = filenamesToDict(indir1, wellNames1)

#wellDict3 = filenamesToDict(indir3, wellNames3)

stitchWellsInRAM(wellDict_broken, indir1, outdir, 0.5)
#stitchWells(wellDict2, indir2, outdir, 0.5)
#stitchWells(wellDict3, indir3, outdir, 0.5)

