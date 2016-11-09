from __future__ import print_function
from stitch_buddy import *





indir1 = r"Q:\Jens\061017_HacaT_mCherryH2B_YFPPML_24hStarve_Serum_16min_t0_110min_1"
indir2 = r"Q:\Jens\161021_HaCaT_YFP-PML_26hStarve_Serum_16MINinterval_t0_2h_1"
indir3 =r"Q:\Jens\061018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h_1"

outdir=r"O:\tempout"


wellNames1 = {3:"161017_HacaT_mCherry-H2B_YFP-PML_24hStarve_Serum_16min-interval_t0_110min"}
wellNames2 = {4:"161021_HaCaT_YFP-PML_26hStarve_Serum_16MINinterval_t0_2h"}
wellNames3 = {1:"161018_HacaT_YFPPML_72hStarve_Serum_16min_t0_4h"}

wellDict1 = filenamesToDict(indir1, wellNames1)
wellDict2 = filenamesToDict(indir2, wellNames2)
wellDict3 = filenamesToDict(indir3, wellNames3)

stitchWells(wellDict1, indir1, outdir, 0.5)
stitchWells(wellDict2, indir2, outdir, 0.5)
stitchWells(wellDict3, indir3, outdir, 0.5)

