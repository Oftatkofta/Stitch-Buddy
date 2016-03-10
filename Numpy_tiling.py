from __future__ import print_function, division, absolute_import
from well_reshape import filenamesToDict
import os

testdir = "/Volumes/HDD/Huygens_SYNC/_SYNC/CollectiveMigrationAnalysis/Examplemovies/160304_well13_128x128"

filenames = os.listdir(testdir)

wellDict = filenamesToDict(filenames)

for k in wellDict.keys():
    for r in range(wellDict[k]['nrows']):
        for c in range(wellDict[k]['ncols']):
            print(wellDict[k]['positions'][(r,c)])