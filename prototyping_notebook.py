
# coding: utf-8

# In[23]:

from __future__ import print_function, absolute_import, division
from PIL import Image, ImageSequence, TiffTags
import tiffile.py
import os


# In[6]:

indir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/2Z_2C_2x2gridx2_2T_2" 
outdir=r"/Users/jens_e/Desktop/Python_laboratory/Stitchhack"                   
#Ingore non-.tif files in indir                                                
filenames = [fname for fname in os.listdir(indir) if ".tif" in fname]          
mm_testfile = os.path.join(indir, filenames[0])                                
ij_testfile = r"/Volumes/HDD/Huygens_SYNC/med_raw.tif"                         


# In[27]:

mm_test = Image.open(mm_testfile)
ij_test = Image.open(ij_testfile)
print(mm_testfile, mm_test.format, mm_test.size, mm_test.mode, mm_test.info)
print(ij_testfile, ij_test.format, ij_test.size, ij_test.mode, ij_test.info)

mm_iterator = ImageSequence.Iterator(mm_test)
for frame in mm_iterator:
    print(frame.info, frame.size)
print(mm_iterator[2])


# In[ ]:




# In[ ]:



