from __future__ import print_function, division, absolute_import
import javabridge as jb
import bioformats as bf
import os
import re
#import xml.etree.ElementTree
from xml.etree import cElementTree as ET
import skimage.external.tifffile as tiffile
import skimage.io
import matplotlib.pyplot as plt

indir=r"/Volumes/HDD/Huygens_SYNC/Raw OME files for test/2Z_2C_2x2gridx2_2T_2"
outdir=r"/Users/jens_e/Desktop/Python_laboratory/Stitchhack"
#Ingore non-.tif files in indir
filenames = [fname for fname in os.listdir(indir) if ".tif" in fname]
mm_testfile = os.path.join(indir, filenames[0])
ij_testfile = r"/Volumes/HDD/Huygens_SYNC/med_raw.tif"
def split_qn(qn):
    '''Split a qualified tag name or return None if namespace not present'''
    m = re.match('\{(.*)\}(.*)', qn)
    return m.group(1), m.group(2) if m else None

NS_RE = r"http://www.openmicroscopy.org/Schemas/(?P<ns_key>.*)/[0-9/-]"

def get_namespaces(node):
    '''Get top-level XML namespaces from a node.'''
    ns_lib = {'ome': None, 'sa': None, 'spw': None}
    for child in node.iter():
        ns = split_qn(child.tag)[0]
        match = re.match(NS_RE, ns)
        if match:
            ns_key = match.group('ns_key').lower()
            ns_lib[ns_key] = ns
    return ns_lib



jb.start_vm(class_path=bf.JARS, max_heap_size="2G", run_headless=True)
xmlmetadata = bf.get_omexml_metadata(ij_testfile).encode(errors='xmlcharrefreplace')

tree = ET.fromstring(xmlmetadata)

#print(tree.tag)
#print(tree.attrib)

print(tree[0][1].attrib)

for node in tree:
    print(node)

#for child in tree.iter():
    #print(split_qn(child.tag), child.attrib)
    #rank = country.find('rank').text
    #name = country.get('name')
    #print name, rank

jb.kill_vm()

with tiffile.TiffFile(mm_testfile) as tif:
    mm_data = tif.asarray()
    print(tif.micromanager_metadata)
    print(tif.pages)
    print(tif.series)
    print(tif.is_imagej, tif.is_micromanager, tif.is_ome, tif.is_bigtiff)


#tiffile.imshow(mm_data[0, 0, 0, :, :])

with tiffile.TiffFile(ij_testfile) as tif:
    ij_data = tif.asarray()

    try:
        print(tif.micromanager_metadata)
    except:
        print("No MicroManager metadata")
    print(tif.is_imagej, tif.is_micromanager, tif.is_ome, tif.is_bigtiff)

#tiffile.imshow(ij_data[0, :, :])

