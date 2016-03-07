__author__ = 'jens_e'
import math
#fin = open("TileConfiguration_gen.txt", 'r').readlines()

fout=open('TileConfiguration.txt', 'w')
basename = "160213_H2B_starve_serumFree_T0_10min_1_MMStack_1-Pos_"

nrows = 10
ncols = 10
xpix=256
ypix=256
alpha=0
pixel_overlap=0

head="# Define the number of dimensions we are working on\ndim = 2\n\n# Define the image coordinates\n"
suffix = '.ome.tif; ; ('
#dy=0
#dx=0

def straightStitch(head, basename, overlap):
    out=head
    for r in range(nrows):
        y=(nrows-1-r)*(ypix-overlap)
        for c in range((ncols-1),-1,-1):
            x=(((ncols-1)-c)*(xpix-overlap))
            out += basename + str(r).zfill(3)+'_'+str(c).zfill(3) + suffix + str(float(x)) + ', ' + str(float(y))+')\n'

        #dy=2

    return out

def shadeStitch(head, overlap):
    out=head
    for r in range(10):
        y=(r*(ypix-overlap))
        for c in range((10-1),-1,-1):
            x=(((10-1)-c)*(xpix-overlap))
            out += "NORM_MED_10x0.3NA_660nm_filterset50-A647_Andor_MaTek.tif; ; (" + str(float(x)) + ', ' + str(float(y))+')\n'

        #dy=2

    return out



def Stitch(head, basename, overlap):
    out=head
    ydiff = (ncols-1)*(xpix-overlap)
    dy = ydiff / float(ncols-1)
    xdiff = (nrows-1)*(ypix-overlap)
    dx = xdiff / float(nrows-1)

    for r in range(0,nrows):
        yoffset=0
        xoffset=r*dx
        for c in range((ncols-1),-1,-1):
            x=int((((ncols-1)-c)*(xpix-overlap)+xoffset))
            y=int((r*(ypix-overlap))+yoffset)
            out += basename + str(c).zfill(3) + '_' + str(r).zfill(3) + '.ome.tif; ; (' + str(float(x)) + ', ' + str(float(y))+')\n'
            yoffset += dy

    return out


def scewedStitch(head, basename, overlap):
    out=head
    if alpha != 0:
        a=math.tan(math.radians(alpha))
    else:
        a=1
    ydiff=(a*(ncols-1)*(xpix-overlap))
    dy = ydiff/float(ncols-1)
    xdiff=a*(nrows-1)*(ypix-overlap)
    dx=xdiff/float(nrows-1)

    for r in range(0,nrows):
        yoffset=0
        xoffset=r*dx
        for c in range((ncols-1),-1,-1):
            x=int((((ncols-1)-c)*xpix)+xoffset)
            y=int((r*ypix)+yoffset)

            out += basename + str(c).zfill(3) + '_' + str(r).zfill(3) + '.ome.tif; ; (' + str(float(x)) + ', ' + str(float(y))+')\n'
            yoffset += dy

    return out

def orcaStraightStitch(head, basename, overlap, x_correction=0, y_correction=0):
    out=head

    for r in range(nrows):
        y=(r*(ypix-overlap))
        for c in range(ncols):

            x=(ncols-1-c)*(xpix-overlap)+r*x_correction
            out += basename + str(r).zfill(3)+'_'+str(c).zfill(3) + suffix + str(float(x)) + ', ' + str(float(y))+')\n'
            y+=y_correction
    return out


test=orcaStraightStitch(head, basename, 0)
fout.write(test)
fout.close()

