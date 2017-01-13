from __future__ import print_function
#import openpiv.tools
import openpiv.process
import openpiv.filters
import openpiv.scaling
import openpiv.validation
import pylab as plt
import numpy as np
import math
from PIL import Image
from matplotlib.colors import LogNorm

try:
    import tiffile as tiffile
except:
    pass

indir = "/Volumes/HDD/Python_laboratory/Vector_visualizer/test"
outdir = "/Volumes/HDD/Python_laboratory/Vector_visualizer/test/out/"
testfile = "/Volumes/HDD/Python_laboratory/Vector_visualizer/test/sub.tif"

def display(array):
    plt.imshow(array, cmap="viridis")
    plt.show()

def polar_bar(theta, radii, color, saveFlag, fname="tst"):


    ax = plt.subplot(111, projection='polar')
    img = plt.scatter(theta, radii, c=color)
    img.set_cmap('viridis')

    #plt.show()
    if saveFlag:
        plt.savefig(outdir+fname+".png", bbox_inches='tight', pad_inches = 0)
        plt.close()
        im = Image.open(outdir+fname+".png")
        #im = im.convert("L")
        im.save(outdir+fname+".jpg")
        im.close()

    if not saveFlag:
        plt.show()


def saveQiver(x,y,u,v, fname, l=None):
    img = plt.quiver(x, y, u, v, l)
    #img.set_cmap('viridis')
    #plt.axis('off')
    #plt.show()
    plt.savefig(outdir+fname+".png", bbox_inches='tight', pad_inches = 0)
    plt.close()
    im = Image.open(outdir+fname+".png")
    im.save(outdir+fname+".jpg")
    im.close()


def saveScatter(arrX, arrY, fname):
    filename =  outdir + fname + ".png"
    img = plt.scatter(arrX, arrY, marker=".", c="black")
    #img.set_cmap('viridis')
    plt.axis('off')
    # plt.show()
    plt.savefig(filename, bbox_inches='tight', pad_inches=0)
    plt.close()
    im = Image.open(filename)
    im = im.convert("L")
    im = im.resize((400,400))
    im.save(filename)
    im.close()

def save2Dhist(x, y, bins, fname):
    plt.hist2d(x, y, bins=bins, norm=LogNorm())
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.savefig(outdir + fname + ".png")
    plt.close()

def arrayProcessor(arr):

    #PIV parameters
    window_size = 8
    overlap = 0
    search_area = 10


    x, y = openpiv.process.get_coordinates(image_size=arr[0].shape, window_size=window_size, overlap=overlap)

    for frame in range(arr.shape[0]):
        if frame == arr.shape[0]-1:
            break

        frame_a = arr[frame]
        frame_b = arr[frame+1]


        u, v, sig2noise = openpiv.process.extended_search_area_piv( frame_a, frame_b, window_size=window_size,
                                                                overlap=overlap,
                                                                dt=1, search_area_size=search_area,
                                                                sig2noise_method='peak2peak' )

        u, v, mask = openpiv.validation.sig2noise_val( u, v, sig2noise, threshold = 1.5)
        u, v = openpiv.filters.replace_outliers(u, v, method='localmean', kernel_size=2)

        size_x = u.shape[-2]
        size_y = u.shape[-1]
        n_pixels = size_x * size_y

        # u, v, mask = openpiv.validation.sig2noise_val(u, v, sig2noise, threshold=noise_threshold)
        vector_lengths = np.sqrt((u ** 2 + v ** 2))
        dot_products = np.zeros(u.shape)
        vector_angles = np.zeros(u.shape)
        dist_out = []
        angles_out = []
        for i in range(n_pixels):  # outer loop for vector a
            a_m = i % size_x
            a_n = i // size_y
            vector_a = np.array((u[a_m][a_n], v[a_m][a_n]))
            distances_to_a = np.sqrt((x - x[a_m][a_n]) ** 2 + (
            y - y[a_m][a_n]) ** 2)  # Maxtrx where each value is its pixel distance to a in the original image
            ab_length_products = vector_lengths[a_m][a_n] * vector_lengths
            for j in range(n_pixels):  # inner loop for vector b
                if j >= i:
                    break  # no need to double count pairs or calculate self-angles

                b_m = j % size_x
                b_n = j // size_y
                vector_b = np.array((u[b_m][b_n], v[b_m][b_n]))

                dot_products[b_m][b_n] = np.dot(vector_a, vector_b)

            vector_angles = np.arccos(dot_products / ab_length_products)
            dist_out.extend(distances_to_a.flatten().tolist())
            angles_out.extend(vector_angles.flatten().tolist())

        save2Dhist(dist_out, angles_out, 30, str(frame))


with tiffile.TiffFile(testfile) as tif:
    arr = tif.asarray()
    arr = arr.astype(np.int32)
    print(arr.shape)

arrayProcessor(arr)

# window_size = 8
# overlap = 0
# search_area = 10
# noise_threshold = 1.5
#
#
# x, y = openpiv.process.get_coordinates(image_size=arr[0].shape, window_size=window_size, overlap=overlap)
#
# frame_a = arr[0]
# frame_b = arr[1]
#
# u, v, sig2noise = openpiv.process.extended_search_area_piv(frame_a, frame_b, window_size=window_size, overlap=overlap,
#                                                            dt=60 * 8, search_area_size=search_area,
#                                                            sig2noise_method='peak2peak')
#
# size_x = u.shape[-2]
# size_y = u.shape[-1]
# n_pixels = size_x * size_y
#
# #u, v, mask = openpiv.validation.sig2noise_val(u, v, sig2noise, threshold=noise_threshold)
# vector_lengths = np.sqrt((u ** 2 + v ** 2))
# dot_products = np.zeros(u.shape)
# vector_angles = np.zeros(u.shape)
# dist_out = []
# angles_out = []
# for i in range(n_pixels): #outer loop for vector a
#     a_m = i % size_x
#     a_n = i // size_y
#     vector_a = np.array( (u[a_m][a_n],v[a_m][a_n]) )
#     distances_to_a = np.sqrt((x - x[a_m][a_n]) ** 2 + (y - y[a_m][a_n]) ** 2) #Maxtrx where each value is its pixel distance to a in the original image
#     ab_length_products = vector_lengths[a_m][a_n] * vector_lengths
#     for j in range(n_pixels): #inner loop for vector b
#         if j >= i:
#             break #no need to double count pairs or calculate self-angles
#
#         b_m = j % size_x
#         b_n = j // size_y
#         vector_b = np.array( (u[b_m][b_n],v[b_m][b_n]) )
#
#         dot_products[b_m][b_n] = np.dot(vector_a, vector_b)
#
#     vector_angles = np.arccos(dot_products / ab_length_products)
#     dist_out.extend(distances_to_a.flatten().tolist())
#     angles_out.extend(vector_angles.flatten().tolist())
#
# plt.hist2d(dist_out, angles_out, bins=40)
# plt.show()

        #polar_bar(vector_angles, vector_lengths, vector_dot_products, False )

#display(vector_angles)
#display(distances)
#plt.scatter(distances, vector_angles)
#plt.show()
#x1, y1 = openpiv.process.get_coordinates(image_size=(20,20), window_size=5, overlap=0)

#display(sig2noise)



#u, v, mask = openpiv.validation.local_median_val(u, v, 1, 1)


#plt.hist(l, 10, (0,0.01))
#plt.scatter(x,y,c=l)

#plt.show()

#u, v = openpiv.filters.replace_outliers( u, v, method='localmean', kernel_size=2)

#x, y, u, v = openpiv.scaling.uniform(x, y, u, v, scaling_factor = 96.52 )

#
#print(vectors[0,2])
#print(vectors.shape)
#
# plt.scatter(u, v)
#plt.show()

##dot = np.dot(x,y)
#>>> x_modulus = np.sqrt((x*x).sum())
#>>> y_modulus = np.sqrt((y*y).sum())
#>>> cos_angle = dot / x_modulus / y_modulus # cosine of angle between x and y
#>>> angle = np.arccos(cos_angle)
#>>> angle
#0.80823378901082499
#>>> angle * 360 / 2 / np.pi # angle in degrees
#46.308384970187326


