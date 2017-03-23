rom __future__ import print_function
import openpiv.process
import openpiv.filters
import openpiv.scaling
import openpiv.validation
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import animation, rc
from IPython.display import HTML
import numpy as np
import math
from PIL import Image
import tiffile as tiffile
import cv2
import time

#PIV parameters
piv_params = dict(window_size = 32,
                overlap = 16,
                dt = 1,
                search_area_size = 36,
                sig2noise_method = "peak2peak")

def openPIV_array_processor_median(arr, stopFrame, startFrame=0, frameSamplingInterval=1, **piv_params):
    assert (startFrame < stopFrame) and (stopFrame <= arr.shape[0])

    n_frames = 1 + (stopFrame - startFrame - 4) // frameSamplingInterval

    x, y = openpiv.process.get_coordinates(image_size=arr[0].shape,
                                           window_size=piv_params["window_size"],
                                           overlap=piv_params["overlap"])

    out_u = np.zeros((n_frames, x.shape[0], x.shape[1]))
    out_v = np.zeros_like(out_u)

    for frame in range(startFrame, stopFrame, frameSamplingInterval):

        if frame >= (stopFrame - 4):
            break

        frame_a = np.median(arr[frame:frame + 3], axis=0).astype(np.int32)
        frame_b = np.median(arr[frame + 1:frame + 4], axis=0).astype(np.int32)

        out_u[frame], out_v[frame], sig2noise = openpiv.process.extended_search_area_piv(frame_a, frame_b,
                                                                                         window_size=piv_params[
                                                                                             "window_size"],
                                                                                         overlap=piv_params["overlap"],
                                                                                         dt=piv_params["dt"],
                                                                                         search_area_size=piv_params[
                                                                                             "search_area_size"],
                                                                                         sig2noise_method=piv_params[
                                                                                             "sig2noise_method"])
    return out_u, out_v, x, y

def alignment_index(u, v, alsoReturnMagnitudes=False):

    """
    Returns an array of the same shape as u and v with the aligmnent index, as in Malinverno et. al 2017
    if returnMagnitudes is set to True, then an additional array with the vector magnitudes is also returned.
    """

    assert (u.shape == v.shape) and (len(u.shape) == 2)  # Only single frames are processed

    vector_0 = np.array((np.mean(u), np.mean(v)))
    v0_magnitude = np.linalg.norm(vector_0)

    vector_magnitudes = np.sqrt((np.square(u) + np.square(v)))
    magnitude_products = vector_magnitudes * v0_magnitude
    dot_products = u * vector_0[0] + v * vector_0[1]

    ai = np.divide(dot_products, magnitude_products)

    if alsoReturnMagnitudes:
        return ai, vector_magnitudes

    else:
        return ai


def msv(u, v): #Mean Square Velocity
    msv = np.mean(np.square(u)+np.square(v))
    return msv

def rms(u, v): #Root Mean Square Velocity
    rms = np.sqrt(np.mean(np.square(u)+np.square(v)))
    return rms

def smvvm(u, v):  # square_mean_vectorial_velocity_magnitude

    return np.square(np.mean(u)) + np.square(np.mean(v))

def instantaneous_order_parameter(u, v):
    return smvvm(u, v) / msv(u, v)



def get_v0_plus_r_coordinates_cardinal(array_shape, v0_cord, r):

    """
    Gets the 4 coordinates of the pixels r away from the coordinate (v0_x, v0_y) in the four cardinal directions.
    The coordinates are returned in Matrix/numpy form (row,col), i.e. (y,x) when compared to traditional image
    coordinate numbering.

    :param array_shape: (tuple)
    :param v0_cord: (tuple) of (ints), (row, col) coordinates of v0 in matrix notation.
    :param r: (int) distance from v0 along the cardinal directions
    :return: List of coordinates in martix notation, if no valid coordinates are found an empty list is returned.

    """

    array_width = array_shape[1]
    array_height = array_shape[0]
    v0_r = v0_cord[0]
    v0_c = v0_cord[1]

    assert r>0, "r needs to be positive and >0 !"
    assert array_width > v0_c, "v0_y needs to be less than array_width!"
    assert v0_c >= 0 , "v0_y needs to be positive!"
    assert array_height > v0_r, "v0_y needs to be < array_height!"
    assert v0_r >= 0, "v0_y needs to be positive and < array_height!"

    top_r = v0_r-r
    right_c = v0_c+r
    bottom_r = v0_r+r
    left_c = v0_c-r

    out = []


    if (top_r >= 0):
        out.append((top_r, v0_c))

    if (right_c < array_width):
        out.append((v0_r, right_c))

    if (bottom_r < array_height):
        out.append((bottom_r, v0_c))

    if (left_c >= 0):
        out.append((v0_r, left_c))

    return out


def get_all_angles(u_array, v_array, v0_coord, resultsDict, r_max, r_step=1, r_min=1):
    """

    :param u_array:
    :param v_array:
    :param v0_coord:
    :param resultsDict:
    :param r_max:
    :param r_step:
    :param r_min:
    :return: updated resultsDict with key:radius val:list of means for cos(theta) v_0-v_r
    """
    assert u_array.shape == v_array.shape, "u and v component arrays have to have identical shapes"

    v0_u = u_array[v0_coord]
    v0_v = v_array[v0_coord]

    magnitudes = np.sqrt(np.square(u_array) + np.square(v_array))
    v0_magnitude = magnitudes[v0_coord]

    dot_products = u_array * v0_u + v_array * v0_v  # Computes ALL the dot products with v0
    magnitudes = magnitudes * v0_magnitude  # Multiplies all magnitudes by the magnitude of v0

    for r in range(r_min, r_max, r_step):
        if not r in resultsDict:
            resultsDict[r] = []

        coords = get_v0_plus_r_coordinates_cardinal(u_array.shape, v0_coord, r)

        if len(coords) == 0:
            break  # stop when we run out of valid coordinates
        for c in coords:
            if magnitudes[c] == 0:
                pass
            else:
                c_vv = dot_products[c] / magnitudes[c]
                resultsDict[r].append(c_vv)

    for k, v in resultsDict.items():
        if len(resultsDict[k]) == 0:
            resultsDict.pop(r, None)  # No need to save empty data lists, it breaks the statistics

    return resultsDict


tmp = openPIV_u2.copy()

results = {}

for t in range(tmp.shape[0]):
    for d in range(0, min(tmp.shape[1], tmp.shape[2])):
        results = get_all_angles(openPIV_u2[t], openPIV_v2[t], (d, d), results, 300, 1)

x = []
y = []

for k, v in results.items():

    if (len(v) == 0):
        print("Empty value at r=%i" % (k))
        break
    mean = np.nanmean(v)
    mean_deg = math.acos(mean) * (180 / math.pi)
    sd = np.nanstd(v)
    sd_deg = math.acos(sd) * (180 / math.pi)
    SEM = sd_deg / math.sqrt(len(v))
    # print(k, mean_deg, len(v), mean_deg+3*SEM)
    if (mean_deg + 3 * SEM > 90):
        print("3-sigma reached at r=%i, last significant distance was %.2f um" % (k, x[-1]))
        break
    x.append(k * 41.8624)
    y.append(mean_deg)

plt.plot(x, y, 'r', label="aCos(v(o)v(r))")
plt.legend()
plt.title("Average angle between velocity vectors")
plt.xlabel("Distance in um")
plt.ylabel("Mean angle (degrees)")
plt.show()
