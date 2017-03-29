from __future__ import print_function
import openpiv.process
import openpiv.filters
import openpiv.scaling
import openpiv.validation
from matplotlib import pyplot as plt
import numpy as np
import math
from stitch_buddy import *
from PIV_good_code import *
import tiffile as tiffile
import time
import os


def openPIV_array_processor_median(arr, stopFrame, startFrame=0, frameSamplingInterval=1, **piv_params):
    """
    The function first runs a gliding 3-frame temporal median on every pixel to smooth out noise and to remove fast
    moving debree that is not migrating cells.
    Then it does PIV analysis between every n frames in the smoothed time lapse.
    It returns the u and v components of the velocity vectors as two (smaller) numpy arrays.
    Two additional arrays with the x and y coordinates corresponding to the centers of the search windows in the
    original input array are also returned.
    This function should not be run on data that has already been smoothed.

    :param arr:
        (3d numpy array) with a shape of (t, y, x) of type np.int32
    :param stopFrame:
        (int) Last frame to analyze
    :param startFrame:
        (int) First frame to analyze
    :param frameSamplingInterval:
        (int) do PIV between every n frames
    :param piv_params:
        (dict) parameters for the openPIV function extended_search_area_piv
    :return:
        u_component_array, v_component_array, original_x_coord_array, original_y_coord_array

    """
    assert (startFrame < stopFrame) and (stopFrame <= arr.shape[0])

    n_frames = 1 + (stopFrame - startFrame - 4) // frameSamplingInterval

    #original x/y coordinates
    x, y = openpiv.process.get_coordinates(image_size=arr[0].shape,
                                           window_size=piv_params["window_size"],
                                           overlap=piv_params["overlap"])

    #Zero-filled output arrays are created beforehand for maximal performance
    out_u = np.zeros((n_frames, x.shape[0], x.shape[1]))
    out_v = np.zeros_like(out_u)

    for frame in range(startFrame, stopFrame, frameSamplingInterval):

        if frame >= (stopFrame - 4):
            break

        frame_a = np.median(arr[frame:frame + 3], axis=0).astype(np.int32) #median of frames n1,n2,n3
        frame_b = np.median(arr[frame + 1:frame + 4], axis=0).astype(np.int32) #median of frames n2,n3,n4

        #the output arrays are modified in-place with the PIV data
        out_u[frame], out_v[frame] = openpiv.process.extended_search_area_piv(frame_a, frame_b,
                                                                                         window_size=piv_params[
                                                                                             "window_size"],
                                                                                         overlap=piv_params["overlap"],
                                                                                         dt=piv_params["dt"],
                                                                                         search_area_size=piv_params[
                                                                                             "search_area_size"],
                                                                                         sig2noise_method=None
                                                                                             )
    return out_u, out_v, x, y


def openPIV_array_processor(arr, stopFrame, startFrame=0, frameSamplingInterval=1, **piv_params):
    """
    The function does PIV analysis between every n frames in input array.
    It returns the u and v components of the velocity vectors as two (smaller) numpy arrays.
    Two additional arrays with the x and y coordinates corresponding to the centers of the search windows in the
    original input array are also returned.

    This function should be run on data that has already been smoothed, or when smoothing is undesirable.

    :param arr:
        (3D numpy array) with a shape of (t, y, x) of type np.int32
    :param stopFrame:
        (int) Last frame to analyze
    :param startFrame:
        (int) First frame to analyze
    :param frameSamplingInterval:
        (int) do PIV between every n frames
    :param piv_params:
        (dict) parameters for the openPIV function extended_search_area_piv
    :return:
        u_component_array, v_component_array, original_x_coord_array, original_y_coord_array

    """
    assert (startFrame < stopFrame) and (stopFrame <= arr.shape[0])

    n_frames = (stopFrame - startFrame - 1) // frameSamplingInterval

    # original x/y coordinates
    x, y = openpiv.process.get_coordinates(image_size=arr[0].shape,
                                           window_size=piv_params["window_size"],
                                           overlap=piv_params["overlap"])

    # Zero-filled output arrays are created beforehand for maximal performance
    out_u = np.zeros((n_frames, x.shape[0], x.shape[1]))
    out_v = np.zeros_like(out_u)

    #pairs the original frame number with the indexing in the output, since frame n in the input might not be n in output
    for frame, i in zip(range(startFrame, stopFrame, frameSamplingInterval), range(n_frames)):

        if frame >= (stopFrame - 1):
            break

        frame_a = arr[frame]
        frame_b = arr[frame + 1]

        out_u[i], out_v[i] = openpiv.process.extended_search_area_piv(frame_a, frame_b,
                                                                      window_size=piv_params["window_size"],
                                                                      overlap=piv_params["overlap"],
                                                                      dt=piv_params["dt"],
                                                                      search_area_size=piv_params["search_area_size"],
                                                                      sig2noise_method=None)

    return out_u, out_v, x, y


def alignment_index(u, v, alsoReturnMagnitudes=False):

    """
    Returns an array of the same shape as u and v with the alignment index (ai), defined as in Malinverno et. al 2017.
    For every frame the ai is the average of the dot products of the mean velocity vector with each individual
    vector, all divided by the product of their magnitudes.

    If alsoReturnMagnitudes is set to True, then an additional array with the vector magnitudes, i.e, speeds in
    pixels/frame is also returned.

    :param u:
        2D numpy array with u component of velocity vectors
    :param v:
        2D numpy array with v component of velocity vectors
    :param alsoReturnMagnitudes:
        (bool) Should the function also return the vector magnitudes
    :return:
        nunpy array with size=input.size where every entry is the alignment index in that pixel

    """
    assert (u.shape == v.shape) and (len(u.shape) == 2)  # Only single frames are processed

    vector_0 = np.array((np.mean(u), np.mean(v)))
    v0_magnitude = np.linalg.norm(vector_0)

    vector_magnitudes = np.sqrt((np.square(u) + np.square(v))) #a^2 + b^2 = c^2
    magnitude_products = vector_magnitudes * v0_magnitude
    dot_products = u * vector_0[0] + v * vector_0[1] #Scalar multiplication followed by array addition

    ai = np.divide(dot_products, magnitude_products)

    if alsoReturnMagnitudes:
        return ai, vector_magnitudes

    else:
        return ai


def msv(u, v): #Mean Square Velocity
    """
     Calculates the mean square velocity of one frame from PIV data

    :param u:
        2D numpy array with the u component of velocity vectors
    :param v:
        2D numpy array with the v component of velocity vectors
    :return:
        (float) the mean square velocity of the velocity vectors in the frame
    """

    msv = np.mean(np.square(u)+np.square(v))

    return msv

def rms(u, v): #Root Mean Square Velocity
    """
    Calculates the root mean square velocity of one frame from PIV data, i.e. speed or vector magnitudes in the unit
    pixels/frame. This is equivalent to taking the square root of the mean square velocity.

    :param u:
        2D numpy array with the u component of velocity vectors
    :param v:
        2D numpy array with the v component of velocity vectors
    :return:
        (float) the root mean square velocity of the velocity vectors in the frame
    """
    rms = np.sqrt(np.mean(np.square(u)+np.square(v))) #sqrt(u^2+v^2)

    return rms

def smvvm(u, v):  # square_mean_vectorial_velocity_magnitude
    """
    Array addition of the squared average vector components, used in calculating the instantaneous order parameter

    :param u:
        2D numpy array with the u component of velocity vectors
    :param v:
        2D numpy array with the u component of velocity vectors
    :return:
        2D numpy array with the u component of velocity vectors

    """

    return np.square(np.mean(u)) + np.square(np.mean(v))

def instantaneous_order_parameter(u, v):
    """
    Calculates the instantaneous order parameter (iop) in one PIV frame see  Malinverno et. al 2017 for a more detailed
    explanation. The iop is a measure of how similar the vectors in a field are, which takes in to account both the
    direcions and magnitudes of the vectors. iop always between 0 and 1, with iop = 1 being a perfectly uniform field
    of identical vectors, and iop = 0 for a perfectly random field.

    :param u:
        2D numpy array with the u component of velocity vectors
    :param v:
        2D numpy array with the u component of velocity vectors
    :return:
        (float) iop of vector field
    """
    return smvvm(u, v) / msv(u, v) #square_mean_vectorial_velocity_magnitude/Mean Square Velocity


def get_v0_plus_r_coordinates_cardinal(array_shape, v0_cord, r):

    """
    Gets the 4 coordinates of the pixels that are r pxels away from the coordinate (v0_x, v0_y) along the four
    cardinal directions. Coordinates are returned in the order: up/top, right, down/bottom, left.

    The coordinates are returned in Matrix/numpy form (row,col), i.e. (y,x) when compared to traditional image
    coordinate numbering. The input coordinate for v0 is also expected to be in this format (y_coord, x_coord).

    If a coordinate falls outside of the image it is not returned.

    :param array_shape:
        (tuple) shape of 2D numpy array
    :param v0_cord:
        (tuple) of (ints), (row, col) coordinates of v0 in matrix notation.
    :param r:
        (int) distance from v0 along the cardinal directions
    :return:
        (list) of (tuple) coordinates in martix notation, if no valid coordinates are found an empty list is returned.

    """

    array_width = array_shape[1] #number of columns
    array_height = array_shape[0] #number of rows
    v0_r = v0_cord[0] #row numer of v0
    v0_c = v0_cord[1] #column number of v0

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
    Gets the vector v0 at coordinate (v0_coord) from a velocity vector field. Grows the distance r from r_min to r_max
    in incremets of r_step. For each distance r calculates the average of the (cos) angles between v0 and v0+r in the four
    cardinal directions. The result is stored in a dictionary where the keys are distances and the values are lists of
    average angles for that distance. The updated resultsDict is returned.


    :param u_array:
        2D numpy array with the u component of velocity vectors
    :param v_array:
        2D numpy array with the v component of velocity vectors
    :param v0_coord:
        (tuple) Coordinates of vector_0
    :param resultsDict:
        (dict)
    :param r_max:
        (int) maximum distance in pixels from v0 to calculate the angels for
    :param r_step:
        (int) step size to grow the distance r by
    :param r_min:
        (int) starting distance from v0 to calculate the angels for
    :return:
        (dict) input resultDict updated with key = radius, value = list of means for cos(theta) v_0-v_r
        {radius:[list of mean angles]}
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
            if magnitudes[c] == 0: #avoid masked or erronious values
                pass
            else:
                c_vv = dot_products[c] / magnitudes[c]
                resultsDict[r].append(c_vv)

    for k, v in resultsDict.items():
        if len(resultsDict[k]) == 0:
            resultsDict.pop(r, None)  # No need to save empty data lists, it breaks the statistics

    return resultsDict


def do_it_all(indir, fname, outdir,
              stopFrame=200,
              px_resolution = 5.466,
              time_resolution = 5,
              r_max=300,
              r_step=1,
              n_sigma = 5,
              intervalWidth=10,
              saveFigFlag = False):
    """

    :param indir:
        Input directory
    :param fname:
        Name of input file
    :param outdir:
        Where to save graphs
    :param stopFrame:
        Last frame to
    :param px_resolution:
        Pixel resolution of input images in um/pixel
    :param time_resolution:
        Time resoluion of input images in frames/hour
    :param r_max:
        maximum distance in pixels to calculate angles for
    :param r_step:
        (int) increment size of distance r during angle calculations
    :param n_sigma:
        significance level that determines the correlation length
    :param intervalWidth:
        (int) width in frames to integrate angle data for, i.e. temporal itegration
    :return:
        (dict) output data for further processing
    """
    # PIV parameters
    piv_params = dict(window_size=32,
                      overlap=16,
                      dt=1,
                      search_area_size=36,
                      sig2noise_method="peak2peak")

    px_scale = px_resolution * piv_params["overlap"]
    piv_scaler =  px_resolution*time_resolution # um/px * frames/hour = um*frames/px*hours

    print("Working on: %s" % (fname))
    t0 = time.time()
    with tiffile.TiffFile(indir+fname) as tif:
        arr = tif.asarray()

    arr = arr.astype(np.int32)

    t1 = time.time()
    print("It took %.2f s to load %s, and it has shape: %s, now processing PIV..." % (t1-t0, fname, arr.shape))

    arr_u, arr_v, arr_x, arr_y = openPIV_array_processor(arr, stopFrame=stopFrame, **piv_params)
    t2 = time.time()
    print("It took %.2f s to process PIV, and output arrays have the shape: %s" % (t2-t1, arr_u.shape))

    inst_order_params, align_idxs, speeds, timepoints = [], [], [], []

    for frame in range(arr_u.shape[0]):

        iop = instantaneous_order_parameter(arr_u[frame], arr_v[frame])
        inst_order_params.append(iop)

        ai = alignment_index(arr_u[frame], arr_v[frame], alsoReturnMagnitudes=True)
        aidx = np.nanmean(ai[0])
        align_idxs.append(aidx)
        speed = np.nanmean(ai[1])*piv_scaler #px/frame * um*frames/px*hours -> um/hour
        speeds.append(speed)

        timepoint = frame*(1/float(time_resolution))
        timepoints.append(timepoint)

        #print(frame, timepoint, aidx, iop, speed)

    if saveFigFlag:

        plt.plot(timepoints, speeds, 'r', label=fname[:7])
        plt.xlabel("Time (h)")
        plt.ylabel("Mean speed (um/h)")
        plt.title("Velocity vector magnitudes (speed)")
        plt.legend()
        plt.savefig(outdir+fname[:-4]+"_speed.pdf", bbox_inches='tight', pad_inches=0)
        #plt.savefig(outdir + fname + "_speed.png", bbox_inches='tight', pad_inches=0)
        plt.close()

        plt.plot(timepoints, align_idxs, 'b', label=fname[:7])
        plt.xlabel("Time (h)")
        plt.ylabel("Align index")
        plt.title("Alignment index")
        plt.legend()
        plt.savefig(outdir+fname[:-4] + "_alignidx.pdf", bbox_inches='tight', pad_inches=0)
        #plt.savefig(outdir + fname + "_alignidx.png", bbox_inches='tight', pad_inches=0)
        plt.close()

        plt.plot(timepoints, inst_order_params, 'g', label=fname[:7])
        plt.xlabel("Time (h)")
        plt.ylabel("$\psi$")
        plt.title("Instantaneous order parameter ($\psi$)")
        plt.legend()
        plt.savefig(outdir+fname[:-4]+ "_iop.pdf", bbox_inches='tight', pad_inches=0)
        #plt.savefig(outdir + fname + "_iop.png", bbox_inches='tight', pad_inches=0)
        plt.close()

    t3=time.time()
    print("It took %.2f s to do Alignemt indexes, order params, and speeds" % (t3-t2) )

    metaResults = {}

    #temporal integration
    for interval in range(0, arr_u.shape[0], intervalWidth):
        results={}
        tmp_u = arr_u[interval:interval + intervalWidth]
        tmp_v = arr_v[interval:interval + intervalWidth]

        for t in range(tmp_u.shape[0]):
            for diagonal in range(0, min(tmp_u.shape[1], tmp_u.shape[2])):
                results = get_all_angles(tmp_u[t], tmp_v[t], (diagonal,diagonal), results, r_max, r_step)
        metaResults[interval] = dict(results)


    lcorrs = {} #interval:correlation length in

    for interv in sorted(metaResults.keys()):
        r = []
        avg_angle = []
        results = metaResults[interv]

        for radius, angles_list in results.items():

            if (len(angles_list) == 0):
                print("Empty value at interal %i and r=%i" % (interv, radius))
                break

            sanitized_angles = [a for a in angles_list if a <= 1.0]
            if len(sanitized_angles) != len(angles_list):
                print("Bad angles at interval %i and radius %i, number ok: %i" % (interv, radius, len(sanitized_angles)))
            mean = np.nanmean(sanitized_angles)
            mean_deg = math.acos(mean) * (180 / math.pi)
            sd = np.nanstd(sanitized_angles)
            sd_deg = math.acos(sd) * (180 / math.pi)
            SEM = sd_deg / math.sqrt(len(angles_list))
            # print(k, mean_deg, len(v), mean_deg+3*SEM)
            if (mean_deg + n_sigma * SEM > 90) and (interv not in lcorrs) and (len(r) != 0):
                #TODO interval_center = interval + intervalWidth/2 * time_resolution
                lcorrs[interv] = r[-1]
                # print("%i-sigma reached at r=%i on interval %i, last significant distance was %.2f um" % (n_sigma, radius, interv, r[-1]))

            r.append(radius * px_scale)
            avg_angle.append(mean_deg)

        label = "interval: " + str(((interv) / time_resolution)) + " - " + \
                str((interv + intervalWidth) / time_resolution) + " h, Lcorr = " + \
                str(int(lcorrs[interv])) + " um"

        plt.plot(r, avg_angle, label=label)


    plt.legend()
    plt.title("Average angle between velocity vectors " + fname[:7])
    plt.xlabel("Distance in um")
    plt.ylabel("Mean angle (degrees)")

    if saveFigFlag:
        plt.savefig(outdir+fname+".pdf")

    print("All done in %.2f s" % (time.time()-t0))
    plt.show()

    return (inst_order_params, align_idxs, speeds, timepoints, lcorrs)

indir = "/Users/jens_e/Python_laboratory/Vector_visualizer/Data for analysis script/"
outdir = "/Users/jens_e/Python_laboratory/Vector_visualizer/test/out/"
analysisfiles = ["test.tif"]
#testfiles = os.listdir(indir)

for fname in analysisfiles:
    t = do_it_all(indir, fname, outdir,
              stopFrame=20,
              px_resolution = 5.466,
              time_resolution = 5,
              r_max=30,
              r_step=1,
              n_sigma = 5,
              intervalWidth=2,
              saveFigFlag=False)

print(t)