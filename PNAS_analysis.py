from __future__ import print_function
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
from stitch_buddy import *
from PIV_good_code import *
import tiffile as tiffile
import cv2
import time



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


def openPIV_array_processor(arr, stopFrame, startFrame=0, frameSamplingInterval=1, **piv_params):
    assert (startFrame < stopFrame) and (stopFrame <= arr.shape[0])

    n_frames = (stopFrame - startFrame - 1) // frameSamplingInterval

    x, y = openpiv.process.get_coordinates(image_size=arr[0].shape,
                                           window_size=piv_params["window_size"],
                                           overlap=piv_params["overlap"])

    out_u = np.zeros((n_frames, x.shape[0], x.shape[1]))
    out_v = np.zeros_like(out_u)

    for frame, i in zip(range(startFrame, stopFrame, frameSamplingInterval), range(n_frames)):

        if frame >= (stopFrame - 1):
            break
        print("PIV analysis on frame %i." %(frame))
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
            resultsDict.pop(k, None)  # No need to save empty data lists, it breaks the statistics

    return resultsDict


def do_it_all(indir, fname, outdir, maxFrameToAnalyze = 201, intervalSize = 10, saveplots = True ):
    """

    Args:
        indir:
        fname:
        outdir:
        maxFrameToAnalyze: run analysis for this long
        intervalSize: No of frames to integrate for the correlation length analysis
        saveplots:

    Returns:

    """
    # PIV parameters
    piv_params = dict(window_size=32,
                      overlap=16,
                      dt=1,
                      search_area_size=36,
                      sig2noise_method="peak2peak")

    px_resolution = 5.466 # um/px
    time_resolution = 5 # frames/h
    r_max = 300 #piv-pixels = overlap*px_resolution um
    piv_scaler =  px_resolution*time_resolution # piv data in px/frame
    n_sigma = 3
    px_scale = px_resolution * piv_params["overlap"]
    print("Working on: %s" % (fname))
    t0 = time.time()

    with tiffile.TiffFile(indir+fname) as tif:
        arr = tif.asarray()

    arr = arr.astype(np.int32)
    arr = arr[0:maxFrameToAnalyze]
    t1 = time.time() - t0
    print("It took %.2f s to load %s, and has shape: %s" % (t1, fname, arr.shape))
    t2 = time.time() - (t0+t1)
    arr_u, arr_v, arr_x, arr_y = openPIV_array_processor(arr, stopFrame=arr.shape[0], **piv_params)

    print("It took %.2f s to process PIV, and has shape: %s" % (t2, arr_u.shape))
    t3 = time.time()-t0+t1+t2

    inst_order_params, align_idxs, speeds, timepoints = [], [], [], []

    for frame in range(arr_u.shape[0]):
        iop = instantaneous_order_parameter(arr_u[frame], arr_v[frame])
        inst_order_params.append(iop)
        ai = alignment_index(arr_u[frame], arr_v[frame], alsoReturnMagnitudes=True)
        aidx = np.nanmean(ai[0])
        align_idxs.append(aidx)
        speed = np.nanmean(ai[1])*piv_scaler
        speeds.append(speed)
        timepoint = frame*(1/float(time_resolution))
        timepoints.append(timepoint)
        print(frame, timepoint, aidx, iop, speed)

    if saveplots:
        plt.plot(timepoints, speeds, 'r', label=fname[:6])
        plt.xlabel("Time (h)")
        plt.ylabel("Mean speed (um/h)")
        plt.title("Velocity vector magnitudes (speed)")
        plt.legend()
        savename = outdir+fname[:-4]+"_speeds.pdf"
        print(savename)
        plt.savefig(savename)
        #plt.savefig(outdir + fname + "_speed.png", bbox_inches='tight', pad_inches=0)
        plt.close()

        plt.plot(timepoints, align_idxs, 'b', label=fname[:6])
        plt.xlabel("Time (h)")
        plt.ylabel("Align index")
        plt.title("Alignment index")
        plt.legend()
        plt.savefig(outdir+fname[:-4]+"_alignidx.pdf", bbox_inches='tight', pad_inches=0)
        #plt.savefig(outdir + fname + "_alignidx.png", bbox_inches='tight', pad_inches=0)
        plt.close()

        plt.plot(timepoints, inst_order_params, 'g', label=fname[:6])
        plt.xlabel("Time (h)")
        plt.ylabel("$\psi$")
        plt.title("Instantaneous order parameter ($\psi$)")
        plt.legend()
        plt.savefig(outdir+fname[:-4]+ "_iop.pdf", bbox_inches='tight', pad_inches=0)
        #plt.savefig(outdir + fname + "_iop.png", bbox_inches='tight', pad_inches=0)
        plt.close()

    print("It took %.2f s to do Alignemt indexes, order params and speeds" % (t3))
    metaResults = {}

    for interval in range(0, arr_u.shape[0], intervalSize):
        results={}
        tmp_u = arr_u[interval:interval+10]
        tmp_v = arr_v[interval:interval + 10]
        for t in range(tmp_u.shape[0]):
            for d in range(0, min(tmp_u.shape[1], tmp_u.shape[2])): #follow the diagonal
                results = get_all_angles(tmp_u[t], tmp_v[t], (d,d), results, r_max, 1)
        metaResults[interval] = dict(results)

    lastrs = [["n_sigma", "radius_px", "max_sign_r_um"]]
    for interv in sorted(metaResults.keys()):
        x = []
        y = []
        results = metaResults[interv]

        for k, v in results.items():

            if (len(v) == 0):
                print("Empty value at interal %i and r=%i" % (interv, k))
                break
            mean = np.nanmean(v)
            try:
                mean_deg = math.acos(mean) * (180 / math.pi)

            except:
                print("Domain error in nterval: %i" % (k))
                break
            sd = np.nanstd(v)
            sd_deg = math.acos(sd) * (180 / math.pi)
            SEM = sd_deg / math.sqrt(len(v))
            # print(k, mean_deg, len(v), mean_deg+3*SEM)
            x.append(k * px_scale)
            y.append(mean_deg)
            if (mean_deg + n_sigma * SEM > 90):
                lastrs.append([n_sigma, k, x[-1]])
                print("%i-sigma reached at r=%i on interval %i, last significant distance was %.2f um" % (
                n_sigma, k, interv, x[-1]))

        lab = "interval: " + str(((interv) / 5.0)) + " - " + str((interv + intervalSize) / 5.0) + " h "
        plt.plot(x, y, label=lab)

    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand")
    plt.title("Average angle between velocity vectors " + fname[:7])
    plt.xlabel("Distance in um")
    plt.ylabel("Mean angle (degrees)")
    plt.savefig(outdir+fname[:-4]+"_cvv.pdf", bbox_inches='tight')
    plt.close()
    print("All done in %.2f s" % (time.time()-t0))
    return lastrs

#indir = "/Users/jens_e/Python_laboratory/Vector_visualizer/test/"
#outdir = "/Users/jens_e/Python_laboratory/Vector_visualizer/test/out/"
#testfile = "sub_gliding_median_3_frames.tif"

outdir = r"Q:\\PNAS revision data\\Analysis\\"
indir= r"Q:\\PNAS revision data\\"
files = ["160124_H2B_serum_3ml_T0_1h_Median_3_frames_RESCALE.tif",
         "160126_H2B_T0_1h_3ml_serum_Median_3_frames_RESCALE.tif",
         "160205_H2B_3ml_T0_1h_RESCALE_Median_3_frames.tif"]

#l=do_it_all(indir, "test.tif", outdir, maxFrameToAnalyze = 201, intervalSize = 2, saveplots = True)
#print(l)
for f in files:
    thefile = open(f[:7]+'rvals.txt', 'w')
    lcorrs = do_it_all(indir, f, outdir, maxFrameToAnalyze = 201, intervalSize = 10, saveplots = True)
    for item in lcorrs:
        thefile.write("%s\n" % item)
        thefile.close()
