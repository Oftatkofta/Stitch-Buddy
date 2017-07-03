from __future__ import print_function
import pickle
import os
from matplotlib import pyplot as plt
import numpy as np
import csv

indir = "/Users/jens_e/Python_laboratory/Vector_visualizer/Data for analysis script/out_diagonal/"
filez = [f for f in os.listdir(indir) if ".pickle" in f]

out={}
for fname in filez:
    with open(indir+fname) as f:
        out[fname] = pickle.load(f)

    inst_order_params, align_idxs, speeds, timepoints, lcorrs = out[fname]
    timepoints = np.array(timepoints) + 1
    correlation_lengths, times = [], []
    for k in sorted(lcorrs.keys()):
        correlation_lengths.append(lcorrs[k])
        times.append((k/5.0)+1)
    times.append(40.8)
    correlation_lengths.append(None)
    plt.plot(times, correlation_lengths, label=fname[:6])

    with open(indir + fname[:6] + "_t_speed_ai_iop_cl.csv", 'wb') as csv_out:
        csvwriter = csv.writer(csv_out)
        csvwriter.writerow(["time_h", "speed_um_h", "alignment_index", "inst_order_param", "correlaton_length_5_sigma_um"])
        for t, s, ai, iop in zip(timepoints, speeds, align_idxs, inst_order_params):
            print([round(t, 2), round(s, 2), ai, iop, lcorrs.get(t, None)])
            csvwriter.writerow([round(t,2),round(s,2), ai, iop, lcorrs.get(t, None)])

    with open(indir + fname[:6] + "_t_speed_ai_iop.csv", 'wb') as csv_out:
        csvwriter = csv.writer(csv_out)
        csvwriter.writerow(["time_h", "speed_um_h", "alignment_index", "inst_order_param"])
        for t, s, ai, iop in zip(timepoints, speeds, align_idxs, inst_order_params):
            print([round(t, 2), round(s, 2), ai, iop])
            csvwriter.writerow([round(t,2),round(s,2), ai, iop])

    with open(indir + fname[:6] + "_correlation_lengths.csv", 'wb') as csv_out:
         csvwriter = csv.writer(csv_out)
         csvwriter.writerow(["time_h", "correlaton_length_5_sigma_um"])
         for t, lcorr in zip(times, correlation_lengths):
             csvwriter.writerow([round(t,2),lcorr])


# plt.title("5-$\sigma$ correlation length of velocity angle")
# plt.xlabel("time (h)")
# plt.ylabel("distance ($\mu$m)")
# plt.legend()
# plt.savefig(indir+"all_correlation_lengths.pdf", bbox_inches='tight', pad_inches=0)
# plt.close()
#
# for fname, v in out.items():
#     plt.plot(timepoints, v[0], label=fname[:6]) #iop
#
# plt.xlabel("time (h)")
# plt.ylabel("$\psi$")
# plt.title("Instantaneous order parameter ($\psi$)")
# plt.legend()
# savename = indir + "all_instantaneous_order_parameters.pdf"
# plt.savefig(savename, bbox_inches='tight', pad_inches=0)
# plt.close()
#
# for fname, v in out.items():
#     plt.plot(timepoints, v[1], label=fname[:6]) #ai
# plt.xlabel("time (h)")
# plt.ylabel("Alignment index")
# plt.title("Alignment index")
# plt.legend()
# savename = indir + "all_alignment_index.pdf"
# plt.savefig(savename, bbox_inches='tight', pad_inches=0)
# plt.close()
#
# for fname, v in out.items():
#     plt.plot(timepoints, v[2], label=fname[:6]) #speed
# plt.xlabel("time (h)")
# plt.ylabel("mean speed ($\mu$m/h)")
# plt.title("Velocity vector magnitudes (speed)")
# plt.legend()
# savename = indir +"all_speeds.pdf"
# plt.savefig(savename, bbox_inches='tight', pad_inches=0)
# plt.close()
