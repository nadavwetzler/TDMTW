
import os
from obspy import Stream, Trace, read
import numpy as np
import datetime
from GF_functions import *
path2folder = os.path.expanduser('~/Dropbox/Moment-tensor/')
path2folder = os.path.expanduser('~/Documents/')
Path2_bin = os.path.expanduser('%~/Dropbox/Moment-tensor/BIN_HighSierra/')# MAC Mojave
#--------------------Set PARAMETERS --------------------
# DEPTH = np.array([5])
# DIST0 = np.array([255])
DEPTH = np.arange(0, 50, 1)  # SET DEPTH RANGE
DIST0 = np.arange(1, 100, 1)  # SET DISTANCES FOR GREEN FUNCTION
dH = 0.05 # Thickness of artificial layer in case Z == layer interface
f_low = 0.0
f_high = 0.0
npts0 = 4096
dt0 = 0.1
ReductionVel = 100



#-------GIL7---------------------------------------
# MT_F = "GIL7_d"
# MODEL = np.zeros((7, 6))
# #                Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [1.000, 3.20, 1.50, 2.28, 600, 300]  # 1.76
# MODEL[1, 0:6] = [3.000, 4.50, 2.40, 2.28, 600, 300]  # 1.79
# MODEL[2, 0:6] = [4.00, 4.80, 2.78, 2.58, 600, 300]  # 1.77
# MODEL[3, 0:6] = [5.00, 5.51, 3.18, 2.58, 600, 300]  # 1.77
# MODEL[4, 0:6] = [17.00, 6.21, 3.40, 2.68, 600, 300]  # 1.73
# MODEL[5, 0:6] = [25.0, 6.89, 3.98, 3.00, 600, 300]  # 1.78
# MODEL[6, 0:6] = [500.0, 7.83, 4.52, 3.26, 600, 300]  # 1.78


#-------Gitterman 2002 Velocity model of Israel ---------------------------------------
MT_F = "Git02_d"
MODEL = np.zeros((7, 6))
#               Depth   Vp   Vs   Dens  Qa   Qb       PR
MODEL[0, 0:6] = [2.010, 4.78, 2.72, 2.10, 600, 300]  # 1.76
MODEL[1, 0:6] = [9.760, 5.72, 3.20, 2.68, 600, 300]  # 1.79
MODEL[2, 0:6] = [16.01, 6.27, 3.54, 2.77, 600, 300]  # 1.77
MODEL[3, 0:6] = [28.01, 6.53, 3.69, 2.85, 600, 300]  # 1.77
MODEL[4, 0:6] = [32.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
MODEL[5, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
MODEL[6, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

#------------- RUN ------------------------
DIST2 = mk_DIST(DIST0)
min_dist = np.min(DIST2[DIST2 > 0])
max_dist = np.max(DIST2)
Path2_Green = set_path_gf(MT_F, path2folder, dt0, f_low, f_high, min_dist, max_dist)
MODELN = np.zeros((MODEL.shape[0] + 1, 6))
mk_b2s(Path2_Green,npts0, dt0)
for kk in range(np.size(DIST2,0)):
    DIST = DIST2[kk]
    mk_perl1(Path2_Green, dt0, npts0, f_high, f_low)

    for ii in range(len(DEPTH)):
        mk_gf1(MODEL, MODELN, DEPTH[ii], npts0,dt0, Path2_bin, Path2_Green, MT_F, DIST, ReductionVel)
        GFlist = GetGreenFiles(Path2_Green, DEPTH[ii])
        GRN1 = Stream()
        mseed_dir = Path2_Green+'MSEED/'
        if not os.path.exists(mseed_dir):
            os.mkdir(mseed_dir)
        for jj in range(len(GFlist)):
            [GRN0, nameG, depth, dist, pathFF] = Green2MSD(Path2_Green + GFlist[jj])
            GRN0.write(mseed_dir + "%s.MSEED" % GFlist[jj], format="MSEED")
            # GRN1 += GRN0
#         GRN1.write(mseed_dir + "S%s.D%s.MSEED" % (kk, DEPTH[ii]), format="MSEED")
# for ii in range(len(DEPTH)):
#     GRD = read(mseed_dir + "S*.D%s.MSEED" % DEPTH[ii])
#     GRD.write(mseed_dir + "SF.D%s.MSEED" % DEPTH[ii], format="MSEED")











