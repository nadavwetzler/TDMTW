# TDMTW GUI
# !/bin/env python
# /**********************************************************************************
# *    Copyright (C) by Nadav Wetzler                                               *
# *                                                                                 *
# *    Automatic Generotor of Green's Functions                                     *
# *                                                                                 *
# *    Requairments - the software used a pre-compiled code                         *
# *    Green’s functions are computed using the frequency-wavenumber integration    *
# *    code (FKPROG) of Saikia (1994) based on the regional velocity model          *
# *     with a 10 Hz sampling rate.                                                 *
# *                                                                                 *
# *   Saikia, C.K., 1994. Modified frequency‐wavenumber algorithm for regional      *
# *      seismograms using Filon’s quadrature: modelling of Lg waves in eastern     *
# *      North America. Geophys. J. Int. 118, 142–158.                              *
# *      https://doi.org/10.1111/j.1365-246X.1994.tb04680.x                         *
# This is free software: you can redistribute it and/or modify                      *
# *                                                                                 *
# *    This program is distributed in the hope that it will be useful,              *
# *    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
# *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         *
# *                                                                                 *
# ***********************************************************************************/

import os
from obspy import Stream, Trace, read
import numpy as np
import datetime
import pathlib
from GF_functions import *


# [1] ----------------------- SET PATHS --------------------------------------
path2folder = pathlib.Path('../GFF').resolve()
Path2_bin = pathlib.Path('../BIN_HighSierra').resolve()
# Path2_bin = pathlib.Path('../BIN_Linux').resolve()
path2model = pathlib.Path('../models').resolve()
# [2] ----------------------- MODEL   ----------------------------------------
# MT_F = 'Gitt02'
MT_F = 'AK135'


# [3] ----------------------- RANGE   ----------------------------------------

DEPTH = np.arange(0, 100, 1)  # SET DEPTH RANGE
DIST0 = np.arange(1, 500, 5)  # SET DISTANCES FOR GREEN FUNCTION
# DEPTH = np.arange(5, 10, 5)  # SET DEPTH RANGE
# DIST0 = np.arange(30, 100, 10)  # SET DISTANCES FOR GREEN FUNCTION



# [4] ----------------------- PARAMETERS   ----------------------------------------
# do not change!
dH = 0.05 # Thickness of artificial layer in case Z == layer interface
f_low = 0.0
f_high = 0.0
npts0 = 4096
# npts0 = 1024
dt0 = 0.1
ReductionVel = 100


# [5] ----------------------- RUN   ----------------------------------------

# Load velocity model
MODEL = loadVmodel(MT_F,path2model)

# Generate distances
DIST2 = mk_DIST(DIST0)

# Make folde name
Path2_Green = set_path_gf(MT_F, path2folder, dt0, f_low, f_high, np.min(DIST2[DIST2 > 0]), np.max(DIST2))

# Make Geen's Functions
MODELN = np.zeros((MODEL.shape[0] + 1, 6))
mk_b2s(Path2_Green,npts0, dt0)
for kk in range(np.size(DIST2,0)):
    DIST = DIST2[kk]
    mk_perl1(Path2_Green, dt0, npts0, f_high, f_low)

    for ii in range(len(DEPTH)):
        mk_gf1(MODEL, MODELN, DEPTH[ii], npts0,dt0, Path2_bin, Path2_Green, MT_F, DIST, ReductionVel)
        GFlist = GetGreenFiles(Path2_Green, DEPTH[ii])
        GRN1 = Stream()
        mseed_dir = Path2_Green+'/MSEED/'
        if not os.path.exists(mseed_dir):
            os.mkdir(mseed_dir)
        for jj in range(len(GFlist)):
            [GRN0, nameG, depth, dist, pathFF] = Green2MSD(Path2_Green +'/'+ GFlist[jj])
            GRN0.write(mseed_dir + "%s.MSEED" % GFlist[jj], format="MSEED")













