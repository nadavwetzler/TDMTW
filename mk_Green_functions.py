
import os
import subprocess, sys
import platform
from obspy import Stream, Trace, read
import numpy as np
import datetime

now = datetime.datetime.now()
tt=now.timetuple()
b2s = "b2s.par"


#--------------------Set PARAMETERS --------------------
# DEPTH = np.array([5])
# DIST0 = np.array([255])
DEPTH = np.arange(0, 50, 1)  # SET DEPTH RANGE
DIST0 = np.arange(1, 100, 1)  # SET DISTANCES FOR GREEN FUNCTION
# DIST = np.arange(101, 201, 1)  # SET DISTANCES FOR GREEN FUNCTION
# DIST = np.arange(201, 301, 1)  # SET DISTANCES FOR GREEN FUNCTION
# DIST = np.arange(301, 401, 1)  # SET DISTANCES FOR GREEN FUNCTION


# DIST0 = np.arange(1,30, 1)  # SET DISTANCES FOR GREEN FUNCTION
# DEPTH = np.array([1])
# DIST = np.array([1])
# DIST = np.array([67, 100, 122, 142, 171, 201, 286, 311, 563])  # SET DISTANCES FOR GREEN FUNCTION
# DIST = np.array([200])  # SET DISTANCES FOR GREEN FUNCTION

dH = 0.05 # Thickness of artificial layer in case Z == layer interface

# f_low = 0.05
# f_high = 0.1
# npts0 = 512
# dt0 = 1.0

# f_low = 0.8
# f_high = 1.8
# npts0 = 2048
# dt0 = 0.1

# f_low = 0.1
# f_high = 1
# npts0 = 2048
# dt0 = 0.1

# ReductionVel = 100

# f_low = 0.008
# f_high = 0.03
# npts0 = 512
# dt0 = 1.0

# f_low = 0.2
# f_high = 0.5
# npts0 = 2048
# dt0 = 0.1

f_low = 0.0
f_high = 0.0
npts0 = 4096
dt0 = 0.1
ReductionVel = 100



# #-----------PREM--------------------------------
# MT_F = "PREMU_d"
# MODEL = np.zeros((23, 6))
# # #                Depth   Vp   Vs   Dens  Qa   Qb
# MODEL[0, 0:6] = [2.010, 4.78, 2.72, 2.10, 600, 300]  # 1.76
# MODEL[1, 0:6] = [9.760, 5.72, 3.20, 2.68, 600, 300]  # 1.79
# MODEL[2, 0:6] = [16.01, 6.27, 3.54, 2.77, 600, 300]  # 1.77
# MODEL[3, 0:6] = [28.01, 6.53, 3.69, 2.85, 600, 300]  # 1.77
# MODEL[4, 0:6] = [32.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
# MODEL[5, 0:6] =[40,8.1012,4.4849,3.3791,1446,600]
# MODEL[6, 0:6] =[60,8.0891,4.4771,3.3769,1447,600]
# MODEL[7, 0:6] =[80,8.0769,4.4695,3.3747,195,80]
# MODEL[8, 0:6] =[115,8.0554,4.4564,3.3709,195,80]
# MODEL[9, 0:6] =[150,8.0337,4.4436,3.3671,195,80]
# MODEL[10, 0:6] =[185,8.0118,4.4311,3.3633,195,80]
# MODEL[11, 0:6] =[220,7.9897,4.4188,3.3595,195,80]
# MODEL[12, 0:6] =[265,8.6455,4.6754,3.4626,365,143]
# MODEL[13, 0:6] =[310,8.7321,4.7069,3.4895,367,143]
# MODEL[14, 0:6] =[355,8.8187,4.7384,3.5164,370,143]
# MODEL[15, 0:6] =[400,8.9052,4.7699,3.5433,372,143]
# MODEL[16, 0:6] =[450,9.3899,5.0784,3.7868,365,143]
# MODEL[17, 0:6] =[500,9.6459,5.2243,3.8498,364,143]
# MODEL[18, 0:6] =[550,9.9018,5.3701,3.9128,363,143]
# MODEL[19, 0:6] =[600,10.158,5.5160,3.9758,362,143]
# MODEL[20, 0:6] =[635,10.212,5.5431,3.9840,362,143]
# MODEL[21, 0:6] =[670,10.266,5.5702,3.9921,362,143]
# MODEL[22, 0:6] =[721,10.910,6.0942,4.4124,744,312]

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


# #-------Gitterman 2005---------------------------------------
# MT_F = "Git05_d"
# MODEL = np.zeros((7, 6))
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [2.010, 4.78, 2.72, 2.10, 600, 300]  # 1.76
# MODEL[1, 0:6] = [9.760, 5.72, 3.20, 2.68, 600, 300]  # 1.79
# MODEL[2, 0:6] = [16.01, 6.27, 3.54, 2.77, 600, 300]  # 1.77
# MODEL[3, 0:6] = [28.01, 6.53, 3.69, 2.85, 600, 300]  # 1.77
# MODEL[4, 0:6] = [32.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
# MODEL[5, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[6, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

# -----------------GFZ Qsis ------------------------------------------
# MT_F = "GFZ_d"
# MODEL = np.zeros((9, 6))
#
# MODEL[0, 0:6] = [20.010, 5.800, 3.360, 2.7200, 600.0, 300.0]
# MODEL[1, 0:6] = [40.010, 5.800, 3.360, 2.7200, 600.0, 300.0]
# MODEL[2, 0:6] = [60.010, 6.500, 3.750, 2.9200, 600.0, 300.0]
# MODEL[3, 0:6] = [95.010, 6.500, 3.750, 2.9200, 600.0, 300.0]
# MODEL[4, 0:6] = [130.01, 8.040, 4.470, 3.3198, 600.0, 300.0]
# MODEL[5, 0:6] = [207.50, 8.045, 4.485, 3.3455, 600.0, 300.0]
# MODEL[6, 0:6] = [327.50, 8.050, 4.500, 3.3713, 600.0, 300.0]
# MODEL[7, 0:6] = [492.50, 8.175, 4.509, 3.3985, 100.00, 250.0]
# MODEL[8, 0:6] = [702.50, 8.300, 4.518, 3.4258, 100.00, 250.0]

# # #-------Kinneret ---------------------------------------
# * TOP - Zvi&Zvi Fig 1.13 + Michal’s model
# 0.0 7.0 13.0 20.0 30
# * VEL
# 2.8 4.7 5.5 6.7 8.0

# MT_F = "KinneretH_d"
# MODEL = np.zeros((8, 6))
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [2.100, 2.1, 1.20, 2.01, 100, 100]  # 1.76
# MODEL[1, 0:6] = [5.100, 2.5, 1.42, 2.10, 200, 100]  # 1.76
# MODEL[2, 0:6] = [7.100, 2.8, 1.60, 2.15, 300, 150]  # 1.76
# MODEL[3, 0:6] = [13.10, 4.72, 2.63, 2.68, 600, 300]  # 1.79
# MODEL[4, 0:6] = [20.10, 5.57, 3.43, 2.87, 600, 300]  # 1.77
# MODEL[5, 0:6] = [30.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
# MODEL[6, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[7, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

# # #-------Kinneret ---------------------------------------
# * TOP - Zvi&Zvi Fig 1.13 + Michal’s model
# 0.0 7.0 13.0 20.0 30
# * VEL
# 2.8 4.7 5.5 6.7 8.0

# MT_F = "GeoAzur_d"
# MODEL = np.zeros((5, 6))
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [0.600, 3.3, 1.90, 2.01, 200, 100]  # 1.76
# MODEL[1, 0:6] = [2.10, 4.5, 2.6, 2.30, 350, 175]  # 1.76
# MODEL[2, 0:6] = [5.100, 5.5, 3.18, 2.5, 500, 250]  # 1.76
# MODEL[3, 0:6] = [30.10, 6.5, 3.75, 2.90, 600, 300]  # 1.79
# MODEL[4, 0:6] = [1000, 8.1, 4.68, 3.30, 1000, 500]  # 1.77

# # #-------Kinneret ---------------------------------------


MT_F = "DSB_d"
MODEL = np.zeros((7, 6))
#               Depth   Vp   Vs   Dens  Qa   Qb       PR
MODEL[0, 0:6] = [2.500, 3.5, 1.94, 2.10, 600, 300]  # 1.8
MODEL[1, 0:6] = [7.500, 4.5, 2.5, 2.50, 600, 300]  # 1.8
MODEL[2, 0:6] = [14.10, 5.00, 2.90, 3.125, 600, 300]  # 1.6
MODEL[3, 0:6] = [20.10, 5.57, 3.43, 2.87, 600, 300]  # 1.77
MODEL[4, 0:6] = [30.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
MODEL[5, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
MODEL[6, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78
#
# # # #-------Kinneret ---------------------------------------
# # * TOP - Zvi&Zvi Fig 1.13 + Michal’s model
# # 0.0 7.0 13.0 20.0 30
# # * VEL
# # 2.8 4.7 5.5 6.7 8.0
#
# MT_F = "Kinneret6_d"
# MODEL = np.zeros((7, 6))
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [3.100, 3.9, 2.2, 2.10, 600, 300]  # 1.9
# MODEL[1, 0:6] = [7.100, 4.2, 2.5, 2.50, 600, 300]  # 1.8
# MODEL[2, 0:6] = [13.10, 4.72, 2.90, 2.68, 600, 300]  # 1.79
# MODEL[3, 0:6] = [20.10, 5.57, 3.43, 2.87, 600, 300]  # 1.77
# MODEL[4, 0:6] = [30.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
# MODEL[5, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[6, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

# #-------Gitterman 2005---------------------------------------
# MT_F = "GitHDUF05_d"
# MODEL = np.zeros((7, 6))
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [2.010, 4.78, 2.72, 2.10, 600, 300]  # 1.76
# MODEL[1, 0:6] = [9.760, 5.72, 3.20, 2.68, 600, 300]  # 1.79
# MODEL[2, 0:6] = [16.01, 6.27, 3.54, 2.77, 600, 300]  # 1.77
# MODEL[3, 0:6] = [28.01, 6.53, 3.69, 2.85, 600, 300]  # 1.77
# MODEL[4, 0:6] = [32.01, 6.97, 4.03, 2.88, 600, 300]  # 1.73
# MODEL[5, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[6, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

# #-------Costa Rica---------------------------------------
# MT_F = "costa_d"
# MODEL = np.zeros((9, 6))
# #                  Depth         Vp          Vs         Dens       Qa      Qb
# MODEL[0, 0:6] = [5.0100, 0.5350E+01, 0.3010E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[1, 0:6] = [13.010, 0.6120E+01, 0.3430E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[2, 0:6] = [16.010, 0.6280E+01, 0.3530E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[3, 0:6] = [20.010, 0.6460E+01, 0.3630E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[4, 0:6] = [25.010, 0.6720E+01, 0.3770E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[5, 0:6] = [30.010, 0.7010E+01, 0.3940E+01, 0.2786E+01, 600.00, 300.00]
# MODEL[6, 0:6] = [45.010, 0.7390E+01, 0.4150E+01, 0.3360E+01, 600.00, 300.00]
# MODEL[7, 0:6] = [50.010, 0.7550E+01, 0.4240E+01, 0.3360E+01, 600.00, 300.00]
# MODEL[8, 0:6] = [149.01, 0.8140E+01, 0.4570E+01, 0.3360E+01, 600.00, 300.00]





#-------INGV---------------------------------------
# MT_F = "INGV_d"
# MODEL = np.zeros((10, 6))
# #               Depth   Vp     Vs   Dens  Qa   Qb
# MODEL[0, 0:6] = [2.0, 3.50, 2.0, 2.00, 600, 300]
# MODEL[1, 0:6] = [8.0, 4.30, 2.43, 2.37, 600, 300]
# MODEL[2, 0:6] = [22.0, 6.05, 3.41, 2.70, 600, 300]
# MODEL[3, 0:6] = [39.0, 6.88, 3.94, 2.95, 600, 300]
# MODEL[4, 0:6] = [103.0, 8.0, 4.50, 3.31, 600, 300]
# MODEL[5, 0:6] = [253.0, 8.15, 4.66, 3.32, 600, 300]
# MODEL[6, 0:6] = [353.0, 8.3, 4.80, 3.35, 600, 300]
# MODEL[7, 0:6] = [453.0, 9.0, 5.14, 3.65, 600, 300]
# MODEL[8, 0:6] = [553.0, 9.6, 5.49, 3.84, 600, 300]
# MODEL[9, 0:6] = [1000.0, 10.3, 5.89, 4.07, 600, 300]

#-------GII JSTAR------------------------------------------
# MT_F = "GII02_d"
# MODEL = np.zeros((5, 6))
#
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [2.5900, 4.36, 2.41, 2.10, 600, 300]  #
# MODEL[1, 0:6] = [9.7900, 5.51, 3.10, 2.77, 600, 300]  #
# MODEL[2, 0:6] = [31.430, 6.23, 3.60, 2.85, 600, 300]  #
# MODEL[3, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[4, 0:6] = [1000.0, 8.15, 4.58, 3.20, 600, 300]  #

#-------GII JSTAR------------------------------------------
# MT_F = "GII03_d"
# MODEL = np.zeros((6, 6))
#
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [1.0100, 3.00, 1.70, 2.10, 600, 300]  #
# MODEL[1, 0:6] = [3.6000, 4.36, 2.41, 2.10, 600, 300]  #
# MODEL[2, 0:6] = [10.800, 5.51, 3.10, 2.77, 600, 300]  #
# MODEL[3, 0:6] = [32.440, 6.23, 3.60, 2.85, 600, 300]  #
# MODEL[4, 0:6] = [131.43, 7.95, 4.45, 2.88, 600, 300]  #
# MODEL[5, 0:6] = [1000.0, 8.15, 4.58, 3.20, 600, 300]  #

#-------DESERT---------------------------
# MT_F = "DESERT_d"
# MODEL = np.zeros((6, 6))
#
# #                Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [2.5, 3.7, 2.1, 2.10, 600, 300]  # 1.76
# MODEL[1, 0:6] = [7.5, 5.1, 2.7, 2.68, 600, 300]  # 1.79
# MODEL[2, 0:6] = [17.1, 6.3, 3.62, 2.77, 600, 300]  # 1.77
# MODEL[3, 0:6] = [30.1, 6.7, 3.75, 2.88, 600, 300]  # 1.73
# MODEL[4, 0:6] = [131.43, 7.9, 4.45, 2.88, 600, 300]  #
# MODEL[5, 0:6] = [1000.1, 8.15, 4.58, 3.20, 600, 300]  # 1.78

#---------- TESTER ------------------------------------
# MT_F = "TEST01_d"
# MODEL = np.zeros((4, 6))
#
# #               Depth   Vp   Vs   Dens  Qa   Qb       PR
# MODEL[0, 0:6] = [3.5, 4.30, 2.44, 2.80, 600, 300]  #
# MODEL[1, 0:6] = [5.0, 5.51, 3.31, 2.68, 600, 300]  #
# MODEL[2, 0:6] = [7.0, 6.23, 3.60, 3.00, 600, 300]  #
# MODEL[3, 0:6] = [50.44, 7.95, 4.60, 3.2, 600, 300]  #




#------------------Set path-------------------------
# Path2_Green = '%s/Dropbox/Moment-tensor/GREEN_%s_%s_%s-%s_fullNNN/' % (os.path.expanduser('~'), MT_F[0 : -1], dt0, f_low, f_high)

def mk_DIST(DIST0):
    maxD = 100
    if len(DIST0) < maxD:
        DIST = np.zeros((1,len(DIST0)))
        DIST[0][0:len(DIST0)] = DIST0
    else:
        ll = len(DIST0)
        nr = int(np.ceil(ll / maxD))
        nr1 = int(np.floor(ll / maxD))
        res = ll - nr1*maxD
        DIST = np.zeros((nr,maxD))
        for ii in range(nr1):
            DIST[ii][0:maxD] = DIST0[ii*maxD:((ii+1)*maxD)]
        if res > 0:
            DIST[nr1][0:res] = DIST0[(ii+1)*maxD:((ii+1)*maxD + res)]
    return DIST


def set_path_gf(name, dt0, f_low, f_high, DIST1, DIST2):
    Path2_Green = '%s/Dropbox/Moment-tensor/GREEN_%s_%s_%s-%s_fullNNN_%s_%s/' % (os.path.expanduser('~'), name[0 : -1], dt0, f_low, f_high, DIST1, DIST2)
    plat = platform.system()
    if plat == 'Darwin':
        v, _, _ = platform.mac_ver()
        v = float('.'.join(v.split('.')[:2]))
        if v == 10.1:
            Path2_MT = '%s/Dropbox/Moment-tensor/BIN_Yosemite/' % (os.path.expanduser('~'))  # MAC Yosemite
        elif v == 10.14:
            Path2_MT = '%s/Dropbox/Moment-tensor/BIN_HighSierra/' % (os.path.expanduser('~'))  # MAC Mojave
    else:
        Path2_MT = '%s/SeisSoft/Moment-tensor/TDMT/BIN/' % (os.path.expanduser('~'))  #  Linux GSI
    
    if os.path.exists(Path2_Green):
        os.chdir(Path2_Green)
    
    else:
        os.mkdir(Path2_Green)
        os.chdir(Path2_Green)
    return Path2_Green, Path2_MT



# ----------build b2s.par file --------------------------

def mk_b2s(Path2_Green):
    b2s_F_name = Path2_Green + b2s

    with open(b2s_F_name, 'w') as b2sf:
        b2sf.write('npts=%d\n' % npts0)
        b2sf.write('dt=%s\n' % dt0)
        b2sf.write('stime=0.0\n')
        b2sf.write('year=%d\n' % now.year)
        b2sf.write('jday=%d\n' % tt.tm_yday)
        b2sf.write('hour=%d\n' % now.hour)
        b2sf.write('min=%d\n' % now.minute)
        b2sf.write('sec=%f\n' % now.second)
        b2sf.write('msec=0\n')
        b2sf.write('ename="Synth"\n')

#-------------- make perl filter files -----------------------
def mk_perl1(Path2_Green, dt0, npts0, f_high, f_low):
    f_name = 'run_filtsyniso'
    if os.path.exists(f_name): os.remove(f_name)
    filt_P_name = Path2_Green + f_name
    filt_P = open(filt_P_name, 'w')
    filt_P.write('#! /bin/csh\n')
    filt_P.write('# Window out the eight vectors for many distances\n')
    filt_P.write('# Remember to set NT and DT CORRECTLY!\n')
    filt_P.write('#\n')
    filt_P.write('set   dt=%s\n' % dt0)
    filt_P.write('set npts=%d\n' % npts0)
    filt_P.write('set hcrn=%s\n' % f_high)
    filt_P.write('set lcrn=%s\n' % f_low)
    filt_P.write('#\n')
    filt_P.write('fromHelm < $1 > tmp2\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=0 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > tss\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=1 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > tds\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=2 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > xss\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=3 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > xds\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=4 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > xdd\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=5 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > zss\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=6 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > zds\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=7 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > zdd\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=8 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > rex\n')
    filt_P.write('window nt=$npts nx=10 nv=1 v0=9 < tmp2 > tmp3\n')
    filt_P.write('bin2sac par=b2s.par npts=$npts < tmp3 > zex\n')
    filt_P.write('sac << sacend\n')
    filt_P.write('setbb HCRN $hcrn\n')
    filt_P.write('getbb HCRN\n')
    filt_P.write('setbb LCRN $lcrn\n')
    filt_P.write('getbb LCRN\n')
    # filt_P.write('cut 0 300\n')
    filt_P.write('read tss tds xss xds xdd zss zds zdd rex zex\n')
    if f_high != 0:
        filt_P.write('bp co %LCRN %HCRN p 2\n')
    
    filt_P.write('write over\n')
    filt_P.write('quit\n')
    filt_P.write('sacend\n')
    filt_P.write('#\n')
    filt_P.write('#\n')
    filt_P.write('sac2bin in=tss out=tmp\n')
    filt_P.write('\mv tmp tss\n')
    filt_P.write('sac2bin in=tds out=tmp\n')
    filt_P.write('\mv tmp tds\n')
    filt_P.write('sac2bin in=xss out=tmp\n')
    filt_P.write('\mv tmp xss\n')
    filt_P.write('sac2bin in=xds out=tmp\n')
    filt_P.write('\mv tmp xds\n')
    filt_P.write('sac2bin in=xdd out=tmp\n')
    filt_P.write('\mv tmp xdd\n')
    filt_P.write('sac2bin in=zss out=tmp\n')
    filt_P.write('\mv tmp zss\n')
    filt_P.write('sac2bin in=zds out=tmp\n')
    filt_P.write('\mv tmp zds\n')
    filt_P.write('sac2bin in=zdd out=tmp\n')
    filt_P.write('\mv tmp zdd\n')
    filt_P.write('sac2bin in=rex out=tmp\n')
    filt_P.write('\mv tmp rex\n')
    filt_P.write('sac2bin in=zex out=tmp\n')
    filt_P.write('\mv tmp zex\n')
    filt_P.write('cat tss tds xss xds xdd zss zds zdd rex zex > tmp2\n')
    filt_P.write('mkHelm ntr=10 nt=$npts dt=$dt format="(6e12.5)" < tmp2 > $2\n')
    filt_P.write('rm tss tds xss xds xdd zss zds zdd rex zex tmp*\n')
    filt_P.close()
    subprocess.call(['chmod', '0777', filt_P_name])


# -----------------RUN OVER DEPTHS! ---------------------------------


def mk_gf1(MODEL, MODELN, DEPTH, npts0,dt0, Path2_MT):
    if os.path.exists('GREEN.1'): os.remove('GREEN.1')

    # Find the closest interface to selected depth
    kk = 0
    for jj in range(len(MODEL)):
        # Get the layer top depth
        if jj == 0:
            Z_top = -0.01
        else:
            Z_top = MODEL[jj-1,0]

        if DEPTH == MODEL[jj,0]:
            if Z_top == 0:
                H1 = MODEL[jj,0] - dH
                H2 = dH
            else:
                H1 = dH
                H2 = MODEL[jj + 1,0] - MODEL[jj,0] - dH
            MODELN[kk, 0:6] = MODEL[jj, 0:6]
            MODELN[kk, 0] = H1
            kk = kk + 1
            Layer_i = kk + 1
            MODELN[kk, 0:6] = MODEL[jj, 0:6]
            MODELN[kk, 0] = H2
            kk = kk + 1
        elif  MODEL[jj,0] > DEPTH and Z_top <= DEPTH:
            H1 = DEPTH - Z_top
            H2 = MODEL[jj,0] - DEPTH
            MODELN[kk, 0:6] = MODEL[jj, 0:6]
            MODELN[kk, 0] = H1
            kk = kk + 1
            Layer_i = kk + 1
            MODELN[kk, 0:6] = MODEL[jj, 0:6]
            MODELN[kk, 0] = H2
            kk = kk + 1

        else:
            MODELN[kk, 0:6] = MODEL[jj, 0:6]
            MODELN[kk, 0] = MODEL[jj,0] - Z_top
            kk = kk + 1



    # write model file
    MT_F_name = Path2_Green + MT_F + str(DEPTH) + ".model.txt"
    mtf = open(MT_F_name, 'w')
    mtf.write('.F.\n')
    mtf.write('    0   64\n')
    mtf.write('GREEN.1\n')
    #

        #              6.0     10.00       1  256  512    1.000     7    1
    mtf.write('    2.0%10.2f       1%5d%5d%9.3f    %s    1\n' % (DEPTH, npts0/2, npts0,dt0,MODELN.shape[0]))
    mtf.write('    1    1    1    1    1    1    1    1    1    1    0\n')


    for jj in range(MODELN.shape[0]):
        mtf.write(' %.4E %.4E %.4E %.4E   %5.2f    %5.2f\n' % (MODELN[jj, 0], MODELN[jj, 1], MODELN[jj, 2], MODELN[jj, 3],MODELN[jj, 4],MODELN[jj, 5]))

    # The layer below the depth!
    mtf.write('%5d\n' % Layer_i)
    mtf.write('  0.4000000E+03  1.500000E+00         0\n')



    # number of distances
    mtf.write('%5d  10000.0     30.0      1.0       0.8\n' % DIST.shape[0])

    # Distances!
    for jj in range(len(DIST)):
        mtf.write('%8.2f      0.0      %4.1f\n' % (DIST[jj], ReductionVel))

    mtf.close()

    # run FKRPROG!
    os.system(Path2_MT+"FKRPROG<%s" %MT_F_name)

    # run wvint9 d
    os.system(Path2_MT+"wvint9d")


    # ---------------Make perl fk file -------------------
    fk_name = 'run_fkrsortiso'
    if os.path.exists(fk_name): os.remove(fk_name)
    fk_P_name = Path2_Green + fk_name
    fk_P = open(fk_P_name, 'w')

    fk_P.write('#! /bin/csh\n')
    fk_P.write('# Window out the eight vectors for many distances\n')
    fk_P.write('# Remember to set NT and DT CORRECTLY!\n')
    fk_P.write('# USE wvint9 (The flipped traces are corrected within this\n')
    fk_P.write('# code rather than by the external program flip\n')
    fk_P.write('#\n')
    fk_P.write('set path=($path ../BIN)\n')
    fk_P.write('#\n')
    fk_P.write('set dt=%s\n' % dt0)
    fk_P.write('set npts=%d\n' % npts0)
    fk_P.write('##\n')
    fk_P.write('set dist=(')
    for jj in range(len(DIST)):
        fk_P.write('%s ' % DIST[jj])
    fk_P.write(')\n')
    fk_P.write('set loopend=%s\n' % len(DIST))
    fk_P.write('set depth=%s\n' % DEPTH)
    fk_P.write('##\n')
    fk_P.write('set count=0\n')
    fk_P.write('set j=1\n')
    fk_P.write('set vshift=0\n')
    fk_P.write('set i=0\n')
    fk_P.write('set nvec=0\n')
    fk_P.write('rehash\n')
    fk_P.write('#\n')
    fk_P.write('#\n')
    fk_P.write('@ nvec=($loopend - $count) * 10\n')
    fk_P.write('while ($count < $loopend)\n')
    fk_P.write('@ vshift=$i + 7\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp1$$\n')
    fk_P.write('@ vshift=$i + 4\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp2$$\n')
    fk_P.write('@ vshift=$i + 6\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp3$$\n')
    fk_P.write('@ vshift=$i + 3\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp4$$\n')
    fk_P.write('@ vshift=$i + 1\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp5$$\n')
    fk_P.write('@ vshift=$i + 5\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp6$$\n')
    fk_P.write('@ vshift=$i + 2\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp7$$\n')
    fk_P.write('@ vshift=$i\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp8$$\n')
    fk_P.write('@ vshift=$i + 8\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp10$$\n')
    fk_P.write('@ vshift=$i + 9\n')
    fk_P.write('window v0=$vshift e0=0 nt=$npts nx=$nvec nv=1 < vec > tmp9$$\n')
    fk_P.write('cat tmp1$$ tmp2$$ tmp3$$ tmp4$$ tmp5$$ tmp6$$ tmp7$$ tmp8$$ tmp9$$ tmp10$$ > junk\n')
    fk_P.write('echo $j\n')
    fk_P.write('mkHelm format="(6e13.5)" ntr=10 dt=$dt nt=$npts < junk > %s{$dist[$j]}d{$depth}.disp\n' % MT_F[0 : -1])
    fk_P.write('./run_filtsyniso %s{$dist[$j]}d{$depth}.disp %s{$dist[$j]}d{$depth}\n' % (MT_F[0 : -1], MT_F[0 : -1]))
    fk_P.write('rm tmp*$$\n')
    fk_P.write('@ i += 10\n')
    fk_P.write('@ count++\n')
    fk_P.write('@ j++\n')
    fk_P.write('end\n')
    fk_P.close()

    subprocess.call(['chmod', '0777', fk_P_name])


    os.system('./%s' % fk_name)

def Green2MSD(pathFF):
    name = pathFF.split('/')
    nameE = name[-1]
    nameG = name[-2]
    nameE = nameE.split('d')
    depth = float(nameE[1])
    dist = nameE[0].split('_')
    dist = float(dist[1])
    GRN0 = Stream()

    f_LOC=open(pathFF)
    lines=f_LOC.readlines()
    line3 = lines[3].split()
    smpls = line3[0]
    dt = line3[1]
    pos = []
    for jj in range(1, len(lines)):
        line = lines[jj].split()
        if line[0] == smpls:
            pos.append(jj)
    nl = pos[1] - pos[0]
    for jj in range(len(pos)):
        trace = []
        for ii in range(pos[jj]+1,pos[jj]+nl-1):
            line = lines[ii]
            dln = np.arange(0,len(line),12)
            for kk in range(len(dln)-1):
                num = float(line[dln[kk]:dln[kk+1]])
                trace.append(num)
        trace1 = Trace()
        trace1.data = np.asarray(trace)
        trace1.stats.delta = dt
        trace1.stats.station = str(dist)
        trace1.stats.location = str(depth)
        trace1.stats.channel = str(jj)
        SEIStmp = Stream(traces = trace1)
        GRN0 += SEIStmp
    return GRN0, nameG, depth, dist, pathFF


def GetGreenFiles(GreePath, depth):
    os.chdir(GreePath)
    disp_dir = GreePath+'d%s.DISP' % depth
    txt_dir = GreePath+'TXT'
    if not os.path.exists(disp_dir):
        os.mkdir(disp_dir)
    if not os.path.exists(txt_dir):
        os.mkdir(txt_dir)

    # Move all *.disp, *.txt files to different folder
    items = os.listdir(os.curdir)
    newlist = []
    for names in items:
        if names.endswith(".disp"):
            newlist.append(names)
            os.rename(GreePath + names, disp_dir+"/"+names)
        elif names.endswith(".txt"):
            os.rename(GreePath + names, txt_dir+"/"+names)
    # Get the name of the bigining of the GF
    os.chdir(disp_dir)
    items = os.listdir(os.curdir)
    name0 = items[0]
    name0 = name0.split('_')
    name0 = name0[0]

    # Get list of all GF
    os.chdir(GreePath)
    items = os.listdir(os.curdir)
    newlist = []
    for names in items:
        if names.startswith(name0):
            newlist.append(names)
    return newlist


############################################################################################################
DIST2 = mk_DIST(DIST0)
min_dist = np.min(DIST2[DIST2 > 0])
max_dist = np.max(DIST2)
Path2_Green, Path2_MT = set_path_gf(MT_F, dt0, f_low, f_high, min_dist, max_dist)
MODELN = np.zeros((MODEL.shape[0] + 1, 6))
mk_b2s(Path2_Green)
for kk in range(np.size(DIST2,0)):
    DIST = DIST2[kk]
    mk_perl1(Path2_Green, dt0, npts0, f_high, f_low)

    for ii in range(len(DEPTH)):
        mk_gf1(MODEL, MODELN, DEPTH[ii], npts0,dt0, Path2_MT)
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

















