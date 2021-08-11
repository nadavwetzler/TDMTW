# TDMTW GUI
# !/bin/env python
# /**********************************************************************************
# *    Copyright (C) by Nadav Wetzler                                               *
# *                                                                                 *
# *    Automatic and GUI Moment tensor inversion solver                             *
# *                                                                                 *
# *    TDMTW is free software: you can redistribute it and/or modify                *
# *                                                                                 *
# *    This program is distributed in the hope that it will be useful,              *
# *    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
# *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
# *                                                                                 *
# ***********************************************************************************/
import numpy as np
import matplotlib
import shapefile as shp
#import pandas as pd
import argparse
import socket
import subprocess
import sys
import copy
import os
import glob

from obspy.core.event import (EventDescription,NodalPlanes, NodalPlane,
                              Origin, Magnitude,
                              Tensor)
from scipy.spatial.transform import Rotation as ROT
from obspy.core.inventory import Inventory
from obspy.clients.fdsn.mass_downloader import RectangularDomain,Restrictions, MassDownloader
from obspy.core.event.base import Comment
from obspy.core.event import FocalMechanism as FocalMechanismEv
from obspy.core.event import Event as EventEv
from obspy.core.event import MomentTensor as MomentTensorEv
from obspy.signal.trigger import classic_sta_lta, trigger_onset
from tkinter import *
from PIL import ImageTk, Image
from obspy.clients.fdsn import Client, RoutingClient
from obspy import UTCDateTime, Stream, Trace, read
from obspy.geodetics.base import gps2dist_azimuth
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from obspy.taup import TauPyModel
from osm import OSM as OSM
from matplotlib import pyplot as plt
from obspy.imaging.beachball import beach, MomentTensor, mt2axes, mt2plane, aux_plane
from mpl_toolkits.basemap import Basemap
from obspy.signal.cross_correlation import xcorr
from obspy.geodetics import FlinnEngdahl


_fe = FlinnEngdahl()
FORMAT = '%(asctime)-15s %(clientip)s %(user)-8s %(message)s'
sys.path.insert(0, '/Users/nadavwetzler/Dropbox/Moment-tensor/Py_MT/TDMTW/')
model = TauPyModel(model="iasp91")
matplotlib.use("TkAgg")
font0 = FontProperties()
import warnings

warnings.filterwarnings("ignore")

# Commandline Arguments
parser = argparse.ArgumentParser(description='''This code performs the full moment tensor calculations based on Dregers TDMT.\n
                                 To run type: python3.6 tdmtw.py --Origintime "2019-05-15T16:53:47.4" --Lat0 32.8107 --Long0 32.7967 --Depth0 23 --AUTOMODE 0 --STN_LIST0 "ALLB" --Mw0 4.6''',
                                 epilog="By Nadav Wetzler 2019, nadav.wetzler@gmail.com")
parser.add_argument('-ot', '--Origintime', type=str, metavar='', help="Origint time: yyyy-mm-ddTHH:MM:SS.F")
parser.add_argument('-am', '--AUTOMODE', type=int, default=0, metavar='', help="0=Solve with GUI, 1=run auto mode")
parser.add_argument('-stn', '--STN_LIST0', type=str, default='ALLB', metavar='', help="Station type from following: ALLB, ALLBS, ALLA")
parser.add_argument('-dr0', '--DepthRang', type=int, default=0, metavar='', help="the range of depths solved by the inversion")
parser.add_argument('-lat', '--Lat0', type=float, metavar='', help="Location, latitude")
parser.add_argument('-lon', '--Long0', type=float, metavar='', help="Location, longitude")
parser.add_argument('-d', '--Depth0', type=float, metavar='', help="Location, depth")
parser.add_argument('-m', '--Mw0', type=float, metavar='', help="Magnitude")
parser.add_argument('-rmax', '--Radious', type=float, default=0, metavar='', help="Maximum distance from epicenter to station")
parser.add_argument('-rmin', '--Radious0', type=float, default=0, metavar='', help="Minimum distance from epicenter to station")
parser.add_argument('-c', '--Client', type=str, metavar='', help="Choose FDSN server like: ISN, IRIS, ORFEUS, GFZ, KOERI, https://www.fdsn.org/webservices/datacenters/")
parser.add_argument('-chnp', '--ChannelPriority', type=str, default='HBSE', metavar='', help="set the channel priority H, B, S, E (B=BH?)")
parser.add_argument('-invl', '--Invlength', type=int, default=-1, metavar='', help="length of the inversion (time)")
parser.add_argument('-gf', '--GreenFunction', type=str, default='0', metavar='', help="The name of the GF folder")
parser.add_argument('-fmax', '--Fmax', type=float, default=0, metavar='', help="Maximum frequency")
parser.add_argument('-fmin', '--Fmin', type=float, default=0, metavar='', help="Minimum frequency")
parser.add_argument('-taper', '--Taper', type=float, default=0, metavar='', help="Taper ratio right after download")
parser.add_argument('-taper2', '--Taper2', type=int, default=0, metavar='', help="Taper the inverted waveforms")

# parser.add_argument('-c', '--Client', type=str, default='ISN', metavar='', help="Choose FDSN server like: ISN, IRIS, ORFEUS, GFZ, KOERI, https://www.fdsn.org/webservices/datacenters/")
####################################################################################################################
# Running the code:
# python3.6 tdmtw_mseed.py -ot 2020-12-05T12:44:40 -lat 36.05 -lon 31.81 -d 98 -m 5.3 -fmin 0.02 -fmax 0.05 -c KOERI
ipSC3 = "10.0.3.15"  # '199.71.138.39'
#NTW = "IS,GE,CQ,BK,IU,KO"  #
NTW = "*"
dTime_search = 50
Dist_max = 80
trim_from_P = 1
dPt2 = 40
time_before_origin = 90  # seconds before t0
n_corners = 4
thresh_s2n = 10
dDays = 30  # number of days for seismicity plot

####################################################################################################################

def _get_resource_id(event_name, res_type, tag=None):
    """
    Helper function to create consistent resource ids.
    """
    res_id = "smi:local/tdmtw/%s/%s" % (event_name, res_type)
    if tag is not None:
        res_id += "#" + tag
    return res_id

def MK_MT_event(MT, ot, latitude, longitude, event_name, MWz, VR, DEPTH, Mo, GFtype, pdc, pclvd, piso, FMSp, STN_LISTN):
    """
    reference from: https://docs.obspy.org/_modules/obspy/io/cmtsolution/core.html

    :param MT:
    :param ot:
    :param latitude:
    :param longitude:
    :param event_name:
    :param MWz:
    :param VR:
    :param DEPTH:
    :param Mo:
    :param GFtype:
    :param pdc:
    :param pclvd:
    :param piso:
    :param FMSp:
    :param STN_LISTN:
    :return:
    """
    STN_LIST = ', '.join(STN_LISTN)
    pos_Mw = np.argmax(VR)
    depth = DEPTH[pos_Mw]
    m_w = MWz[pos_Mw]
    MT1 = MT[pos_Mw]
    FMS = FMSp[pos_Mw]
    FMS2 = aux_plane(FMS[0],FMS[1],FMS[2])
    tensor = Tensor(
        m_rr=MT1[0],
        m_pp=MT1[1],
        m_tt=MT1[2],
        m_rt=MT1[3],
        m_rp=MT1[4],
        m_tp=MT1[5]
    )

    # Create event with Moment tensor data.

    origin_time = UTCDateTime(ot)

    cat_origin = Origin(
        resource_id=_get_resource_id(event_name, "origin", tag="catalog"),
        time=origin_time,
        longitude=longitude,
        latitude=latitude,
        # Depth is in meters.
        depth=depth * 1000.0,
        origin_type="hypocenter",
        region=_fe.get_region(longitude=longitude, latitude=latitude),
        evaluation_status="preliminary"
    )


    tdmtw_mag = Magnitude(
        resource_id=_get_resource_id(event_name, "magnitude", tag="mw"),
        mag=round(m_w, 2),
        magnitude_type="mw",
        origin_id=cat_origin.resource_id
    )

    foc_mec = FocalMechanismEv(
        resource_id=_get_resource_id(event_name, "focal_mechanism"),
        # The preliminary origin most likely triggered the focal mechanism
        # determination.
        triggering_origin_id=cat_origin.resource_id,
        nodal_planes = NodalPlanes(nodal_plane_1=NodalPlane(strike=FMS[0], dip=FMS[1], rake=FMS[2]),
                                   nodal_plane_2=NodalPlane(strike=FMS2[0], dip=FMS2[1], rake=FMS2[2])) #  ,
                                   #  preferred_plane=[])
    )


    mt = MomentTensorEv(
        resource_id=_get_resource_id(event_name, "moment_tensor"),
        # Convert to Nm.
        scalar_moment=Mo / 1E7,
        tensor=tensor,
        double_couple = round(pdc[pos_Mw]) / 100,
        clvd = round(pclvd[pos_Mw]) / 100,
        iso = round(piso[pos_Mw]) / 100,
        variance_reduction=VR[pos_Mw],
        greens_function_id=GFtype,
        #comments = Comment(text = STN_LIST), #  resource_id=_get_resource_id(event_name, "moment_tensor", tag="stations")),
        method_id = "TDMTW"
    )

    foc_mec.moment_tensor = mt

    ev = EventEv(resource_id=_get_resource_id(event_name, "event"),
               event_type="earthquake")
    ev.event_descriptions.append(EventDescription(text=event_name,
                                                  type="earthquake name"))

    ev.origins.append(cat_origin)
    ev.magnitudes.append(tdmtw_mag)
    ev.focal_mechanisms.append(foc_mec)

    # Set the preferred items.
    ev.preferred_origin_id = cat_origin.resource_id.id
    ev.preferred_magnitude_id = tdmtw_mag.resource_id.id
    ev.preferred_focal_mechanism_id = foc_mec.resource_id.id

    return ev


def GetIP():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    ip0 = s.getsockname()[0]
    return ip0

def ReplacetoHome(ff):
    path0 = (os.path.expanduser('~'))
    ff0 = ff.split('"')
    ff=ff0[1]
    if ff[0] == '~':
        ff = path0 + ff[1:]
    return ff

def GetPathfile(pathfilename):

    path_2_iso = []
    Dir = []
    MT_G = []
    path2Localfaults = []
    PB_boundaries = []
    pathFF = []
    LOGO = []
    f=open(pathfilename)
    lines=f.readlines()

    for jj in range(1, len(lines)):
        line = lines[jj].split('=')
        if 'path_2_iso' == line[0]:
            path_2_iso = line[1]
            path_2_iso = ReplacetoHome(path_2_iso)
        elif 'MT_G' == line[0]:
            MT_G = line[1]
            MT_G = ReplacetoHome(MT_G)
        elif 'Dir' == line[0]:
            Dir = line[1]
            Dir = ReplacetoHome(Dir)
        elif 'path2Localfaults' == line[0]:
            path2Localfaults = line[1]
            path2Localfaults = ReplacetoHome(path2Localfaults)
        elif 'PB_boundaries' == line[0]:
            PB_boundaries = line[1]
            PB_boundaries = ReplacetoHome(PB_boundaries)
        elif 'pathFF' == line[0]:
            pathFF = line[1]
            pathFF = ReplacetoHome(pathFF)
        elif 'LOGO' == line[0]:
            LOGO = line[1]
            LOGO = ReplacetoHome(LOGO)
    return path_2_iso, Dir, MT_G, path2Localfaults, PB_boundaries, pathFF, LOGO

path_2_iso, Dir, MT_G, path2Localfaults, PB_boundaries, pathFF, LOGO =  GetPathfile('Pathsfile.txt')

def GetPstn(cat, STN):
    tp = 0
    cat1 = cat[0]
    for kk in range(len(cat1.picks)):
        if cat1.picks[kk].phase_hint == 'P':
            stn1 = cat1.picks[kk].waveform_id.station_code
            if stn1 == STN:
                tp = cat1.picks[kk].time
    return tp


def Mw2a(M):
    '''based on Eshelby model of circular fault and constant stress drop'''
    SD = 3 # [MPa]
    a = np.power(((7*np.power(10,(1.5*M+10.73)))/(16*SD*1000000)),(1/3))
    a = a/1000.0 # to km
    return a

def FMS2fType(Rake):
    if Rake >= -45 and Rake <= 45:
        iFault = 1 # Strike-slip
    elif Rake >= (180-45) and Rake <=180:
        iFault = 1 # Strike-slip
    elif Rake >= -180 and Rake < (-180+45):
        iFault = 1 # Strike-slip
    elif Rake > 90-45 and Rake <= 90+45:
        iFault = 3 # Reverse
    elif Rake >= -90-45 and Rake <= -90+45:
        iFault = 2 # Normal
    else:
        iFault = 1
    return iFault

def iFault2LW(ifault, m):
    if ifault == 1:  # strike-slip
        al = 4.33
        bl = 1.49
        aw = 3.80
        bw = 2.59
    elif ifault == 3:  # reverse
        al = 4.49
        bl = 1.49
        aw = 4.37
        bw = 1.95
    elif ifault == 2:  # normal
        al = 4.34
        bl = 1.54
        aw = 4.04
        bw = 2.11
    else:
        al = 4.38
        bl = 1.49
        aw = 4.06
        bw = 2.25


    l = np.power(10,((m - al)/bl))
    w = np.power(10,((m - aw)/bw))
    return l, w


def AUTOPICK(tr):
    # Get P use STA/LTA:
    df = tr.stats.sampling_rate
    cft1 = classic_sta_lta(tr.data, int(5 * df), int(10 * df))
    #cft2 = recursive_sta_lta(tr.data, int(5 * df), int(10. * df))
    on_of = trigger_onset(cft1, 1.5, 0.5)
    Pt = on_of[0, 0] / df
    PtT = tr.stats.starttime + Pt
    # plotSTALTA = 1
    # if plotSTALTA == 1: # Plotting the sta/lta results
    #     fig500 = plb.figure(500)
    #     ax = fig500.add_subplot(N_STN*2,1,ii+1)
    #     ax.plot(tr0.data, 'k')
    #     ax.set_ylabel(STNC[nn])
    #     ymin, ymax = ax.get_ylim()
    #     plb.vlines(on_of[:, 0], ymin, ymax, color='r', linewidth=2)
    #     plb.vlines(on_of[:, 1], ymin, ymax, color='b', linewidth=2)
    #     plb.subplot(N_STN*2,1,2*(ii+1), sharex=ax)
    #     plb.plot(cft2, 'k')
    #     plb.hlines([1.5, 0.5], 0, len(cft2), color=['r', 'b'], linestyle='--')
    #     plb.axis('tight')
    #     plb.show()
    return Pt, PtT



def GetClient(Lat0,Lon0,fdsn):
    CLIENT_ISN = "http://82.102.143.46:8181"


    fdsn_names = fdsn.split(',')
    fdsn_l = len(fdsn_names)
    if fdsn == 'AUTO':
        Latmin = 28
        Latmax = 36
        Lonmin = 31
        Lonmax = 37
        if Lat0 > Latmin and Lat0 < Latmax and Lon0 > Lonmin and Lon0 < Lonmax:

            out = 0
            faultsT = 1
            client = Client(CLIENT_ISN,user='wetzerna',password='Venus1')
            print('Using ISN')
            CLIENT = CLIENT_ISN
        elif Lon0 < 0:
            out = 1
            faultsT = 2
            # CLIENT = "http://service.iris.edu/"
            CLIENT = "IRIS"
            client = Client(CLIENT)
            print('Using IRIS')
        else:
            CLIENT = "http://geofon.gfz-potsdam.de"
            out = 1
            faultsT = 2
            client = Client(CLIENT)
            print('Using GEOFON')
    else:
        print('Using %d servers' % fdsn_l)
        servers = []

        for ii in range(fdsn_l):
            fdsn_name = fdsn_names[ii]
            if fdsn_name == 'ISN':
                out = 0
                faultsT = 1
                client = Client(CLIENT_ISN,user='wetzerna',password='Venus1')
                CLIENT = CLIENT_ISN
            elif fdsn_name == 'GFZ':
                CLIENT = "http://geofon.gfz-potsdam.de"
                out = 1
                faultsT = 2
                client = Client(CLIENT)
            elif fdsn_name == 'RC':
                client = RoutingClient("iris-federator")
                #client = RoutingClient("eida-routing")
                out = 1
                faultsT = 2
            elif fdsn_name == 'GNSr':
                CLIENT = "http://service-nrt.geonet.org.nz"
                # CLIENT = 'http://service.geonet.org.nz'
                out = 1
                faultsT = 2
                client = Client(CLIENT)
            elif fdsn_name == 'NOA':
                CLIENT = "http://eida.gein.noa.gr/"
                out = 1
                faultsT = 2
                client = Client(CLIENT)
            elif fdsn_name == 'OVSICORI':
                CLIENT = "http://10.10.128.91:8080"
                out = 1
                faultsT = 2
                client = Client(CLIENT)
            elif fdsn_name == 'GSI':
                # CLIENT = 'http://172.16.46.102:8181'
                CLIENT = 'http://172.16.46.102:8181'
                client = Client(CLIENT)
                out = 1
                faultsT = 2
            elif fdsn_name == 'SYNTH':
                CLIENT = 'http://172.16.46.140:8181'
                client = Client(CLIENT,user='test',password='test')
                out = 1
                faultsT = 2
            else:
                CLIENT = fdsn_name
                out = 1
                faultsT = 2
                client = Client(CLIENT)
            servers.append(client)
            print('Using FDSN: %s %s' % (client.base_url, CLIENT))
        client = servers

    # else:
    #     print('Using %d servers' % fdsn_l)
    #     servers = []
    #     for ii in range(fdsn_l):
    #         if fdsn_names[ii] == 'ISN':
    #             CLIENT = CLIENT_ISN  # ISRAEL fdsn
    #             out = 0
    #             faultsT = 1
    #             client = Client(CLIENT,user='wetzerna',password='Venus1')
    #             servers.append(client)
    #         else:
    #             client = Client(fdsn_names[ii])
    #             servers.append(client)
    #             out = 1
    #             faultsT = 1
    #     client = servers



    return client, out, faultsT, fdsn_l

def INV_LENGTH_T(Mw0, useVal):
    if useVal > 0:
        Sampl4inv = useVal
    else:
        Sampl4inv    = 25 # in seconds
        if Mw0 >=3.2 and Mw0 < 3.5:
            Sampl4inv    = 50
        elif Mw0 >= 3.5 and Mw0 < 4:
            Sampl4inv    = 120
        elif Mw0 >= 4 and Mw0 < 5:
            Sampl4inv = 150
        elif Mw0 >= 5 and Mw0 < 6:
            Sampl4inv = 350
        elif Mw0 >= 6:
            Sampl4inv = 350
        # For now 350 sec. is the maximum!
        # elif Mw0 > 7:
        #     Sampl4inv = 380
    dPt1 = np.round(Sampl4inv / 3)
    return Sampl4inv, dPt1

def DISP_ROT_FILT(SEIS1, az, cutlow, fcuthigh, n_corners, dt0, dPt1):
    SEIS = SEIS1.copy()

    location = SEIS[0].stats.location
    if cutlow > 0.002:
        pre_filt0 = [0.001, 0.005, 45, 50]
    TimeT = SEIS[0].stats.endtime-SEIS[0].stats.starttime
    if location == '21':
        wl = 120
    else:
        wl = 60
    SEIS.detrend("linear")
    #SEIS.detrend("demean")
    SEIS.taper(max_percentage=dPt1/TimeT, type='cosine')
    try:
        SEIS.remove_response(output="DISP", taper=True, water_level= wl, pre_filt=pre_filt0)
    except:

        print('Removing response manually!!!')
        if location == '21':
            SEIS.integrate()
            SEIS.integrate()
        elif location == '22':
            SEIS.integrate()
        for ii in range(len(SEIS)):
            SEIS[ii].data = SEIS[ii].data * 10**-7


    for tr in SEIS: # from M to cm:
        tr.data = tr.data * 100

    SEIS.rotate(method='NE->RT', back_azimuth=az)
    SEIS.detrend("linear")

    #SEIS.taper(max_percentage=0.1, type='cosine')
    SEIS.taper(max_percentage=dPt1/TimeT, type='cosine')
    SEISRD = SEIS.copy()
    SEIS.filter('bandpass', freqmin=cutlow, freqmax=fcuthigh, corners=n_corners, zerophase=True)
    SEIS.interpolate(sampling_rate=int(1/dt0), npts=int(TimeT / dt0))
    return SEISRD, SEIS

def MSEED2GF(GRN0, depth, dist0, name0):
    '''

    '''

    # Find closest distance from availble GF
    Distv = []
    Depthv = []
    for ii in range(len(GRN0)):
        Distv.append(float(GRN0[ii].stats.station))
        Depthv.append(float(GRN0[ii].stats.location))
    Distv = set(Distv)
    Depthv = set(Depthv)
    DIST = np.asarray(Distv)
    pos = np.argmin(np.abs(DIST - dist0))
    dist = GRN0.DIST[pos]
    if np.abs(dist-dist0) > 5:
        print("GF is different than stn dist by 5 km!!")

    name = "%s%sd%sN" %(name0, dist, depth)
    SEIS = GRN0.select(station=dist, location=depth, channel="?")
    S2H = open(name, 'w')
    format = '6e12.5'
    S2H.write('      10\n')
    S2H.write('(%s)\n' % format)
    for gg in range(10):
        tr = SEIS.select(component = (gg +1))
        tr = tr[0]
        ltr = len(tr)
        S2H.write('     0.0000e+00     0.0000e+00      0  0  0.00\n')
        S2H.write('%d %f   0.0000e+00\n' % (ltr, tr.stats.delta))
        ntr = int(np.floor(ltr / 6))
        ntrl = ltr - ntr*6
        for jj in range(ntr):
            S2H.write('%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n'
                      % (tr[jj*6], tr[jj*6+1], tr[jj*6+2], tr[jj*6+3], tr[jj*6+4], tr[jj*6+5]))
        if ntrl > 0:
            for jj in range(ntrl):
                S2H.write('%12.5e ' % tr[6*ntr + jj])
            S2H.write('\n')
    S2H.close()

def decimate2lowsmpl(SEIS):
    ll = len(SEIS)
    smpls = np.zeros(ll)
    for ii in range(ll):
        smpls[ii] = SEIS[ii].stats.sampling_rate

    for ii in range(ll):
        if SEIS[ii].stats.sampling_rate > min(smpls):
            print('Decimating %s from %d to %d' % (SEIS[0].stats.station, max(smpls), min(smpls)))
            SEIS[ii].decimate(int(SEIS[ii].stats.sampling_rate / min(smpls)))
    return SEIS

def GetGreenFiles(GreePath):
    os.chdir(GreePath)
    disp_dir = GreePath+'DISP'
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
    dt = float(line3[1])
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


    smpls = int(smpls)
    return GRN0, nameG, depth, dist, pathFF, smpls, dt

def MakeMSEED(GreePath):
    GFlist = GetGreenFiles(GreePath)
    GRN1 = Stream()

    for ii in range(len(GFlist)):
        print(GFlist[ii])
        [GRN0, nameG, depth, dist, pathFF] = Green2MSD(GreePath + GFlist[ii])
        GRN0.write(pathFF + ".MSEED", format="MSEED")
        #GRN1 += GRN0

    return GRN1, nameG

def Green2gg(GRN, vol_p):
    N_STN = len(GRN)
    # Load all Green Functions
    GRN0 = Stream()
    dist0 = []
    for ii in range(N_STN):
        [GRN1, nameG, depth, dist, pathFF, smpls, dt] = Green2MSD(GRN[ii])
        GRN0 += GRN1
        dist0.append(dist)
    dist0 = np.asarray(dist0)

    # Make gg (green function matrix)

    gg = np.zeros((N_STN, 10, smpls))
    for ii in range(N_STN):
        GRN1 = GRN0.select(station=str(dist0[ii]))
        gg[ii][0][:] = GRN1[0].data
        gg[ii][1][:] = GRN1[1].data
        gg[ii][2][:] = GRN1[2].data
        gg[ii][3][:] = GRN1[3].data
        gg[ii][4][:] = GRN1[4].data
        gg[ii][5][:] = GRN1[5].data * -1.0
        gg[ii][6][:] = GRN1[6].data * -1.0
        gg[ii][7][:] = GRN1[7].data * -1.0
        if vol_p == 1:
            gg[ii][8][:] = GRN1[8].data * -1.0
            gg[ii][9][:] = GRN1[9].data * -1.0
    return gg

def Green2ggMSEED(GRN, vol_p, cutlow, fcuthigh, dt0):
    N_STN = len(GRN)
    # Load all Green Functions

    GRN0 = Stream()
    dist0 = []
    for ii in range(N_STN):
        # [GRN1, nameG, depth, dist, pathFF, smpls, dt] = Green2MSD(GRN[ii])
        GRN1 = read(GRN[ii])
        pathFF = GRN[ii]
        name = pathFF.split('/')
        nameE = name[-1]
        nameE = nameE.split('d')
        depth = float(nameE[1].split('.')[0])
        dist = float(nameE[0].split('_')[1])
        GRN0 += GRN1
        dist0.append(dist)
    dist0 = np.asarray(dist0)

    # Make gg (green function matrix)
    n_corners = 4
    GRN0.filter('bandpass', freqmin=cutlow, freqmax=fcuthigh, corners=n_corners, zerophase=True)
    TimeT = GRN0[0].stats.endtime - GRN0[0].stats.starttime
    GRN0.interpolate(sampling_rate=int(1/dt0), npts=int(TimeT / dt0))
    smpls = len(GRN0[0].data)

    gg = np.zeros((N_STN, 10, smpls))
    for ii in range(N_STN):
        GRN1 = GRN0.select(station=str(dist0[ii]))
        gg[ii][0][:] = GRN1[0].data
        gg[ii][1][:] = GRN1[1].data
        gg[ii][2][:] = GRN1[2].data
        gg[ii][3][:] = GRN1[3].data
        gg[ii][4][:] = GRN1[4].data
        gg[ii][5][:] = GRN1[5].data * -1.0
        gg[ii][6][:] = GRN1[6].data * -1.0
        gg[ii][7][:] = GRN1[7].data * -1.0
        if vol_p == 1:
            gg[ii][8][:] = GRN1[8].data * -1.0
            gg[ii][9][:] = GRN1[9].data * -1.0
    return gg

def whichGreen(Mw0, out, GF0, Fmin0, Fmax0):
    min_depth = 1
    if GF0 == '0':
        Green = 'GREEN_GitHDUF05__0.1_0.0-0.0_fullNNN_5.0_495.0/MSEED/'
        max_depth = 99
        dt = 0.1

    else:
        Green = GF0
        max_depth = 49
        dt = 0.1

    if Mw0 < 3.2:
        f1 = 0.2
        f2 = 0.6
    elif Mw0 >= 3.2 and Mw0 < 3.5:
        f1 = 0.2
        f2 = 0.6
    elif Mw0 >= 3.5 and Mw0 < 3.9:
        f1 = 0.05
        f2 = 0.1
    elif Mw0 >= 3.9 and Mw0 < 5.0:
        f1 = 0.05
        f2 = 0.1
    elif Mw0 >= 5.0 and Mw0 < 6.5:
        f1 = 0.02
        f2 = 0.05
        dt = 1.0
    elif Mw0 >= 6.5:
        f1 = 0.01
        f2 = 0.035
        dt = 1.0

    if dt < float(Green.split('_')[3]):
        dt  = float(Green.split('_')[3])
        print('sample rate of GF is now: %f' % float(Green.split('_')[3]))


    print(Green)
    if Fmin0 > 0:
        f1 = Fmin0
    if Fmax0 > 0:
        f2 = Fmax0
        if 1 / Fmax0 > 10:
            dt = 1.0
        else:
            dt = 0.1

    return Green, min_depth, max_depth, f1, f2, dt


def MW2DIST(Mw):
    Dist_min = 0
    Dist_max_0 = 20
    # Dist_max = Mw * 100.0 # km
    Dist_max = 0.024 * np.power(Mw, 6.7)
    if Dist_max> 500:
        Dist_max = 500
    elif Dist_max < Dist_max_0:
        Dist_max = Dist_max_0
    if Mw > 5:
        Dist_min = 10
    return Dist_max, Dist_min

def SIGNAL2NOISE(tp, SEIS, timew, plotter):
    ''' return the signal 2 noise ratio befor and after the P'''
    # timew = 30 # in sec.
    S2N = np.zeros(len(SEIS))
    for ii in range(len(SEIS)):
        tr = SEIS[ii]
        smpl = tr.stats.sampling_rate

        # Length of time window for Signal
        n1 = int(timew * smpl)

        # position of P
        ntp = int(tp * smpl)

        # Make noise
        noise = tr.data[0:ntp]

        # Make Signal
        length1 = len(tr.data)
        pos_end = int(np.min([length1, n1]))
        #signal = tr.data[ntp:ntp+pos_end]
        signal = tr.data[ntp:-1]

        # Calc S2N
        # Anoise  = np.mean(np.abs(noise))
        # Asignal = np.mean(np.abs(signal))
        Anoise  = np.max(np.abs(noise))
        Asignal = np.max(np.abs(signal))
        S2N[ii] = Asignal / Anoise
        S2N[ii] = np.round(S2N[ii])
        if plotter == 1:
            f = Figure(figsize=(2,5), dpi=100)
            ax = f.add_subplot(111)
            ax.plot(tr.data)
            ax.plot([n1,n1],[min(tr.data),max(tr.data)])
            ax.plot([ntp,ntp],[min(tr.data),max(tr.data)])
            ax.plot([ntp+n1,ntp+n1],[min(tr.data),max(tr.data)])
            plt.show()
    S2N = np.max(S2N)
    return S2N

def HH1_HHN(SEIS):
    for ii in range(len(SEIS)):
        if SEIS[ii].stats.channel == 'HH1':
            SEIS[ii].stats.channel = 'HHN'
        elif SEIS[ii].stats.channel == 'HH2':
            SEIS[ii].stats.channel = 'HHE'



def rot_to_north(SEIS, rot_ang):
    if rot_ang != 0:
        st = SEIS.copy()
        rot_ang2 = make360(-(rot_ang+180))
        st.rotate(method="NE->RT", back_azimuth=rot_ang2)
        n = st.select(component="R")[0]
        e = st.select(component="T")[0]

        n.stats.channel = n.stats.channel[:2] + "N"
        e.stats.channel = e.stats.channel[:2] + "E"

        for tr in SEIS.select(channel='??E'):
            SEIS.remove(tr)
        for tr in SEIS.select(channel='??N'):
            SEIS.remove(tr)

        SEIS+=n
        SEIS+=e

    return SEIS


def LOAD_STN_SEIS(stn_name, chn, Pt, Sampl4inv, client, dPt1, Taper0, rot_angl,INT):
    STNC = 'NEZ'
    STNC2 = '12Z'
    isData = 1
    rotated = 1
    SEIS = Stream()
    for nn in range(len(STNC)):
        try:
            if stn_name == 'EIL':
                chn = 'HH'
                SEIStmp = client.get_waveforms(NTW, stn_name, "*",
                                               "%s%s" % (chn, STNC[nn]),
                                               Pt - time_before_origin,
                                               Pt + Sampl4inv + dPt2 + 20,
                                               attach_response=True)

            # elif stn_name == 'SLTI':
            #     isData = 0

            else:
                try:
                    SEIStmp = client.get_waveforms(NTW, stn_name, "*",
                                                   "%s%s" % (chn, STNC[nn]),
                                                   Pt - time_before_origin,
                                                   Pt + Sampl4inv + dPt2 + 20,
                                                   attach_response=True)
                except:
                    try:
                        SEIStmp = client.get_waveforms(NTW, stn_name, "*",
                                                   "%s%s" % (chn, STNC2[nn]),
                                                   Pt - time_before_origin,
                                                   Pt + Sampl4inv + dPt2 + 20,
                                                   attach_response=True)
                        rotated = 0
                    except:
                        try:
                            chn = 'BH'
                            SEIStmp = client.get_waveforms(NTW, stn_name, "*",
                                                           "%s%s" % (chn, STNC[nn]),
                                                           Pt - time_before_origin,
                                                           Pt + Sampl4inv + dPt2 + 20,
                                                           attach_response=True)
                        except:
                            chn = 'SH'
                            SEIStmp = client.get_waveforms(NTW, stn_name, "*",
                                                           "%s%s" % (chn, STNC[nn]),
                                                           Pt - time_before_origin,
                                                           Pt + Sampl4inv + dPt2 + 20,
                                                           attach_response=True)

            # SEIStmp.plot()
            


            if Taper0 > 0:
                SEIStmp.detrend()
                SEIStmp.taper(max_percentage=0.2, type='cosine')
            SEIStmp.merge()
            #
            # SEIStmp.plot()
            SEIStmp = SEIStmp[0]
            SEIS += SEIStmp

            # SEIS.plot()
        except:
            isData = 0
            # print('No data for %s %s' % (stn_name, STNC[nn]))
    # if len(SEIS) > 3:
    #     isData = 0
    #     print('More than 3 chn for %s' % stn_name)
    if len(SEIS) == 3:
        SEIS = decimate2lowsmpl(SEIS)
        initc = Pt - dPt1  #2018-10-25T22:53:59.991054Z
        enddc = initc + Sampl4inv + dPt2 # 2018-10-25T23:00:29.991054Z
        SEIS = SEIS.trim(starttime=initc, endtime=enddc, nearest_sample=True, pad=True, fill_value=0)
        SEIS.attach_response(INT)
        ## Rotate HH1/HH2 --> HHE/HHN
        channels1 = [tr.stats.channel for tr in SEIS]
        check_rot = False
        for chni in channels1:
            if chni[2] == '1':
                check_rot = True
                chan1 = chni
        if check_rot == True:
            print("%s Rotating %s1/%s2 to %sE/%sN" % (stn_name, chan1[0:2], chan1[0:2],chan1[0:2],chan1[0:2]))
            #az1 = inv.get_channel_metadata('%s.%s..HH1' % (network, station))
            #az2 = inv.get_channel_metadata('%s.%s..HH2' % (network, station))
            SEIS.rotate('->ZNE', inventory=INT)


        if np.abs(rot_angl) > 0:
            SEIS = rot_to_north(SEIS, rot_angl)
            print('rot %s with %d deg.' % (stn_name, rot_angl))

        # SEIS.plot()
        SEIS.detrend()
        if stn_name == 'ZNT':
            print(' flip channels for ZNT stn')
            for ii in range(len(SEIS)):
                SEIS[ii].data = -SEIS[ii].data
        if stn_name == 'DOR':
            print(' flip channels for DOR stn')
            for ii in range(len(SEIS)):
                SEIS[ii].data = -SEIS[ii].data
        if stn_name == 'ATR':
            print(' flip channels for ATR stn')
            for ii in range(len(SEIS)):
                SEIS[ii].data = -SEIS[ii].data
        if stn_name == 'DOR':
            print(' flip channels for DOR stn')
            for ii in range(len(SEIS)):
                SEIS[ii].data = -SEIS[ii].data
        if stn_name == 'ATZ':
            print(' flip channels for ATZ stn')
            for ii in range(len(SEIS)):
                SEIS[ii].data = -SEIS[ii].data

    else:
        isData = 0
    return SEIS, isData, chn

def MSEED2HELM(path, stn_name, SEIST):
    '''
    ---- Py version for sac2helm ----------------------
    os.system(path_2_iso+"sac2helm out=%s" %STN_LIST[ii])
    :param path: the path of the exe file of the TDMT
    :param stn_name: station name
    :param SEIST: Streams - three channels seismic data
    :return:
    '''
    S2H = open(path +'/'+ stn_name, 'w')
    numtype = '(7e14.5)'
    order = 'TRZ'
    S2H.write('       3\n')
    S2H.write('%s\n' % numtype)
    for gg in range(len(SEIST)):
        tr = SEIST.select(component = order[gg]).copy()
        tr = tr[0]
        ltr = len(tr)
        S2H.write('     0.0000e+00     0.0000e+00      0  0  0.00\n')
        S2H.write('%d %f   0.0000e+00\n' % (ltr, tr.stats.delta))
        ntr = int(np.floor(ltr / 7))
        ntrl = ltr - ntr*7
        for jj in range(ntr):
            S2H.write('%14.5e %14.5e %14.5e %14.5e %14.5e %14.5e %14.5e\n'
                      % (tr[jj*7], tr[jj*7+1], tr[jj*7+2], tr[jj*7+3], tr[jj*7+4], tr[jj*7+5], tr[jj*7+6]))
        if ntrl > 0:
            for jj in range(ntrl):
                S2H.write('%14.5e ' % tr[7*ntr + jj])
            S2H.write('\n')
    S2H.close()


def MKPLOTDIR(pathplt0, STN_LISTN):
    #------Creat the folder by STN names----------
    Folder_name_stn = '_'.join(sorted(STN_LISTN))
    if os.path.exists(pathplt0 + '/' + Folder_name_stn):
        pathplt = pathplt0 + '/' + Folder_name_stn
    else:
        pathplt = pathplt0 + '/' + Folder_name_stn
        os.mkdir(pathplt)
    return(pathplt)


def MKFOLDERS(Event, Dir):
    Event = Event[0]
    EventF = Event.replace(':', '-')
    EventN = Event.replace(':', '.')
    MT_F_name = "mt_inv.in"
    # Make SAC folder
    MT_F_name = Dir + EventF + "/SAC/" + MT_F_name
    if os.path.exists(Dir + EventF + "/SAC/"):
        pathdir = Dir+EventF + "/SAC"
    else:
        pathdir = Dir + EventF
        os.mkdir(pathdir)
        pathdir = pathdir + "/SAC"
        os.mkdir(pathdir)

    # Make PLOTs folder
    if os.path.exists(Dir + EventF + "/PLT/"):
        pathplt = Dir+EventF + "/PLT"
    else:
        pathplt = Dir + EventF + "/PLT"
        os.mkdir(pathplt)

    # Make Output folder
    if os.path.exists(Dir + EventF + "/OUT/"):
        pathOut = Dir+EventF + "/OUT"
    else:
        pathOut = Dir + EventF + "/OUT"
        os.mkdir(pathOut)

    # Make green-functions folder
    if os.path.exists(Dir + EventF + "/GF/"):
        pathGF = Dir+EventF + "/GF"
    else:
        pathGF = Dir + EventF + "/GF"
        os.mkdir(pathGF)


    return pathOut, pathplt, pathdir, MT_F_name, EventF, EventN, pathGF

def SetDepthRang(Depth0, Min_Depth, Max_Depth, dr0):
    if dr0 == 0:
        dr1 = 10
        dr2 = 15
    else:
        dr1 = dr0
        dr2 = dr0


    D1 = Depth0 - dr1
    if D1 < Min_Depth:
        D1 = Min_Depth
    D2 = Depth0 + dr2

    if D2 >= Max_Depth:
        D2 = Max_Depth

    if Depth0 > Max_Depth:
        D1 = Max_Depth - 3
        D2 = Max_Depth
    return D1, D2

def FindIDinFile(pathFF, eventID):
    f_LOC=open(pathFF)
    lines=f_LOC.readlines()
    pos = -1
    for jj in range(1, len(lines)):
        line = lines[jj].split(',')
        if eventID == line[0]:
            pos = jj
    return pos

def ADDeventF(pathFF, EVENT):
    names0 = ['eventID', 'T0', 'lat', 'lon', 'depth', 'strike', 'dip', 'rake', 'Mw', 'VR', 'MODE']
    if not os.path.exists(pathFF):
        FF = pd.DataFrame(columns = names0)
        FF.to_csv(path_or_buf=pathFF, encoding='utf-8', index=False)
    FF = pd.read_csv(pathFF, delimiter = ',')
    pos = FindIDinFile(pathFF, EVENT[0])
    FF2 = pd.DataFrame([EVENT], columns = names0)
    if pos == -1:
        FF = FF.append(FF2)
    else:
        FF.iloc[pos-1] = FF2.iloc[0]
    os.remove(pathFF)
    FF.to_csv(path_or_buf=pathFF, encoding='utf-8', index=False)



def Get_Reg_EQ(client, dDays, Origine_time):

    dsec = dDays*24*3600
    okcat = 1
    try:
        cat = client.get_events(starttime=Origine_time-dsec, endtime=Origine_time+dsec, minmagnitude=1)
    except:
        okcat = 0
        cat = 0
    return okcat, cat

def Calc_eq_reg_R(cat, Lat0, Lon0, Origine_time):
    """
    The function returns the distances and time gap between a set of  earthquakes and a main event
    :param cat: Catalog Class from Obspy e.g. cat = client.get_catalog(???)
    :param Lat0: latitude of the main event
    :param Lon0: longitude of the main event
    :param Origine_time: origin time of the main event (via utctime)
    :return: vetors contaning the distances and the catalog
    """
    ll = len(cat)
    R = []
    Lat = []
    Lon = []
    Mag = []
    Tt0 = []
    for ii in range(ll):
        # print(ii)
        try:
            Lat.append(cat.events[ii].preferred_origin().latitude)
            Lon.append(cat.events[ii].preferred_origin().longitude)
            Mag.append(cat.events[ii].preferred_magnitude().mag)
            [Ri, az, baz] = gps2dist_azimuth(Lat[ii], Lon[ii], Lat0, Lon0, a=6378137.0, f=0.0033528106647474805)
            R.append(Ri)
            Tt0.append(cat.events[ii].preferred_origin().time - Origine_time)
        except:
            pass

    R = np.array(R)
    Lat = np.array(Lat)
    Lon = np.array(Lon)
    Mag = np.array(Mag)
    Tt0 = np.array(Tt0)
    return R, Lat, Lon, Mag, Tt0

def Make_circle_R(R, Lat0, Long0):
    """
    The function returns two vectors with coordinated: Lon, Lat of a circle based on an input location and radius
    :param R: Radius length in KM
    :param Lat0: center of circle: Latitude
    :param Long0: center of circle: Longitude (positive E)
    :return: x, y coordinates of a circle
    """
    R = R * 0.009
    degr = np.arange(0, 360, 1)
    degr = np.radians(degr)
    x_c = Long0 + R *np.cos(degr)
    y_c = Lat0 + R*np.sin(degr)
    return x_c, y_c

def write2SAC(SEIS0,pathdir):
    # Write original to SAC file
    for tr in SEIS0:
        tr.write(pathdir+"/"+tr.id + ".O.SAC", format="SAC")

def Seis2ss(SEIS_DISP_ROT_F, STN_LISTN):
    N_STN = len(STN_LISTN)
    order = 'TRZ'
    ltr = len(SEIS_DISP_ROT_F[0].data)

    ss = np.zeros((N_STN, 3, ltr))
    for ii in range(N_STN):
        SEIST = SEIS_DISP_ROT_F.select(station=STN_LISTN[ii]).copy()
        for kk in range(len(order)):
            tr = SEIST.select(component = order[kk]).copy()
            ss[ii][kk][:] = tr[0].data
    return ss


def PrintCC(dat1,dat2,Pos,CCC,figid, subpid):
    dat1 = dat1 / np.max(dat1)
    dat2 = dat2 / np.max(dat2)
    fig0 = plt.figure(figid)
    tt1 =  np.arange(0,len(dat1),1)
    tt2 =  np.arange(Pos,Pos+len(dat2),1)
    ax = fig0.add_subplot(8,1,subpid)
    ax.plot(tt1,dat1)
    ax.plot(tt2,dat2)
    ax.set_title('CCC = %f   lag = %d' % (CCC, Pos))


def Correlate3(ss,gg,PlotCC):
    '''
    Using the correlation from the maximum value between the three channels
    :param ss:
    :param gg:
    :param PlotCC:
    :return:
    '''
    #
    [n1,m1,k1] = ss.shape
    [n2,m2,k2] = gg.shape
    CClag = np.zeros(n1)
    CCC = np.zeros(n1)
    wid = int(k1/4)

    for ii in range(n1):

        dataT = ss[ii][0][:]
        dataR = ss[ii][1][:]
        dataZ = ss[ii][2][:]

        maxAmpT = np.max(np.abs(ss[ii][0][:]))
        maxAmpR = np.max(np.abs(ss[ii][1][:]))
        maxAmpZ = np.max(np.abs(ss[ii][2][:]))

        maxChan = np.argmax([maxAmpT, maxAmpR, maxAmpZ])
        # maxChan = 2


        data_u1 = gg[ii][0][:]
        data_u2 = gg[ii][1][:]
        data_u3 = gg[ii][2][:]
        data_u4 = gg[ii][3][:]
        data_u5 = gg[ii][4][:]
        data_u6 = gg[ii][5][:]
        data_u7 = gg[ii][6][:]
        data_u8 = gg[ii][7][:]

        x = np.arange(16.).reshape(8,2)
        t = np.arange(4.).reshape(2,2)
        r = np.arange(6.).reshape(3,2)
        v = np.arange(6.).reshape(3,2)

        # Tangential

        a,b = xcorr(dataT, data_u1, wid)
        x[0][0] = b #abs(b)
        x[0][1] = a
        t[0][0] = b #abs(b)
        t[0][1] = a
        a,b = xcorr(dataT, data_u2, wid)
        x[1][0] = b # abs(b)
        x[1][1] = a # a
        t[1][0] = b # abs(b)
        t[1][1] = a # a


        # Radial
        a,b = xcorr(dataR, data_u3, wid)
        x[2][0] = b  #abs(b)
        x[2][1] = a
        r[0][0] = b  # abs(b)
        r[0][1] = a
        a,b = xcorr(dataR, data_u4, wid)
        x[3][0] = b  #abs(b)
        x[3][1] = a
        r[1][0] = b  #abs(b)
        r[1][1] = a
        a,b = xcorr(dataR, data_u5, wid)
        x[4][0] = b  # abs(b)
        x[4][1] = a
        r[2][0] = b  # abs(b)
        r[2][1] = a

        # Vertical
        a,b = xcorr(dataZ, data_u6, wid)
        x[5][0] = b  #abs(b)
        x[5][1] = a
        v[0][0] = b  #abs(b)
        v[0][1]=a
        a,b = xcorr(dataZ, data_u7, wid)
        x[6][0] = b  #abs(b)
        x[6][1] = a
        v[1][0] = b  #abs(b)
        v[1][1] = a
        a,b = xcorr(dataZ, data_u8, wid)
        x[7][0] = b  #abs(b)
        x[7][1] = a
        v[2][0] = b  #abs(b)
        v[2][1] = a

        # sort for zcor
        X = np.array(sorted(sorted(x,key=lambda e:e[1]),key=lambda e:e[0]))
        T = np.array(sorted(sorted(t,key=lambda e:e[1]),key=lambda e:e[0]))
        R = np.array(sorted(sorted(r,key=lambda e:e[1]),key=lambda e:e[0]))
        V = np.array(sorted(sorted(v,key=lambda e:e[1]),key=lambda e:e[0]))
        Zco = X[-1][1]
        Tco = T[-1][1]
        Rco = R[-1][1]
        Vco = V[-1][1]

        if maxChan == 0:
            CClag[ii] = Tco
        elif maxChan == 1:
            CClag[ii] = Rco
        else:
            CClag[ii] = Zco



        if PlotCC == 1:
            PrintCC(dataT, data_u1,Tco, Vco,50+ii,1)
            PrintCC(dataT, data_u2,Tco, Vco,50+ii,2)
            PrintCC(dataR, data_u3,Rco, Vco,50+ii,3)
            PrintCC(dataR, data_u4,Rco, Vco,50+ii,4)
            PrintCC(dataR, data_u5,Rco, Vco,50+ii,5)
            PrintCC(dataZ, data_u6,Zco, Vco,50+ii,6)
            PrintCC(dataZ, data_u7,Zco, Vco,50+ii,7)
            PrintCC(dataZ, data_u8,Zco, Vco,50+ii,8)
            plt.show()


    return CClag


def Dist2Weight(DIST):
    mindist=100000.0
    ll = len(DIST)
    cormax = np.zeros(ll)
    for ii in range(ll):
        cormax[ii] = DIST[ii] / mindist
    return cormax


def column(matrix, i):
   return [row[i] for row in matrix]


def minvdbl(A,v):

    v=column(v, 0)
    N = len(v)

    # Gaussian elimination
    for m in range(N):

        # Divide by the diagonal element
        div = A[m,m]
        A[m,:] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1,N):
            mult = A[i,m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]

    # Backsubstitution

    #define length of x
    x = np.zeros(N+1,float)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m,i]*x[i]

    return x

def toRtf(k,M):

    Aki = np.zeros(shape=(3,3))
    Aki[0][0] = ( 1.0* k * M[0])
    Aki[0][1] = (-1.0* k * M[2])
    Aki[0][2] = ( 1.0* k * M[3])
    Aki[1][0] = (-1.0* k * M[2])
    Aki[1][1] = ( 1.0* k * M[1])
    Aki[1][2] = (-1.0* k * M[4])
    Aki[2][0] = ( 1.0* k * M[3])
    Aki[2][1] = (-1.0* k * M[4])
    Aki[2][2] = ( 1.0* k * M[5])

    return Aki

def toXyz(k,M):

    Aki = np.zeros(shape=(3,3))
    Aki[0][0] = ( 1.0* k * M[0])
    Aki[0][1] = ( 1.0* k * M[2])
    Aki[0][2] = ( 1.0* k * M[3])
    Aki[1][0] = ( 1.0* k * M[2])
    Aki[1][1] = ( 1.0* k * M[1])
    Aki[1][2] = ( 1.0* k * M[4])
    Aki[2][0] = ( 1.0* k * M[3])
    Aki[2][1] = ( 1.0* k * M[4])
    Aki[2][2] = ( 1.0* k * M[5])
    return Aki

def mt2parameters(MTx,M,gfscale):

    #Iso Mo
    MoIso = (MTx[0][0] + MTx[1][1] + MTx[2][2])/3

    c = MomentTensor(MTx[2][2],MTx[0][0],MTx[1][1],MTx[0][2],-MTx[1][2],-MTx[0][1],gfscale)

    # Principal axes
    (T, N, P) = mt2axes(c)

    # Nodal planes
    np0 = mt2plane(c)
    np2 = aux_plane(np0.strike,np0.dip,np0.rake)
    # Convention rake: up-down
    if (np0.rake>180):
       np0.rake= np0.rake-360
    if (np0.rake<-180):
       np0.rake= np0.rake+360
    np1 = np.zeros(3)
    np1 = [np0.strike,np0.dip,np0.rake]

    # Compute Eigenvectors and Eigenvalues
    # Seismic Moment and Moment Magnitude
    (EigVal, EigVec) = np.linalg.eig(MTx)
    b = copy.deepcopy(EigVal)
    b.sort()
    Mo = (abs(b[0]) + abs(b[2]))/2.
    Mw = np.log10(Mo)/1.5-10.73

    # Compute double-couple, CLVD & iso
    d = copy.deepcopy(EigVal)
    d[0]=abs(d[0])
    d[1]=abs(d[1])
    d[2]=abs(d[2])
    d.sort()
    eps=abs(d[0])/abs(d[2])
    pcdc=100.0*(1.0-2.0*eps)
    pcclvd=200.0*eps
    pcdc=pcdc/100.0
    pcclvd=pcclvd/100.0
    pciso=abs(MoIso)/Mo
    pcsum=pcdc+pcclvd+pciso
    pcdc=100.0*pcdc/pcsum
    pcclvd=100.0*pcclvd/pcsum
    pciso=100.0*pciso/pcsum

    Pdc   = pcdc
    Pclvd = pcclvd

    return (Mo, Mw, Pdc, Pclvd, pciso, EigVal, EigVec, T, N, P, np.round(np1), np.round(np2), c)



def MakeSynth(mti, gg, BAZ):

    shp = gg.shape
    N_STA = shp[0]
    ngrn = shp[1]
    Np = shp[2]
    sg = np.zeros((N_STA, 3, Np))

    mxx = mti[0]
    myy = mti[1]
    mxy = mti[2]
    mxz = mti[3]
    myz = mti[4]
    mzz = mti[5]

    for i in range(N_STA):

        Az = np.deg2rad(BAZ[i])

        for j in range(Np):
            # TAN Component
            sg[i][0][j] = \
                   (mxx*0.5*gg[i][0][j]*np.sin(2*Az) \
                    - myy*0.5*gg[i][0][j]*np.sin(2*Az) \
                    - mxy*gg[i][0][j]*np.cos(2*Az) \
                    - mxz*gg[i][1][j]*np.sin(Az) \
                    + myz*gg[i][1][j]*np.cos(Az)) *-1

            # RAD Component
            sg[i][1][j] = \
                   (mxx*1/6*gg[i][4][j] \
                    - mxx*0.5*gg[i][2][j]*np.cos(2*Az) \
                    + mxx*1/3*gg[i][8][j] \
                    + myy*1/6*gg[i][4][j] \
                    + myy*0.5*gg[i][2][j]*np.cos(2*Az) \
                    + myy*1/3*gg[i][8][j] \
                    + mzz*1/3*gg[i][8][j] \
                    - mzz*1/3*gg[i][4][j] \
                    - mxy*gg[i][2][j]*np.sin(2*Az) \
                    + mxz*gg[i][3][j]*np.cos(Az) \
                    + myz*gg[i][3][j]*np.sin(Az)) *-1

            # VER Component
            sg[i][2][j] = \
                   (mxx*1/6*gg[i][7][j] \
                    - mxx*0.5*gg[i][5][j]*np.cos(2*Az) \
                    + mxx*1/3*gg[i][9][j] \
                    + myy*1/6*gg[i][7][j] \
                    + myy*0.5*gg[i][5][j]*np.cos(2*Az) \
                    + myy*1/3*gg[i][9][j] \
                    + mzz*1/3*gg[i][9][j] \
                    - mzz*1/3*gg[i][7][j] \
                    - mxy*gg[i][5][j]*np.sin(2*Az) \
                    + mxz*gg[i][6][j]*np.cos(Az) \
                    + myz*gg[i][6][j]*np.sin(Az)) *-1
    return sg


def plotSynth(ss,sg,Zcor):
    N_STA = len(Zcor)
    fig1 = plt.figure(44)



    for ii in range(N_STA):

        maxVal1 = np.zeros(3)
        maxVal2 = np.zeros(3)
        for jj in range(3):
            maxVal1[jj] = np.max(np.abs(ss[ii][jj][:]))
            maxVal2[jj] = np.max(np.abs(sg[ii][jj][:]))
        maxVal1 = np.max(maxVal1)
        maxVal2 = np.max(maxVal2)
        maxVal = np.max([maxVal1, maxVal2])

        txss = np.arange(0,len(ss[ii][0][:]),1) - Zcor[ii]
        txsg = np.arange(0,len(sg[ii][0][:]),1)

        ax1 = fig1.add_subplot(N_STA, 3, 1 + 3*ii)
        ax1.plot(txss, ss[ii][0][:],'-k')
        ax1.plot(txsg, sg[ii][0][:],'--r')
        ax1.set_ylim([-maxVal, maxVal])
        ax1.set_xlim([0, len(ss[ii][0][:]) - Zcor[ii]])

        ax1 = fig1.add_subplot(N_STA, 3, 2 + 3*ii)
        ax1.plot(txss, ss[ii][1][:],'-k')
        ax1.plot(txsg, sg[ii][1][:],'--r')
        ax1.set_ylim([-maxVal, maxVal])
        ax1.set_xlim([0, len(ss[ii][0][:]) - Zcor[ii]])

        ax1 = fig1.add_subplot(N_STA, 3, 3 + 3*ii)
        ax1.plot(txss, ss[ii][2][:],'-k')
        ax1.plot(txsg, sg[ii][2][:],'--r')
        ax1.set_ylim([-maxVal, maxVal])
        ax1.set_xlim([0, len(ss[ii][0][:]) - Zcor[ii]])

    plt.show()


def fitcheck(ss,sg,W,N_STA,isoflag,Zcor):

    # mo  = 1.0*mo
    # MM  = mo/(1e20)
    # MM  = 1.0
    # Mscl = Mo/1.0e+20
    WSUM=VAR=DVAR=Dtot=Etot=0.0
    cnt = 0
    VRstn = np.zeros(N_STA)
    for i in range(N_STA):
        Dpower=0.0
        E     =0.0
        Zg=int(Zcor[i])
        Np=len(ss[i][0][:])
        for j in range(Np):
            if Zg+j < Np:

                Etmp = ss[i][0][Zg+j] - sg[i][0][j]
                E += Etmp*Etmp


                Etmp = ss[i][1][Zg+j] - sg[i][1][j]
                E += Etmp*Etmp


                Etmp = ss[i][2][Zg+j] - sg[i][2][j]
                E += Etmp*Etmp


                Dpower += ss[i][0][Zg+j]**2
                Dpower += ss[i][1][Zg+j]**2
                Dpower += ss[i][2][Zg+j]**2
                cnt = cnt + 1
        WSUM += W[i] # 3*cnt-1]
        Etot += E
        VAR += W[i]*E # 3*cnt-1]*E
        Dtot += Dpower
        DVAR += W[i]*Dpower # 3*cnt-1]*Dpower
        E /= Dpower
        VRstn[i] = (1.0 - E)*100.0
        # if VRstn[i] < 0:
        #     print(VRstn[i])

    var0 = Etot/(3.0*float(cnt) - float(isoflag) - 1.0)
    Etot /= Dtot
    vred = (1.0-Etot)*100.0
    VAR /= WSUM
    DVAR /= WSUM
    VAR /= DVAR
    VAR = (1.0-VAR)*100.0

    return VRstn, VAR


def SetQuality(VAR):
    # Set UQality
    if (VAR < 20.0):
       qual = 0
    elif (VAR >= 20.0 and VAR < 40.0):
       qual = 1
    elif (VAR >= 40.0 and VAR < 60.0):
       qual = 2
    elif (VAR >= 60.0 and VAR < 80.0):
       qual = 3
    else:
       qual = 4
    return qual

def vol_p2isoflag(vol_p):
    if vol_p == 0:
        isoflag = 5
    else:
        isoflag = 6
    return isoflag

def MakeWeight(N_STN, Np, W0):
    WCNT = 3*(N_STN * Np)
    W = np.ones(WCNT)
    l = 0
    if N_STN > 1:
        for ii in range(N_STN):
            for jj in range(3 * Np):
                W[l] = W0[ii]
                l = l + 1
    return W, WCNT

def MakeAIV(N_STN, isoflag, ZCORS, ZcorX, BAZ, gg, WCNT, Np, W):
    ZCOR = np.zeros(N_STN)
    for ii in range(N_STN):
        if ZCORS[ii] == 0:
            ZCOR[ii] = ZcorX[ii]
        else:
            ZCOR[ii] = ZCORS[ii]

    AIV = np.zeros((isoflag,isoflag))
    AJ = np.zeros((isoflag,WCNT))
    cnt1 = 0
    cnt2 = 0
    cnt3 = 0

    for ii in range(N_STN):
        BAZR = float(BAZ[ii]) * (np.pi / 180.0)
        # Z = ZCOR[ii]
        Z = 0
        cnt1 = cnt2 = cnt3
        cnt2 += Np
        cnt3 += 2*Np
        for j in range(Np):
            jj = int(j + Z)
            AJ[0][cnt1] = 0.5 * np.sin(2*BAZR)*gg[ii][0][jj]
            if isoflag == 6:
                AJ[0][cnt2] = 1/6*gg[ii][4][jj] - 0.5*np.cos(2.*BAZR)*gg[ii][2][jj] + 1/3*gg[ii][8][jj]
                AJ[0][cnt3]	= 1/6*gg[ii][7][jj] - 0.5*np.cos(2.*BAZR)*gg[ii][5][jj] + 1/3*gg[ii][9][jj]

            if isoflag==5:
                AJ[0][cnt2] = 0.5*gg[ii][4][jj] - 0.5*np.cos(2.*BAZR)*gg[ii][2][jj]
                AJ[0][cnt3]	= 0.5*gg[ii][7][jj] - 0.5*np.cos(2.*BAZR)*gg[ii][5][jj]

            AJ[1][cnt1] =(-0.5)*np.sin(2.*BAZR)*gg[ii][0][jj]
            if isoflag==6:
                AJ[1][cnt2] = 1/6*gg[ii][4][jj] + 0.5*np.cos(2.*BAZR)*gg[ii][2][jj] + 1/3*gg[ii][8][jj]
                AJ[1][cnt3] = 1/6*gg[ii][7][jj] + 0.5*np.cos(2.*BAZR)*gg[ii][5][jj] + 1/3*gg[ii][9][jj]

            if isoflag==5:
                AJ[1][cnt2]	= 0.5*gg[ii][4][jj] + 0.5*np.cos(2.*BAZR)*gg[ii][2][jj]
                AJ[1][cnt3]	= 0.5*gg[ii][7][jj] + 0.5*np.cos(2.*BAZR)*gg[ii][5][jj]

            AJ[2][cnt1] = (-1.0)*np.cos(2.*BAZR)*gg[ii][0][jj]
            AJ[2][cnt2] = (-1.0)*np.sin(2.*BAZR)*gg[ii][2][jj]
            AJ[2][cnt3] = (-1.0)*np.sin(2.*BAZR)*gg[ii][5][jj]

            AJ[3][cnt1] = (-1.0)*np.sin(BAZR)*gg[ii][1][jj]
            AJ[3][cnt2] = np.cos(BAZR)*gg[ii][3][jj]
            AJ[3][cnt3] = np.cos(BAZR)*gg[ii][6][jj]

            AJ[4][cnt1] = np.cos(BAZR)*gg[ii][1][jj]
            AJ[4][cnt2] = np.sin(BAZR)*gg[ii][3][jj]
            AJ[4][cnt3] = np.sin(BAZR)*gg[ii][6][jj]

            if isoflag==6:
                AJ[5][cnt1] = 0.0
                AJ[5][cnt2] = 1/3*gg[ii][8][jj]-1/3*gg[ii][4][jj]
                AJ[5][cnt3] = 1/3*gg[ii][9][jj]-1/3*gg[ii][7][jj]

            cnt1 = cnt1+1
            cnt2 = cnt2+1
            cnt3 = cnt3+1


    for ii in range(isoflag):
        for jj in range(isoflag):
            for kk in range(cnt3):
                AIV[ii][jj] += AJ[ii][kk]* AJ[jj][kk] * W[kk]

    return AIV, ZCOR, AJ

def MakeB(isoflag, N_STN, Np, ss, AJ, ZCOR, gg, W):

    B = np.zeros((isoflag, N_STN))
    nn = len(gg[0][0][:])

    cnt1=0
    cnt2=0
    cnt3=0
    tmp=np.zeros(10 * Np * nn)
    for jj in range(N_STN):
        Z = ZCOR[jj]
        cnt1=cnt2 = cnt3
        cnt2 += Np
        cnt3 += 2*Np
        for ii in range(Np):
            i = int(ii+Z)
            if i < Np:
                tmp[cnt1] = ss[jj][0][i]
                tmp[cnt2] = ss[jj][1][i]
                tmp[cnt3] = ss[jj][2][i]
            else:
                tmp[cnt1] = 0
                tmp[cnt2] = 0
                tmp[cnt3] = 0

            cnt1 = cnt1 + 1
            cnt2 = cnt2 + 1
            cnt3 = cnt3 + 1

    for ii in range(isoflag):
        for jj in range(cnt3):
            B[ii][0] += AJ[ii][jj] * tmp[jj] * W[jj]
    return B

def tdmtw_invc_iso(ss,gg,ZcorX,ZCORS,W0, isoflag,BAZ):

    N_STN = len(BAZ)
    Np = len(ss[0][0][:])

    W, WCNT = MakeWeight(N_STN, Np, W0)

    # Set and normalize AtA and AIV matrix
    AIV, ZCOR, AJ = MakeAIV(N_STN, isoflag, ZCORS, ZcorX, BAZ, gg, WCNT, Np, W)

    # Calculate Righthand Side
    B = MakeB(isoflag, N_STN, Np, ss, AJ, ZCOR, gg, W)

    # Solve for MT
    M = minvdbl(AIV,B)

    # gf scaling factor
    gfscale = 1.0e+20 # Dyn * cm

    # if isotropic constrain, set Mzz
    if(isoflag == 5):
       zz = -1*(M[0]+M[1])
       M[5]=zz

    M *= -1.

    # *Convert deviatoric moment tensor to AKI convention*
    MTx = toXyz(gfscale,M)
    MTr = toRtf(gfscale,M)

    return MTx, M, gfscale, ZCOR


def make360(a):
    if a < 0:
        a = a + 360
    if a > 360:
        a = a - 360
    if a < 0:
        a = a + 360
    return a

def GetGRN(Path2_Green, N_STN, depth, DIST):
    ListGF = glob.glob(Path2_Green+"*d%s.MSEED" % depth)
    if len(ListGF) == 0:
        sys.exit("No available Green Functions for depth: %d  " % depth)
    GRN = list()
    for ii in range(N_STN):
        DistGF = np.zeros(len(ListGF))
        for jj in range(len(ListGF)):
            for kk in range(len(ListGF[jj])):
                if ListGF[jj][kk] == "d":
                    posD2 = kk
                elif ListGF[jj][kk] == "_":
                    posD1 = kk
            DistGF[jj] = ListGF[jj][posD1+1:posD2]
        Opt_Dist = abs(DistGF - DIST[ii])
        GRN += {ListGF[np.argmin(Opt_Dist)]}
    return GRN


def MakeInfile(GRN, MT_F_name, vol_p, N_STN, depth, ZCORS, STN_LISTN, DIST, BAZ,Sampl4inv):
    mtf = open(MT_F_name, 'w')
    if vol_p == 0:
        mtf.write('%s %s 1 5 1\n' %(N_STN, depth))
    else:
        mtf.write('%s %s 1 6 1\n' %(N_STN, depth))

    for ii in range(N_STN):

        if N_STN > 1 and len(ZCORS) > 1:
            Zcor = ZCORS[ii]
        else:
            Zcor = ZCORS[0]

        mtf.write('%s %s %s %d %d\n' % (STN_LISTN[ii], DIST[ii], BAZ[ii], Zcor, Sampl4inv))

    for ii in range(len(GRN)):
        mtf.write('%s\n' % GRN[ii])
    f_out1 = 'mt_inv.d%d.out' % depth
    f_out2 = 'mt_inv_out_d%d' % depth
    mtf.write('%s\n' % f_out2)
    mtf.close()
    return f_out1, f_out2


def RunSysIsoInv(pathdir, f_out1, path_2_iso, f_out2, pathOut):
    if os.path.exists(pathdir+"/"+f_out1):
        os.remove(pathdir+"/"+f_out1)
    os.system(path_2_iso+"tdmt_invc_iso")
    os.rename(pathdir+"/mt_inv.out", pathOut+"/"+f_out1)
    os.rename(pathdir+"/"+f_out2, pathOut+"/"+f_out2)


def GetCatSolJ(DEPTH, Depth0):
    Vz = np.abs(DEPTH - Depth0)
    pos = np.argmin(Vz)
    depth1 = DEPTH[pos]
    depth1d = Vz[pos]
    if depth1d > 2:
        print('Solution for Catalog depth is solved for depth=%d' % depth1)
    return pos

def NrowsP(nsta):
    # Plot py MT
    if nsta < 5:
        nRows = 5
    else:
        nRows = nsta
    return nRows

def PlotMomentTensor(idfig, DEPTH, Depth0, VR, Sampl4inv, FMS,
                     MWz, data_T, data_R, data_Z,
                     synth_T, synth_R, synth_Z,
                     stn_F, dist_F, Az_F, Zcor_F, npts_F, VR_F, dt_F, MT,DIST,
                     Lat0, Long0, LON_STN, LAT_STN, faults1,
                     Green, Event, EventN, pathplt, GFtype,
                     MoL, pdc, pclvd, piso, MODE, AUTOMODE, Fmin, Fmax):

    # Get catalog solution
    pos_Depth0 = GetCatSolJ(DEPTH, Depth0)
    pos_Mw = np.argmax(VR)
    nsta = len(stn_F)

    fig1000 = plt.figure(idfig+pos_Mw, figsize=(23.5, 13.0), dpi=50)
    posSTN_txt = int(Sampl4inv * 0.3)
    Rake_Final1 = FMS[pos_Mw, 2]
    Mw0_final = MWz[pos_Mw]

    iFault = FMS2fType(Rake_Final1)
    L_wc, W_wc = iFault2LW(iFault, Mw0_final)
    A_r = Mw2a(Mw0_final)
    nRows = NrowsP(nsta)

    for jj in range(nsta):
        amp_Max = max([max(abs(data_T[:, jj,pos_Mw])), max(abs(data_R[:, jj, pos_Mw])), max(abs(data_Z[:, jj, pos_Mw]))])
        # Tangential
        ax1 = fig1000.add_subplot(nRows,6,jj*6+1)
        ax1.plot(data_T[:, jj, pos_Mw],'-k')
        ax1.plot(synth_T[:, jj, pos_Mw],'--r')
        ax1.plot(synth_T[:, jj, pos_Depth0],'--g')
        ax1.set_ylim([-amp_Max, amp_Max])
        ax1.set_axis_off()
        ax1.text(-(Sampl4inv/5),0,stn_F[jj], fontsize=14)
        #  postxtY = min([min(data_T[:, jj, pos_Mw]), min(synth_T[:, jj, pos_Mw]), min(synth_R[:, jj, pos_Mw]), min(synth_Z[:, jj, pos_Mw]), min(data_Z[:, jj, pos_Mw])])
        #  postxtY = min([min(data_T[:, jj, pos_Mw]), min(synth_T[:, jj, pos_Mw])])
        ax1.text(0,amp_Max,'Dist = %dkm, Azimuth = %d, Zcor = %d, VR = %d'
                 % (dist_F[jj], Az_F[jj], Zcor_F[jj,pos_Mw], VR_F[jj,pos_Mw]))

        if jj == 0:
            # ax1.set_title('Tangential')
            ax1.text(0,amp_Max*1.2, 'Tangential', fontsize=16)
        if jj == nsta-1:
            font = font0.copy()
            font.set_style('italic')
            font.set_weight('bold')
            # ax1.text(-posSTN_txt,-amp_Max*1.5, 'Waveform data (solid line) and synthetic data (dashed line: Red - max(VR), Green - cat. depth) from the moment tensor inversion',fontproperties=font)
            ax1.text(-posSTN_txt,-amp_Max*1.5,
                     'Waveform data (solid line) and synthetic data (dashed red line) '
                     'from the moment tensor inversion',fontproperties=font)

        # Radial
        ax1 = fig1000.add_subplot(nRows,6,jj*6+2)
        ax1.plot(data_R[:, jj, pos_Mw],'-k')
        ax1.plot(synth_R[:, jj, pos_Mw],'--r')
        ax1.plot(synth_R[:, jj, pos_Depth0],'--g')
        ax1.set_ylim([-amp_Max, amp_Max])
        ax1.set_axis_off()
        if jj == 0:
            # ax1.set_title('Radial')
            ax1.text(0,amp_Max*1.2, 'Radial', fontsize=14)

        # Vertical
        ax1 = fig1000.add_subplot(nRows,6,jj*6+3)
        ax1.plot(data_Z[:, jj, pos_Mw],'-k')
        ax1.plot(synth_Z[:, jj, pos_Mw],'--r')
        ax1.plot(synth_Z[:, jj, pos_Depth0],'--g')
        ax1.set_ylim([-amp_Max, amp_Max])
        ax1.set_axis_off()

        # set scale bar:
        Ybar = amp_Max/2.0
        dybar = amp_Max/10.0
        barp = Sampl4inv/4
        Xbar1 = npts_F-int(barp)
        Xbar2 = npts_F
        ax1.plot([Xbar1, Xbar2], [-Ybar, -Ybar],'-b')
        ax1.plot([Xbar1, Xbar1],[-Ybar-dybar, -Ybar+dybar],'-b')
        ax1.plot([Xbar2, Xbar2],[-Ybar-dybar, -Ybar+dybar],'-b')
        ax1.text(Xbar1+1, -Ybar-(2*dybar), '%s sec.' % round(barp*dt_F,2))


        if jj == 0:
            #ax1.set_title('Vertical')
            ax1.text(0,amp_Max*1.2, 'Vertical', fontsize=14)


    LR = 1.5
    ax1 = fig1000.add_subplot(3,6,4) # Plot Beach-ball according to the best VR
    b1 = beach(FMS[pos_Mw, 0:3], xy=(0,0.2), width=LR,linewidth=1.1, facecolor='r', nofill=True , edgecolor='r')
    b2 = beach(MT[pos_Mw, 0:6], xy=(0,0.2), width=LR,linewidth=0.2, facecolor='k')
    ax1.add_collection(b2)
    ax1.set_ylim([-1,1])
    ax1.set_xlim([-1,1])
    ax1.set_aspect("equal")
    ax1.set_axis_off()
    ax1.add_collection(b2)
    ax1.add_collection(b1)
    ax1.set_title('Solution Max VR: Depth= %d km, VR= %d'
                  % (DEPTH[pos_Mw], np.round(VR[pos_Mw])),family="monospace", weight="bold", color="r", fontsize=14)


    LR = 1.1
    ax1 = fig1000.add_subplot(3,6,10) # Plot Beach-ball according to the Catalog depth
    b1 = beach(FMS[pos_Depth0, 0:3], xy=(0,0.2), width=LR,linewidth=1.1, facecolor='g', nofill=True , edgecolor='g')
    b2 = beach(MT[pos_Depth0, 0:6], xy=(0,0.2), width=LR,linewidth=0.2, facecolor='k')
    ax1.add_collection(b2)
    ax1.set_ylim([-1,1])
    ax1.set_xlim([-1,1])
    ax1.set_aspect("equal")
    ax1.set_axis_off()
    ax1.add_collection(b2)
    ax1.add_collection(b1)
    ax1.set_title('Cat Solution: Depth= %d km, VR= %d'
                  % (DEPTH[pos_Depth0], np.round(VR[pos_Depth0])),family="monospace", weight="bold", color="g", fontsize=14)


    ax1 = fig1000.add_subplot(3,6,16)

    ax1.text(0,1.0,'GF file: %s\nM11 = %5.3e\nM22 = %5.3e\nM33 = %5.3e\nM12 = %5.3e\nM13 = %5.3e\nM23 = %5.3e\n'\
             % (Green, MT[pos_Mw,0], MT[pos_Mw,1], MT[pos_Mw,2], MT[pos_Mw,3], MT[pos_Mw,4], MT[pos_Mw,5]))

    ax1.text(0,0.4,'Centriod Depth = %s km \nStrike (VR)= %d ; %d \nDip (VR)= %d ; %d \nRake (VR)= %d ; %d \n'
                   'Mo (VR)= %5.2e [dyn/cm]\nMw (VR)= %4.2f \nPrecent DC (VR)= %5.2f\n'
                   'Precent CLVD (VR)= %5.2f\nPrecent ISO (VR)= %5.2f\nMax. Var. Red. (VR)= %5.2f\n' \
             % (DEPTH[pos_Mw], FMS[pos_Mw, 0], FMS[pos_Mw, 3], FMS[pos_Mw, 1], FMS[pos_Mw, 4], FMS[pos_Mw, 2],
                FMS[pos_Mw, 5], MoL, MWz[pos_Mw], float(pdc[pos_Mw]), float(pclvd[pos_Mw]),
                float(piso[pos_Mw]), float(VR[pos_Mw])), color="r")

    ax1.text(0,-0.2,'Strike (D)= %d ; %d \nDip (D)= %d ; %d \nRake (D)= %d ; '
                    '%s \nMo (D)= %5.2e [dyn/cm]\nMw (D)= %4.2f \nPrecent DC (D)= %5.2f'
                    '\nPrecent CLVD (D)= %5.2f\nPrecent ISO (D)= %5.2f\nMax. Var. Red. (D)= %5.2f\n' \
             % (FMS[pos_Depth0, 0], FMS[pos_Depth0, 3], FMS[pos_Depth0, 1], FMS[pos_Depth0, 4],
                FMS[pos_Depth0, 2], FMS[pos_Depth0, 5], MoL, MWz[pos_Depth0], pdc[pos_Depth0],
                pclvd[pos_Depth0], piso[pos_Depth0], VR[pos_Depth0]),color="g")

    ax1.text(0,-0.4, 'FAULT LENGTH (W & C): %6.1f \nFAULT WIDTH (W & C): %6.1f \nFault length Crack: %6.1f\n' \
             % (L_wc, W_wc, A_r), color="b")

    ax1.set_axis_off()             #'RES/Pdc'



    ax = fig1000.add_subplot(2,3,3)
    for zz in range(len(DEPTH)-1):
        try:
            b = beach(MT[zz, 0:6], xy=(DEPTH[zz], VR[zz]), width=1,linewidth=0.1,facecolor='k')
            # b = beach(FMS[zz, 0:3], xy=(DEPTH[zz], VR[zz]), width=1.0,linewidth=0.1, facecolor='k')
            ax.add_collection(b)
        except:
            continue
    # b1 = beach(FMS[pos_Depth0, 0:3], xy=(DEPTH[pos_Depth0], VR[pos_Depth0]), width=1.0,linewidth=0.1, facecolor='g')
    b1 = beach(MT[pos_Depth0, 0:6], xy=(DEPTH[pos_Depth0], VR[pos_Depth0]), width=1.0,linewidth=0.1, facecolor='g')
    ax.add_collection(b1)
    # b2 = beach(FMS[pos_Mw, 0:3], xy=(DEPTH[pos_Mw], VR[pos_Mw]), width=1.0,linewidth=0.1, facecolor='r')
    b2 = beach(MT[pos_Mw, 0:6], xy=(DEPTH[pos_Mw], VR[pos_Mw]), width=1.0,linewidth=0.1, facecolor='r')
    ax.add_collection(b2)
    Mw = MWz[pos_Mw]
    if 'FMS0' in globals():
        # b = beach(FMS0[0:3], xy = (DEPTH[pos_Mw],VR[pos_Mw] - 2), width=1.0,linewidth=0.1, facecolor='black')
        b = beach(FMS0[0:3], xy=(Depth0, max(VR) - 2), width=1.0,linewidth=0.1, facecolor='black')
        # b = beach(MT[pos_Mw, 0:6], xy=(Depth0, max(VR) - 2), width=1.0,linewidth=0.1, facecolor='black')
        ax.add_collection(b)


    ax.set_ylim(np.min(VR)-1, np.max(VR)+1)
    ax.set_xlabel('Depth [km]')
    ax.set_ylabel('Variance reduction [%]')

    ax.set_xlim(np.min(DEPTH)-0.1, np.max(DEPTH)+0.1)
    n_ticks = 5
    dx_tick = ((np.max(DEPTH)+5) - (np.min(DEPTH)-5))/n_ticks
    dy_tick = ((np.ceil(np.max(VR))+5) - (np.floor(np.min(VR))-5))/n_ticks
    plt.xticks(np.arange(np.min(DEPTH)-5, np.max(DEPTH)+5, dx_tick))
    plt.yticks(np.arange(np.floor(np.min(VR)) - 5, np.ceil(np.max(VR)) + 5, 5))
    ax.set_aspect("equal")


    plt.subplot(2,3,6)
    dll = (max(DIST) * 1.3) * 0.0090 # km to deg
    llat = Lat0 - dll # 29.0
    ulat = Lat0 + dll #34.0
    llon = Long0 - dll #32.0
    ulon = Long0 + dll # 37.0

    plot_map_basemap = 0

    Long01 = Long0 + np.cos(45) * (dll/2.5)
    Lat01 = Lat0 + np.sin(45) * (dll/2.5)

    if plot_map_basemap == 1:
        m = Basemap(projection='cyl', lon_0=Long0, lat_0=Lat0,
                    llcrnrlon=llon,llcrnrlat=llat,urcrnrlon=ulon,urcrnrlat=ulat, resolution='i')
        m.fillcontinents(color='wheat', lake_color='skyblue')
        m.drawmapboundary(fill_color='skyblue')
        m.drawparallels(np.arange(llat, ulat, (ulat - llat) / 4.0), labels=[1,0,0,0],fontsize=10, fmt="%.2f")
        m.drawmeridians(np.arange(llon, ulon, (ulon - llon) / 4.0), labels=[0,0,0,1],fontsize=10, fmt="%.2f")
        m.drawcoastlines()

        ax1 = plt.gca()
        m.scatter(Long0, Lat0, 50, color="r", marker=(5, 1), edgecolor="k")
        # b = beach(FMS[pos_Mw, 0:3], xy=(Long0, Lat0),width=dll/10, linewidth=0.1, facecolor='b')
        b = beach(MT[pos_Mw, 0:6], xy=(Long0, Lat0),width=dll/10, linewidth=0.1, facecolor='k')
        ax1.add_collection(b)
        x0, y0 = m(LON_STN, LAT_STN)
        m.scatter(x0, y0, 50, color="g", marker="v", edgecolor="k")
        for i in range(nsta):
            plt.text(x0[i], y0[i], stn_F[i], va="top", family="monospace", weight="bold")

        for shape in faults1.shapeRecords():
            x = [i[0] for i in shape.shape.points[:]]
            y = [i[1] for i in shape.shape.points[:]]
            plt.plot(x,y, color="k", linewidth=0.8)
    else:
        ax1 = plt.gca()
        for shape in faults1.shapeRecords():
            x = [i[0] for i in shape.shape.points[:]]
            y = [i[1] for i in shape.shape.points[:]]
            ax1.plot(x,y, color="k", linewidth=0.8)
        ax1.scatter(LON_STN,LAT_STN,s=150,marker='v', color='g')
        for i in range(len(LON_STN)):
            ax1.text(LON_STN[i], LAT_STN[i], stn_F[i], va="top", family="monospace", weight="bold")
        ax1.plot([Long0, Long01], [Lat0, Lat01], color='k')
        ax1.set_xlim([llon, ulon])
        ax1.set_ylim([llat, ulat])
        # b = beach(FMS[pos_Mw, 0:3], xy=(Long01, Lat01),width=dll/5, linewidth=0.1, facecolor='b')
        b = beach(MT[pos_Mw, 0:6], xy=(Long01, Lat01),width=dll/5, linewidth=0.1, facecolor='k')
        ax1.add_collection(b)
        ax1.osm = OSM(ax1)



    ax50 = fig1000.add_subplot(6,11,66)
    ilogo = Image.open(LOGO)
    ax50.imshow(ilogo)
    ax50.set_axis_off()


    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.92, wspace=0.08, hspace=0.17)


    Event = Event[0]
    Event = Event[0]
    Event.replace(':', '_')
    plt.suptitle('%s   Lat: %7.4f   Lon: %7.4f Depth: %d   Mw %3.1f   FMS (s/d/r) %d/%d/%d   %2.2f - %2.2f Hz' %
                 (Event, Lat0, Long0, Depth0, MWz[pos_Mw], FMS[pos_Mw, 0], FMS[pos_Mw, 1], FMS[pos_Mw, 2], Fmin, Fmax), fontsize=20)
    Path2pdf = pathplt+'/%s_d%s_%s.pdf' % (EventN, DEPTH[pos_Mw], GFtype)
    ax1.osm.draw()
    plt.savefig(Path2pdf, dpi=400)

    # Make small figure for FMS display
    plotMiniFMS = 0
    if plotMiniFMS == 1:
        fig1001 = plt.figure(edgecolor=None)
        ax11 = fig1001.add_subplot(1,1,1) # Plot Beach-ball according to the best VR
        b2 = beach(MT[pos_Mw, 0:6], xy=(0,0.2), width=LR,linewidth=0.1, facecolor='k')
        ax11.add_collection(b2)
        ax11.set_ylim([-1,1])
        ax11.set_xlim([-1,1])
        ax11.set_aspect("equal")
        ax11.set_axis_off()
        ax11.add_collection(b2)
        plt.savefig(pathplt+'/%s_d%s_%s_mini.png' % (EventN, DEPTH[pos_Mw], GFtype), dpi=400)
        plt.close(fig1001)


    print('----------- INVERSION ENDED SUCCESSFULLY!!! ------------')
    print('Date-Time  Lon  Lat  Depth ')
    print('%s %7.4f %7.4f %d %d %d %d %3.1f' % (Event, Long0, Lat0, Depth0, FMS[pos_Mw, 0], FMS[pos_Mw, 1], FMS[pos_Mw, 2], MWz[pos_Mw]))

    EVENT = [Event, Event, Long0, Lat0, Depth0, FMS[pos_Mw, 0], FMS[pos_Mw, 1], FMS[pos_Mw, 2], '%2.1f' % MWz[pos_Mw], '%2.1f' % VR[pos_Mw], MODE[AUTOMODE]]
    # ADDeventF(pathFF, EVENT)

    if AUTOMODE == 0:
        plt.show()
    else:
        return Path2pdf
        plt.close('all')

def GetFault(faultsT):
    if faultsT == 1:
        faults1 = shp.Reader(path2Localfaults)
    else:
        faults1 = shp.Reader(PB_boundaries)
    return faults1


def ReadMTfile(N_STN, pathOut, f_out1, f_out2,VR, FMS, MT, pdc, pclvd, piso,
               Zcor_F, VR_F, zz, data_T, data_R, data_Z, synth_T, synth_R, synth_Z, MWz):
    f_MT=open(pathOut+"/"+f_out2)
    lines=f_MT.readlines()
    VR_line     = lines[1].split()
    nsta = int(VR_line[3])

    VR[zz]      = float(VR_line[2])
    FMS_line    = lines[5].split()
    FMS[zz, 0]  = int(FMS_line[0]) # STR 1
    FMS[zz, 2]  = int(FMS_line[1]) # rake 1
    FMS[zz, 1]  = int(FMS_line[2]) # dip 1

    FMS[zz, 3]  = int(FMS_line[3]) # STR 2
    FMS[zz, 5]  = int(FMS_line[4]) # rake 2
    FMS[zz, 4]  = int(FMS_line[5]) # dip 2

    MT_line = lines[3].split()
    MT_linen = np.zeros(len(MT_line))
    for ii in range(len(MT_line)):
        MT_linen[ii] = float(MT_line[ii])
    MT[zz,1] =  MT_linen[0] # mxx M22
    MT[zz,5] = -MT_linen[1] # mxy M23
    MT[zz,3] =  MT_linen[2] # mxz M12
    MT[zz,2] =  MT_linen[3] # myy M33
    MT[zz,4] = -MT_linen[4] # myz M13
    MT[zz,0] =  MT_linen[5] # mzz M11

    DC_line     = lines[7].split()
    pdc[zz] = int(DC_line[0])
    pclvd[zz] = int(DC_line[1])
    piso[zz] = int(DC_line[2])


    nptsn = 0
    stn_F = list()
    dist_F = np.zeros(nsta)

    Az_F = np.zeros(nsta)
    nstai = 0
    for jj in range(len(lines)):
        if lines[jj] == '#filename\n':
            line = lines[jj+1].split()
            stn_F += {line[0]}
            line = lines[jj+3].split()
            dt_F = float(line[0])
            npts_F = int(line[1])
            dist_F[nstai] = float(line[2])
            Az_F[nstai]   = float(line[3])
            Zcor_F[nstai,zz] = int(line[4])
            if line[5] == 'inf':
                VRF = -999
            elif line[5] == '-inf':
                VRF = -999
            else:
                VRF = float(line[5])
            VR_F[nstai,zz] = VRF
            nstai = nstai + 1
            ln0 = jj + 5


            for nn in range(npts_F):
                line = lines[ln0 + nn].split()
                data_T[nn, nptsn, zz] = float(line[0])
                data_R[nn, nptsn, zz] = float(line[1])
                data_Z[nn, nptsn, zz] = float(line[2])

                synth_T[nn, nptsn, zz] = float(line[3])
                synth_R[nn, nptsn, zz] = float(line[4])
                synth_Z[nn, nptsn, zz] = float(line[5])
            nptsn = nptsn + 1



    f_MT    = open(pathOut+"/"+f_out1)
    lines   = f_MT.readlines()
    MWL     = lines[15 + N_STN - 1]
    MWz[zz] = float(MWL[3:7])
    MoL     = lines[14 + N_STN - 1].split()
    MoL = MoL[0]; MoL = float(MoL[3:])

    return VR, FMS, MT, pdc, pclvd, piso, stn_F, dist_F, Az_F, Zcor_F, VR_F, \
           npts_F, dt_F, data_T, data_R, data_Z, synth_T, synth_R, synth_Z, MWz, MoL


def BashPlot(plot_each_ps, path_2_iso, pathOut, f_out2, pathdir, pathplt, DEPTH):
    if plot_each_ps == 1:
        subprocess.run([path_2_iso+"tdmt_plot2", pathOut+"/"+f_out2])
        os.rename(pathdir+"/plot_d%02d_1.ps" % DEPTH, pathplt+"/plot_d%02d_1.ps" % DEPTH)



def MakeHelmSTN(STN_LISTN, SEIS_DISP_ROT_F, pathdir):
    for ii in range(len(STN_LISTN)):
        SEIST = SEIS_DISP_ROT_F.select(station=STN_LISTN[ii]).copy()
        MSEED2HELM(pathdir, STN_LISTN[ii], SEIST)

def AllocateFilesInvPy(l_DEPTH):
    FMS1    = np.zeros((l_DEPTH,3))
    FMS2    = np.zeros((l_DEPTH,3))
    VAR     = np.zeros(l_DEPTH)
    MWzp    = np.zeros(l_DEPTH)

    return FMS1, FMS2, VAR, MWzp


def AllocateFilesInvIso(l_DEPTH, Sampl4inv, N_STN):
    VR    = np.zeros(l_DEPTH)
    MWz   = np.zeros(l_DEPTH)
    FMS   = np.zeros((l_DEPTH,6))
    MT    = np.zeros((l_DEPTH,6))
    pdc   = np.zeros(l_DEPTH)
    pclvd = np.zeros(l_DEPTH)
    piso  = np.zeros(l_DEPTH)


    data_T = np.zeros((Sampl4inv,N_STN, l_DEPTH))
    data_R = np.zeros((Sampl4inv,N_STN, l_DEPTH))
    data_Z = np.zeros((Sampl4inv,N_STN, l_DEPTH))

    synth_T = np.zeros((Sampl4inv,N_STN, l_DEPTH))
    synth_R = np.zeros((Sampl4inv,N_STN, l_DEPTH))
    synth_Z = np.zeros((Sampl4inv,N_STN, l_DEPTH))

    Zcor_F  = np.zeros((N_STN, l_DEPTH))
    VR_F    = np.zeros((N_STN, l_DEPTH))

    return VR, MWz, FMS, MT, pdc, pclvd, piso, data_T, data_R, data_Z, synth_T, synth_R, synth_Z, Zcor_F, VR_F
def sum_gg(gg):
    n1,n2,n3 = gg.shape
    ggN = np.zeros((n1,n3))
    for ii in range(n1):
        for jj in range(n2):
            ggN[ii,:] = ggN[ii,:] + np.abs(gg[ii,jj,:])
        ggN[ii,:] = ggN[ii,:] / np.max(ggN[ii,:])
    return ggN

def calc_ggN_mean(ggN, dmn):
    n1,n2 = ggN.shape
    ggNM = np.zeros((n1,n2))
    for ii in range(n1):
        for jj in range(n2-dmn):
            ggNM[ii,jj] = np.mean(ggN[ii,jj:jj+dmn])
    return ggNM
        
def find_len_taper(ggNM,minA):
    n1,n2 = ggNM.shape
    len_taper = np.zeros(n1)
    first_L = np.zeros(n1)
    for ii in range(n1):
        v = np.arange(0,n2)
        len_taper[ii] = len(v[ggNM[ii,:] > minA])
        first_L[ii] = next(x[0] for x in enumerate(ggNM[ii,:]) if x[1] > minA)
    return len_taper, first_L

def taperer(dat,n,direction):
    if n < len(dat):
        n=int(n)
        if direction == 1:
            dn = int(n/4)
            norm = np.zeros(len(dat))
            norm[0:n-dn] = 1
            norm[n-dn:n] = np.cos(np.linspace(0,np.pi/2,dn))
        elif direction == 2:
            dn = int(n/4)
            norm = np.zeros(len(dat))
            norm[0:n] = 1
            norm[n:n+dn] = np.cos(np.linspace(0,np.pi/2,dn))

        else:
            norm = np.ones(len(dat))
            norm[0:n] = np.cos(np.linspace(-np.pi/2,0,n))
        dat2 = dat * norm
    else:
        dat2 = dat

    return dat2

def taper_ss(ss,len_taper, first_L):
    n1,n2,n3 = ss.shape
    for ii in range(n1):
        for jj in range(n2):
            ss[ii,jj,:] = taperer(ss[ii,jj,:], len_taper[ii],1)
            if first_L[ii] > 0:
                ss[ii,jj,:] = taperer(ss[ii,jj,:], first_L[ii],-1)

    return ss

def MakeMTpyFiles(Zcor, Sampl4inv, zz,sg,ss,MM, N_STN, data_Tp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, MTp, Zcor_Fp, VR_Fp, VRstn):
    MTp[zz][0] = MM.mt[0,0]
    MTp[zz][1] = MM.mt[1,1]
    MTp[zz][2] = MM.mt[2,2]
    MTp[zz][3] = MM.mt[0,1]
    MTp[zz][4] = MM.mt[0,2]
    MTp[zz][5] = MM.mt[2,1]

    for ii in range(N_STN):
        Z = int(Zcor[ii])
        for nn in range(Sampl4inv):
            data_Tp[nn, ii, zz] = ss[ii][0][nn]
            data_Rp[nn, ii, zz] = ss[ii][1][nn]
            data_Zp[nn, ii, zz] = ss[ii][2][nn]
            if nn > Z:
                synth_Tp[nn, ii, zz] = sg[ii][0][nn-Z]
                synth_Rp[nn, ii, zz] = sg[ii][1][nn-Z]
                synth_Zp[nn, ii, zz] = sg[ii][2][nn-Z]

        Zcor_Fp[ii][zz] = Zcor[ii]
        VR_Fp[ii][zz] = VRstn[ii]

    return data_Tp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, MTp, Zcor_Fp, VR_Fp



def TDMTRUN(Depth0, Lat0, Long0, GFtype, Green, STN_LISTN, DIST, BAZ, LAT_STN, LON_STN,
            SEIS_DISP_ROT_F, Sampl4inv, smpl, D1, D2, vol_p, ZCORS, faultsT, AUTOMODE, Event, SEIS0, fcutlow, fcuthigh,Taper2):
    """

    :param Depth0: Catalog depth
    :param Lat0: origine latitude location
    :param Long0: origine longitude location
    :param GFtype: Green function code
    :param Green: Green function model
    :param STN_LISTN: list of selected stations
    :param DIST: distances from earthquake to station
    :param BAZ: azimuth from earthquake to station
    :param LAT_STN: station latitude location
    :param LON_STN: station longitude location
    :param SEIS_DISP_ROT_F: seismic data (displacement, rotated and filtered)
    :param Sampl4inv: length of the inversion
    :param smpl: sample length in seconds
    :param D1: starting depth
    :param D2: ending depth
    :param vol_p: use isotropic component yes/no (1 or 0)
    :param ZCORS: the lag from the cross correlation
    :param faultsT: faults file
    :param AUTOMODE: use automatic mode
    :param Event: event name
    :param SEIS0: seismic data
    :return:
    """

    MODE = 'MA'

    faults1 = GetFault(faultsT)

    pathOut, pathplt, pathdir, MT_F_name, EventF, EventN, pathGF = MKFOLDERS(Event[0], Dir)
    pathplt = MKPLOTDIR(pathplt, STN_LISTN)
    Path2_Green = MT_G + Green

    os.chdir(pathdir)

    # Write SAC files
    write2SAC(SEIS0,pathdir)

    # Convert STN seis to Helm ascii
    MakeHelmSTN(STN_LISTN, SEIS_DISP_ROT_F, pathdir)

    # Set depth range for inversion
    DEPTH = np.arange(int(D1), int(D2), 1)

    l_DEPTH = len(DEPTH)
    N_STN  = len(STN_LISTN)
    Sampl4inv = int(Sampl4inv / smpl)

    Origine_time = UTCDateTime(Event[0][0])

    print('Running Python solver')
    # Allocate for Py
    VRp, MWzp, FMSp, MTp, pdcp, pclvdp, pisop, data_Tp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, Zcor_Fp, VR_Fp = AllocateFilesInvIso(l_DEPTH, Sampl4inv, N_STN)

    for zz in range(l_DEPTH):

        # Obtain the best Green Function from data
        GRN = GetGRN(Path2_Green, N_STN, DEPTH[zz], DIST)

        # make input file for inversion
        f_out1, f_out2 = MakeInfile(GRN, MT_F_name, vol_p, N_STN, DEPTH[zz], ZCORS, STN_LISTN, DIST, BAZ,Sampl4inv)

        isoflag = vol_p2isoflag(vol_p)

        # Make seis files
        ss = Seis2ss(SEIS_DISP_ROT_F, STN_LISTN)

        # Make green function files
        gg = Green2ggMSEED(GRN, vol_p, fcutlow, fcuthigh, smpl)
        
        # make a proxy for taper length from gf
        if Taper2 == 1:
            ggN = sum_gg(gg)
            dmn = int(Sampl4inv / 10)
            ggNM = calc_ggN_mean(ggN, dmn)
            len_taper, first_L = find_len_taper(ggNM,0.01)
            ss = taper_ss(ss,len_taper, first_L)
        
        # Calculate cross correlation
        ZcorX1 = Correlate3(ss,gg,0)

        # Weight by distance
        W = Dist2Weight(DIST)

        # Solve MT
        MTx, M, gfscale, ZcorX2 = tdmtw_invc_iso(ss,gg,ZcorX1,ZCORS,W, isoflag,BAZ)

        # Here compute Planes, Axes, Mo, Mw
        (Mo, MWzp[zz], pdcp[zz], pclvdp[zz], pisop[zz], EigVal, EigVec, T, N, P, FMSp[zz][0:3], FMSp[zz][3:6], MM) = mt2parameters(MTx,M,gfscale)

        # make synt
        sg = MakeSynth(M, gg, BAZ)

        # plot synth + seis
        # plotSynth(ss,sg,ZcorX2)

        # compute Variance
        VRstn,VRp[zz] = fitcheck(ss,sg,W,N_STN,isoflag,ZcorX2)

        # Set quality
        Qual = SetQuality(VRp[zz])

        # Orgenize file to plot figure
        dataTp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, MTp, Zcor_Fp, VR_Fp = MakeMTpyFiles(ZcorX2, Sampl4inv, zz,sg,ss,MM, N_STN, data_Tp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, MTp, Zcor_Fp, VR_Fp, VRstn)

        print('Depth %d | Quality = %d | VR = %5.2f | Mw %3.2f' % (DEPTH[zz], Qual, VRp[zz], MWzp[zz]))

    QmlEvent = MK_MT_event(MTp, Origine_time, Lat0, Long0, EventN, MWzp, VRp, DEPTH, Mo, Green, pdcp, pclvdp, pisop, FMSp, STN_LISTN)
    Event1 = Event[0][0].replace(':', '_')

    try:
        QmlEvent.write(pathplt+'/'+Event1+'.xml', format="SC3ML")
    except:
        QmlEvent.write(pathplt+'/'+Event1+'.xml', format="QUAKEML")


    # Plot final figure
    Path2pdf = PlotMomentTensor(2000, DEPTH, Depth0, VRp, Sampl4inv, FMSp, MWzp, data_Tp, data_Rp, data_Zp, synth_Tp, synth_Rp, synth_Zp, STN_LISTN, DIST, BAZ, Zcor_Fp, Sampl4inv, VR_Fp, smpl, MTp,DIST, Lat0, Long0, LON_STN, LAT_STN, faults1, Green, Event, EventN, pathplt, GFtype, Mo, pdcp, pclvdp, pisop, MODE, AUTOMODE,fcutlow, fcuthigh)

    return Path2pdf



def GETEVENT4TDMT(Event, dTime_search):

    #-------------READ ANTELOPE CATALOG -----------------------------

    file_LOC = '/Users/nadavwetzler/Dropbox/Public/DataSet/AllIsrael/All_Israel_1900-082018/Gitterman2005/Reloc_GII_1900-082018_Gitt05_Origins.txt'
    LatC = []
    LonC = []
    DepC = []
    MagC = []
    DTVC = []
    f_LOC=open(file_LOC)
    lines=f_LOC.readlines()
    for jj in range(1, len(lines)):
        line = lines[jj].split()
        LatC.append(line[2])
        LonC.append(line[3])
        DepC.append(line[4])
        MagC.append(line[10])
        YMD = line[5].split('/')
        HMS = line[7]
        date = line[5].split()
        datetime = UTCDateTime(YMD[2]+'-'+YMD[0]+'-'+YMD[1]+'T'+HMS)
        DTVC.append(datetime)

    #CLIENT = "http://82.102.143.46:8181"
    #client = Client(CLIENT)

    #------------------ Get event info and stations------------------

    Origine_time = UTCDateTime(Event[0])

    dTt = np.asarray(DTVC) - Origine_time
    pos = np.argmin(abs(dTt))
    mint = min(abs(dTt))
    if mint > 3:
        print('Could not fine event in catalog')
        cat = client.get_events(starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, includearrivals=True)
        posCat = 0
        if len(cat) > 1:
            print("None single event is retrieved! reducing catalog to the nearest...")
            t3 = []
            for jj in range(len(cat)):
                t3.append(cat.events[jj].origins[0].time)
            t3v = np.asarray(t3)
            posCat = np.argmin(abs(t3v - Origine_time))
            cat = cat[posCat]
            # sys.exit("None single event is retrieved! use smaller dTime_search")
        elif len(cat) == 0:
            sys.exit("Event is not found! use larger dTime_search or specify location ")
        Lat0 = cat.events[0].origins[posCat].latitude
        Long0 = cat.events[0].origins[posCat].longitude
        Time0 = cat.events[0].origins[posCat].time
        Depth0 = cat.events[0].origins[posCat].depth/1000.0
        M0 = cat.events[0].magnitudes[0].mag
        print("Event found in relocataed catalog, using relocation coordinates")
    else:
        Lat0 = float(LatC[pos])
        Long0=float(LonC[pos])
        Depth0=round(float(DepC[pos]))
        M0 = float(MagC[pos])

    if Depth0 < 1:
        Depth0 = 1

    return Lat0, Long0, Depth0, M0


def MK_SEIS_FILES(Lat0,Long0, Depth0, Mw0,STN_LIST0S, Event, fdsn, MaxDist2stn, MinDist2stn, CHNP, invlength0, GF0, Fmin0, Fmax0, Taper0):
    if len(CHNP) < 4:
        print('length of ChanelPriority needs to be 4 like: "HBSH"')

    clients, out, faultsT, fdsn_l = GetClient(Lat0,Long0,fdsn)

    Dist_max0,Dist_min0 = MW2DIST(Mw0)

    if MaxDist2stn > 0:
        Dist_max = MaxDist2stn
    else:
        Dist_max = Dist_max0

    if MinDist2stn > 0:
        Dist_min = MinDist2stn
    else:
        Dist_min = Dist_min0


    Event = Event[0]
    if len(Event) == 1:
        Origine_time = UTCDateTime(Event[0])
    else:
        Origine_time = UTCDateTime(Event)

    Green, Min_Depth, Max_Depth, fcutlow, fcuthigh, dt0 = whichGreen(Mw0, out, GF0, Fmin0, Fmax0)
    GFtype0 =  Green.split('_')
    GFtype = GFtype0[1]
    Sampl4inv, dPt1 = INV_LENGTH_T(Mw0, invlength0)

    STN_LISTS        = np.zeros(0)
    STN_CHN_LISTS    = np.zeros(0)
    DistS            = np.zeros(0)
    bazS             = np.zeros(0)
    S2NS             = np.zeros(0)
    Longs_stnS       = np.zeros(0)
    Lats_stnS        = np.zeros(0)
    North_dirS       = np.zeros(0)
    N_STNS           = 0
    SEIS0S           = Stream()
    SEIS_DISP_ROTS   = Stream()
    SEIS_DISP_ROT_FS = Stream()

    for si in range(fdsn_l):
        STN_LIST0 = STN_LIST0S
        if fdsn_l > 1:
            client = clients[si]
        else:
            client = clients[0]
        INT_R_check = 0

        if STN_LIST0 == 'ALLA':
            STN_LIST0 = []
            try:
                INT_A = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="21", channel="??Z", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
                try:
                    INT_B = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="?HZ", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
                    INT_R = INT_B + INT_A
                    INT_R_check = 1
                except:
                    INT_R = INT_A
                    INT_R_check = 1
            except:
                INT_B = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="B?Z", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
                INT_R = INT_B
                INT_R_check = 1

        elif STN_LIST0 == 'ALLB':

            STN_LIST0 = []

            # INT_H = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="H?Z", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
            # INT_B = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="B?Z", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
            # INT_S = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="S?Z", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
            INT_R = client.get_stations(starttime=Origine_time - dTime_search, endtime=Origine_time + dTime_search, network=NTW, sta="*", loc="*", channel="HHZ,BHZ,SHZ,EHZ", latitude=Lat0, longitude=Long0, minradius=Dist_min/100, maxradius=Dist_max/100, level="response")
            # INT_R = INT_H + INT_B + INT_S
            INT_R_check = 1
            # print(INT_R)
        else:
            stntmp = STN_LIST0
            STN_LIST0 = []
            STN_LIST0.append(stntmp)

        if INT_R_check ==1:
            for kk in range(len(INT_R)):
                for ii in range(len(INT_R[kk])):
                    stn_name = INT_R.networks[kk].stations[ii].code
                    STN_LIST0.append(stn_name)
        # Get stn info
        STN_LIST0 = list(dict.fromkeys(STN_LIST0))
        STN_LIST0 = np.array(STN_LIST0)
        STN_LIST0 = np.unique(STN_LIST0)
        N_STN = len(STN_LIST0)
        Dist      = np.zeros(N_STN)
        az        = np.zeros(N_STN)
        baz       = np.zeros(N_STN)
        Lats_stn  = np.zeros(N_STN)
        Longs_stn = np.zeros(N_STN)
        S2N       = np.zeros(N_STN)
        S2NF      = np.zeros(N_STN)
        isDataV   = np.zeros(N_STN)
        North_dir = np.zeros(N_STN)
        STN_LIST_F = []
        STN_LIST = []
        STN_CHN_LIST = []
        SEIS0 = Stream()
        SEIS_DISP_ROT = Stream()
        SEIS_DISP_ROT_F = Stream()

        if out == 0:
            try:
                cat = client.get_events(starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, includearrivals=True)
            except:
                cat = 0
                okcat = 0
        print(" Running over %s stations" % N_STN)
        for ii in range(N_STN):
            # print(" Read %s" %STN_LIST0[ii])
            ok_stn = 1
            if  STN_LIST0[ii] not in STN_LISTS:
                try:
                    INT = client.get_stations(network=NTW, channel="%sH?" % CHNP[0], station=STN_LIST0[ii], starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, level="response")
                except:
                    try:
                        INT = client.get_stations(network=NTW, channel="%sH?" % CHNP[1], station=STN_LIST0[ii], starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, level="response")
                    except:
                        try:
                            INT = client.get_stations(network=NTW, channel="%sH?" % CHNP[2], station=STN_LIST0[ii], starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, level="response")
                        except:
                            try:
                                INT = client.get_stations(network=NTW, channel="%sN?" % CHNP[3], station=STN_LIST0[ii], starttime=Origine_time-dTime_search, endtime=Origine_time+dTime_search, level="response")
                            except:
                                ok_stn = 0
                                print('Station %s not found in inventory' % STN_LIST0[ii])
            else:
                ok_stn = 0 
            # if STN_LIST0[ii] == 'MMA0B':
            #     print('Excluding MMA0B')
            #     ok_stn = 0
            if ok_stn == 1:
                if len(INT.networks[0].stations[0].channels) >= 3:
                    Lats_stn[ii]  = INT.networks[0].stations[0].latitude
                    Longs_stn[ii] = INT.networks[0].stations[0].longitude
                    [Dist[ii], az[ii], baz[ii]] = gps2dist_azimuth(Lats_stn[ii], Longs_stn[ii], Lat0, Long0, a=6378137.0, f=0.0033528106647474805)
                    Dist[ii] = round(Dist[ii]/1000.0, 0)  # convert m to km
                    az[ii] = round(az[ii], 0)
                    baz[ii] = round(baz[ii], 0)

                    chan_north = '??N'
                    code = INT.networks[0].stations[0].channels[0].code

                    for jj in range(len(INT.networks[0].stations[0].channels)):
                        if INT.networks[0].stations[0].channels[jj].code[2] == '1':
                            chan_north = '??1'
                            code = INT.networks[0].stations[0].channels[jj].code

                    invN = INT
                    # invN = invN.select(component='N')
                    invN = invN.select(channel=chan_north)
                    try:
                        stn_north = invN.networks[0].stations[0].channels[0].azimuth
                    except:
                        stn_north = 0

                    # print('%s %d ' % (STN_LIST0[ii], stn_north))
                    North_dir[ii] = stn_north

                    # print(STN_LIST0[ii])
                    # print('Station %s not found in inventory' % STN_LIST0[ii])

                    if 'cat' in globals():
                        Pt = GetPstn(cat,STN_LIST0[ii])
                    else:
                        arrivals = model.get_ray_paths(source_depth_in_km=Depth0, distance_in_degree=Dist[ii]*0.009,phase_list=["P","p","Pn"])
                        len_arrivals = len(arrivals)
                        if len_arrivals > 0:
                            t_arrivals = np.zeros(len_arrivals)
                            for jj in range(len_arrivals):
                                t_arrivals[jj] = arrivals[jj].time
                            Pt = Origine_time + np.min(t_arrivals)
                        else:
                            Pt = Origine_time

                    SEIST, isData, chn = LOAD_STN_SEIS(STN_LIST0[ii], code[0:2], Pt, Sampl4inv, client, dPt1, Taper0, North_dir[ii],INT)

                    if isData == 1:
                        S2N[ii] = SIGNAL2NOISE(dPt1, SEIST, 30, 0)
                        SEIS0 += SEIST
                        try:
                            S1, S2 = DISP_ROT_FILT(SEIST, az[ii], fcutlow, fcuthigh, n_corners, dt0, dPt1)
                            SEIS_DISP_ROT += S1
                            SEIS_DISP_ROT_F += S2
                            S2NF[ii] = SIGNAL2NOISE(dPt1, SEIS_DISP_ROT_F, 30, 0)
                            isDataV[ii] = 1
                            STN_CHN_LIST.append(chn)
                            STN_LIST.append(STN_LIST0[ii])
                            print('Station %s %s ok!' % (STN_LIST0[ii], chn))
                        except:
                            print('Unable to remove response for %s' % STN_LIST0[ii])
                    else:
                        print('No data for %s' % STN_LIST0[ii])

        SEIS0.merge()
        SEIS_DISP_ROT.merge()
        SEIS_DISP_ROT_F.merge()

        #SEIS0.select(component='Z')
        #SEIS0.plot(equal_scale=False)

        I         = isDataV > 0
        Dist      = Dist[I]
        baz       = baz[I]
        az        = az[I]
        Longs_stn = Longs_stn[I]
        Lats_stn  = Lats_stn[I]
        S2N       = S2N[I]
        North_dir = North_dir[I]



        STN_CHN_LIST = np.asarray(STN_CHN_LIST)
        STN_LIST     = np.asarray(STN_LIST)
        N_STN        = len(STN_LIST)
        N_STNS       = N_STNS + N_STN

        STN_CHN_LISTS = np.concatenate((STN_CHN_LISTS,STN_CHN_LIST))
        STN_LISTS     = np.concatenate((STN_LISTS,STN_LIST))
        DistS         = np.concatenate((DistS,Dist))
        bazS          = np.concatenate((bazS,baz))
        S2NS          = np.concatenate((S2NS,S2N))
        Longs_stnS    = np.concatenate((Longs_stnS,Longs_stn))
        Lats_stnS     = np.concatenate((Lats_stnS,Lats_stn))
        North_dirS    = np.concatenate((North_dirS,North_dir))

        SEIS0S+=SEIS0
        SEIS_DISP_ROTS+=SEIS_DISP_ROT
        SEIS_DISP_ROT_FS+=SEIS_DISP_ROT_F

    STN_LISTS = np.asarray(STN_LISTS)
    STN_CHN_LISTS = np.asarray(STN_CHN_LISTS)

    I_dist = np.argsort(DistS)
    DistS = DistS[I_dist]
    bazS = bazS[I_dist]
    Longs_stnS = Longs_stnS[I_dist]
    Lats_stnS = Lats_stnS[I_dist]
    S2NS = S2NS[I_dist]
    STN_CHN_LISTS = STN_CHN_LISTS[I_dist]
    STN_LISTS = STN_LISTS[I_dist]

    return STN_LISTS, STN_CHN_LISTS, DistS, bazS, S2NS, Longs_stnS, Lats_stnS, N_STNS, SEIS0S, SEIS_DISP_ROTS, SEIS_DISP_ROT_FS, North_dirS, Green, Min_Depth, Max_Depth, Sampl4inv, dt0, faultsT, GFtype, fcutlow, fcuthigh, client, Origine_time, dPt1
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################




def run_TDMTW(Origintime, AUTOMODE, STN_LIST0, Lat0, Long0, Depth0, Mw0, fdsn,dr0, MaxDist2stn, MinDist2stn, CHNP, invlength0, GF0, Fmin0, Fmax0, Taper0, Taper2):
    Event = [Origintime]

    STN_LIST, STN_CHN_LIST, Dist, baz, S2N, Longs_stn, Lats_stn, N_STN, SEIS0, SEIS_DISP_ROT, SEIS_DISP_ROT_F, North_dir, Green, Min_Depth, Max_Depth, Sampl4inv, smpl, faultsT, GFtype, fcutlow, fcuthigh, client, Origine_time, dPt1 = MK_SEIS_FILES(Lat0,Long0, Depth0, Mw0,STN_LIST0, Event, fdsn,MaxDist2stn, MinDist2stn, CHNP, invlength0,GF0,Fmin0, Fmax0, Taper0)

    if N_STN > 0:
        if AUTOMODE == 0:

            dLatLon = (np.max(Dist) / 100) * 1.3
            MinMaxLon = np.zeros(2)
            MinMaxLat = np.zeros(2)
            MinMaxLon[0:2] = [Long0-dLatLon, Long0+dLatLon]
            MinMaxLat[0:2] = [Lat0-dLatLon, Lat0+dLatLon]
            FIG_DPI = 100
            FIGSIZE = 6
            COLOR_OFF_STN = [0.6,0.6,0.6]
            COLOR_ON_STN = 'g'
            LARGE_FONT = 12

            N_Plot_max = 10

            VLENGTH = 30 * N_Plot_max
            root = Tk()
            root.title("GSI TDMTW")
            if VLENGTH < 900:
                root.geometry("1900x900")
                VLENGTH = 900
            else:
                root.geometry("1900x%s" % VLENGTH)
            nrows = np.max([N_Plot_max + 2, 10])
            ncols = 25

            for ii in range(ncols):
                root.rowconfigure(ii, weight=0)
            for ii in range(nrows):
                root.columnconfigure(ii,weight=0)

            class Application(Frame):


                def __init__(self, master):
                    Frame.__init__(self,master)
                    # Frame.__init__(self,master, width = 2000, height = 2000)
                    self.grid()
                    self.create_widgets()


                def create_widgets(self):
                    """Create buttons """
                    #Label(self, text = "EVENT: %s" % Event[0]).grid(row = 0, column = 0, sticky = W)

                    Color0 = 'lightblue'

                    dxPos = 0.01
                    Xplace = [0.25,0.4,0.5]
                    Yplace = [0.1,0]
                    Yplace[1] =  1 - (2*dxPos + Yplace[0])
                    Xplace[2] =  1-4*dxPos - (Xplace[0] + Xplace[1])


                    self.frame0 = Frame(root,bg=Color0)
                    self.frame0.place(relx=dxPos, rely=dxPos, relwidth=1.0-2*dxPos,relheight=Yplace[0])
                    self.frame1 = Frame(root)  #, bg='#80c1ff')
                    self.frame1.place(relx=dxPos, rely=Yplace[0]+dxPos, relwidth=Xplace[0],relheight=Yplace[1])
                    self.frame2 = Frame(root)  #, bg='#60c1ff')
                    self.frame2.place(relx=dxPos + Xplace[0], rely=Yplace[0]+dxPos, relwidth=Xplace[1], relheight=Yplace[1])
                    self.frame3 = Frame(root)  #, bg='#20c1ff')
                    self.frame3.place(relx=dxPos + Xplace[0]+ Xplace[1], rely=Yplace[0]+dxPos, relwidth=Xplace[2], relheight=0.9)

                    # Isotropic option
                    tkvar = StringVar(root)

                    # Dictionary with options
                    choices = { 'No','Yes'}
                    tkvar.set('No') # set the default option
                    self.vol_p = 0
                    self.iconPath = LOGO
                    IMG = Image.open(self.iconPath)
                    logo_size = int(VLENGTH * Yplace[0] *1.0)
                    IMG = IMG.resize((logo_size,logo_size), Image.ANTIALIAS)
                    photoImg =  ImageTk.PhotoImage(IMG)
                    self.icon = photoImg
                    self.icon_size = Label(self.frame0)
                    # self.icon_size.image = self.icon  # <== this is were we anchor the img object
                    self.icon_size.create_image = self.icon
                    self.icon_size.configure(image=self.icon)
                    self.icon_size.grid(row=0, column=12, columnspan=3, rowspan=3)

                    #Button(self.frame0, image=logo).grid(row=0, column=12)

                    OrigintimeTXT = Event[0]
                    OrigintimeTXT = OrigintimeTXT[0]
                    OrigintimeTXT = OrigintimeTXT.replace('T',' ')

                    Label(self.frame0, text = "%s" % OrigintimeTXT, bg=Color0).grid(row = 1, column = 0)
                    Label(self.frame0, text = "Origin Time", bg=Color0).grid(row = 0, column = 0)
                    popupMenu = OptionMenu(self.frame0, tkvar, *choices)
                    Label(self.frame0, text="ISO", bg=Color0).grid(row = 0, column = 4)
                    popupMenu.grid(row = 1, column = 4)
                    Label(self.frame0, text = "Latitude", bg=Color0).grid(row = 0, column = 5, sticky = W)
                    Label(self.frame0, text = "%s" % Lat0, bg=Color0).grid(row = 1, column = 5, sticky = W)
                    #Entry(self, textvariable = .grid(row=1, column = 9, sticky = W)
                    Label(self.frame0, text = "Longitude", bg=Color0).grid(row = 0, column = 6, sticky = W)
                    Label(self.frame0, text = "%s" % Long0, bg=Color0).grid(row = 1, column = 6, sticky = W)
                    Label(self.frame0, text = "Magnitude", bg=Color0).grid(row = 0, column = 7, sticky = W)
                    Label(self.frame0, text = "%s" % Mw0, bg=Color0).grid(row = 1, column = 7)

                    Button(self.frame0, text='RUN!!!', bg=Color0, command=self.plotfinalstn).grid(row=1, column=9)
                    Button(self.frame0, text='SEISMIC', bg=Color0, command=self.plotseis_stn).grid(row=1, column=8,sticky=W)


                    def textupdates(root, txt):
                        T = Text(root, height=2, width=30).grid(row=0, column=9)
                        T.insert(END, txt)

                    def change_dropdown(*args):
                        print( tkvar.get() )
                        if tkvar.get() == 'No':
                            self.vol_p = 0
                        else:
                            self.vol_p = 1


                    # link function to change dropdown
                    tkvar.trace('w', change_dropdown)


                    # Plot Regional seismicity
                    okcat, cat = Get_Reg_EQ(client, dDays, Origine_time)

                    # Add scrollbar to frame1 (stations vertical list)
                    # Add a canvas in that frame.
                    canvas1 = Canvas(self.frame1, bg="White")
                    canvas1.grid(row=0, column=0)

                    # Create a vertical scrollbar linked to the canvas.
                    vsbar = Scrollbar(self.frame1, orient=VERTICAL, command=canvas1.yview)
                    vsbar.grid(row=0, column=1, sticky=NS)
                    canvas1.configure(yscrollcommand=vsbar.set)

                    # Create a horizontal scrollbar linked to the canvas.
                    hsbar = Scrollbar(self.frame1, orient=HORIZONTAL, command=canvas1.xview)
                    hsbar.grid(row=1, column=0, sticky=EW)
                    canvas1.configure(xscrollcommand=hsbar.set)

                    # Create a frame on the canvas to contain the buttons.
                    self.buttons_frame = Frame(canvas1, bd=2)


                    # Instructions
                    self.button = []
                    self.buttonSeis = []
                    self.Zcor_e = []
                    for ii in range(N_STN):
                        stn_button = BooleanVar()
                        Checkbutton(self.buttons_frame, text='%s  (%s)  %d km S2N: %3.2f' % (STN_LIST[ii], STN_CHN_LIST[ii], Dist[ii], S2N[ii]), variable = stn_button, command = self.update_stn).grid(row = 2 + ii, column = 0, columnspan=1, rowspan=1, sticky=W)
                        self.button.append(stn_button)
                        seis_button = Button(self.buttons_frame, text="Seismic", command=lambda c=ii: self.plotseis(c)).grid(row = 2 + ii, column = 1)
                        seis_R_button = Button(self.buttons_frame, text="DRF", command=lambda c1=ii: self.plotseisR(c1)).grid(row = 2 + ii, column = 2)

                        zcor_txt = 0.0
                        self.Zcor_e.append(Entry(self.buttons_frame, width = 3))
                        self.Zcor_e[ii].grid(row = 2 + ii, column = 3)
                        self.Zcor_e[ii].focus_set()
                        # Zcor_button = Button(self.frame1, text="get", width=3, command=lambda c1=ii: self.zcorseis(c1)).grid(row = 2 + ii, column = 4)


                    # Create canvas window to hold the buttons_frame.
                    canvas1.create_window((0,0), window=self.buttons_frame, anchor=NW)

                    self.buttons_frame.update_idletasks()  # Needed to make bbox info available.
                    bbox = canvas1.bbox(ALL)  # Get bounding box of canvas with Buttons.
                    #print('canvas.bbox(tk.ALL): {}'.format(bbox))

                    # Define the scrollable region as entire canvas with only the desired
                    # number of rows and columns displayed.
                    ROWS, COLS = N_STN, 6  # Size of grid.
                    ROWS_DISP = 20  # Number of rows to display.
                    COLS_DISP = 6  # Number of columns to display.

                    w, h = bbox[2]-bbox[1], bbox[3]-bbox[1]
                    dw, dh = int((w/COLS) * COLS_DISP), int((h/ROWS) * ROWS_DISP)
                    canvas1.configure(scrollregion=bbox, width=dw, height=dh)

                    f = Figure(figsize=(FIGSIZE, 2*FIGSIZE), dpi=FIG_DPI)
                    ax = f.add_subplot(211)
                    ax.plot(Long0,Lat0,marker='o', color='r',markeredgecolor='g')

                    for ii in range(N_STN):
                        ax.plot(Longs_stn[ii],Lats_stn[ii],marker='o', color=COLOR_OFF_STN,label=STN_LIST[ii])
                        ax.text(Longs_stn[ii],Lats_stn[ii],STN_LIST[ii], color='k')
                    ax.set_xlim(MinMaxLon)
                    ax.set_ylim(MinMaxLat)



                    if okcat == 1:
                        eqR, eqLat, eqLon, eqM, eqT = Calc_eq_reg_R(cat, Lat0, Long0, Origine_time)
                        R0 = Mw2a(Mw0)
                        eqR = eqR / 1000
                        R0r = 3*R0
                        Iin = eqR < R0r
                        eqT = eqT / (24*3600)
                        eqTin = eqT[Iin]
                        eqMin = eqM[Iin]
                        if len(eqMin) > 0:
                            x_c, y_c = Make_circle_R(R0r, Lat0, Long0)
                            Rc = 20/Mw0
                            ax2 = f.add_subplot(212,aspect = Rc)
                            ax2.scatter(eqTin, eqMin,5,color='r')
                            for ii in range(len(eqTin)):
                                ax2.plot([eqTin[ii], eqTin[ii]], [0, eqMin[ii]],color='k', linewidth=0.5)
                            # ax2.stem(eqT[Iin], eqM[Iin], linewidth=0.5)
                            ax2.set_xlim([-dDays, dDays])
                            ax2.set_ylim([0, np.max(eqMin) + 0.2])
                            ax2.set_xlabel('Days from event, R=%d km' % np.round(R0))
                            ax2.set_ylabel('Magnitude')
                            ax.scatter(eqLon,eqLat,eqM,color='k')
                            ax.scatter(eqLon[Iin],eqLat[Iin],eqM[Iin],color='r')
                            ax.plot(x_c, y_c)

                    canvas = FigureCanvasTkAgg(f, self.frame2)

                    canvas.get_tk_widget().grid(row=0, column=0,columnspan=FIGSIZE, rowspan= FIGSIZE)
                    self.frame2.canvas = canvas
                    self.frame2.canvas.draw()
                    ax.osm = OSM(ax)
                    self.frame2.canvas.figure.axes[0].osm.draw()
                    self.frame2.canvas.draw()

                def plotseis_stn(self):
                    ''' Plots all vertical channels'''
                    for widget in self.frame3.winfo_children():
                        widget.destroy()
                    SEISZ = SEIS0.select(component="Z").copy()
                    lw = 0.3
                    times1 = np.zeros(N_STN)
                    times2 = np.zeros(N_STN)
                    AMPS = np.zeros(N_STN)

                    seis1 = Stream()
                    for jj in range(N_STN):
                        tr = SEISZ.select(station=STN_LIST[jj]).copy()
                        tr[0].stats.distance = Dist[jj] * 1000
                        seis1.append(tr[0])
                    t1 = seis1[0].stats.starttime
                    t2 = seis1[0].stats.endtime
                    seis1.trim(t1+50,t2-150)
                    seis1.filter('highpass',freq=1,zerophase=True)
                    for tr in seis1.select(station='K10B'):
                        seis1.remove(tr)
                    # seis1.plot(equal_scale=False, type='section')

                    seis1.write('/Users/nadavwetzler/Desktop/eventZ.mseed',format="MSEED")

                    orderSTN = np.argsort(Dist)
                    for jj in range(N_STN):
                        tr = SEISZ.select(station=STN_LIST[jj]).copy()
                        times1[jj] = tr[0].stats.starttime
                        times2[jj] = tr[0].stats.endtime
                        AMPS[jj] = np.max([abs(np.min(tr[0].data)), abs(np.max(tr[0].data))])
                    t1 = np.min(times1)
                    t2 = np.max(times2)
                    FIGSIZEy = N_STN
                    if FIGSIZEy > FIGSIZE*2:
                        FIGSIZEy = FIGSIZE*2
                    f4 = Figure(figsize=(FIGSIZE, FIGSIZEy), dpi=FIG_DPI)
                    f4.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)
                    for jj in range(N_STN):
                        ii = orderSTN[jj]
                        Ymax = AMPS[ii]
                        tr = SEISZ.select(station=STN_LIST[ii]).copy()
                        dt = tr[0].stats.delta
                        ln = len(tr[0].data)
                        tend = dt*ln
                        dtx = times1[ii] - t1
                        tx = np.arange(0,tend,dt) + dtx
                        ax1 = f4.add_subplot(N_STN,1,jj+1)
                        ax1.plot(tx, tr[0].data,linewidth = lw)
                        ax1.text(dPt1*2,Ymax/2,'%s N' % STN_LIST[ii])
                        ax1.set_xlim([dPt1-15, t2-t1])
                        ax1.set_ylim([-Ymax,Ymax])
                        canvas3 = FigureCanvasTkAgg(f4, self.frame3)
                        canvas3.get_tk_widget().grid(row=0, column=0)
                        canvas3.draw()
                        canvas3.show()


                def plotseis(self, ii):
                    for widget in self.frame3.winfo_children():
                        widget.destroy()
                    trN = SEIS0.select(component="N", station=STN_LIST[ii]).copy()
                    trE = SEIS0.select(component="E", station=STN_LIST[ii]).copy()
                    trZ = SEIS0.select(component="Z", station=STN_LIST[ii]).copy()
                    Ymax = np.max(np.abs([np.min(trN[0].data), np.min(trE[0].data), np.min(trZ[0].data),np.max(trN[0].data), np.max(trE[0].data), np.max(trZ[0].data)]))
                    lw = 0.2
                    ln = len(trN[0].data)
                    dt = trN[0].stats.delta
                    tend = dt*ln
                    tx = np.arange(0,tend,dt)
                    f3 = Figure(figsize=(FIGSIZE, FIGSIZE), dpi=FIG_DPI)
                    ax1 = f3.add_subplot(311)
                    ax1.plot(tx, trN[0].data,linewidth = lw)
                    ax1.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax1.text(dPt1/2,Ymax/2,'N')
                    ax1.set_xlim([0, tend])
                    ax1.set_ylim([-Ymax,Ymax])
                    ax1.yaxis.tick_right()
                    ax1.set_title('%s' % STN_LIST[ii])
                    ax2 = f3.add_subplot(312)
                    ax2.plot(tx, trE[0],linewidth = lw)
                    ax2.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax2.text(dPt1/2,Ymax/2,'E')
                    ax2.set_xlim([0, tend])
                    ax2.set_ylim([-Ymax,Ymax])
                    ax2.yaxis.tick_right()
                    ax3 = f3.add_subplot(313)
                    ax3.plot(tx, trZ[0].data,linewidth = lw)
                    ax3.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax3.text(dPt1/2,Ymax/2,'Z')
                    ax3.set_xlabel('SEC.')
                    ax3.set_xlim([0, tend])
                    ax3.set_ylim([-Ymax,Ymax])
                    ax3.yaxis.tick_right()
                    f3.subplots_adjust(left=0.05, right=0.8, top=0.9, bottom=0.1)
                    self.frame3.config(bg="white")
                    canvas3 = FigureCanvasTkAgg(f3, self.frame3)
                    canvas3.get_tk_widget().grid(row=0, column=0)
                    canvas3.draw()
                    canvas3.show()

                def plotseisR(self, ii):
                    for widget in self.frame3.winfo_children():
                        widget.destroy()
                    trN = SEIS_DISP_ROT.select(component="R", station=STN_LIST[ii]).copy()
                    trNf = SEIS_DISP_ROT_F.select(component="R", station=STN_LIST[ii]).copy()
                    trE = SEIS_DISP_ROT.select(component="T", station=STN_LIST[ii]).copy()
                    trEf = SEIS_DISP_ROT_F.select(component="T", station=STN_LIST[ii]).copy()
                    trZ = SEIS_DISP_ROT.select(component="Z", station=STN_LIST[ii]).copy()
                    trZf = SEIS_DISP_ROT_F.select(component="Z", station=STN_LIST[ii]).copy()
                    Ymax = np.max(np.abs([np.min(trN[0].data), np.min(trE[0].data), np.min(trZ[0].data),np.max(trN[0].data), np.max(trE[0].data), np.max(trZ[0].data)]))
                    Ymax2 = np.max(np.abs([np.min(trNf[0].data), np.min(trEf[0].data), np.min(trZf[0].data),np.max(trNf[0].data), np.max(trEf[0].data), np.max(trZf[0].data)]))
                    lw = 0.2
                    lwf = 1.2
                    ln = len(trN[0].data)
                    lnf = len(trNf[0].data)
                    dt = trN[0].stats.delta
                    dtf = trNf[0].stats.delta
                    tx = np.arange(0,ln)*dt
                    tend = np.max(tx)
                    txf = np.arange(0,lnf)*dtf
                    f3 = Figure(figsize=(FIGSIZE, FIGSIZE), dpi=FIG_DPI)
                    ax1 = f3.add_subplot(311)
                    ax1.plot(tx, trN[0].data,linewidth = lw)
                    ax1.plot(txf, trNf[0].data/Ymax2*Ymax,linewidth = lwf)
                    ax1.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax1.text(dPt1/2,Ymax/2,'R')
                    ax1.set_xlim([0, tend])
                    ax1.set_ylim([-Ymax,Ymax])
                    ax1.yaxis.tick_right()
                    ax1.set_title('BP Filter: %s - %s %s' % (fcutlow, fcuthigh, STN_LIST[ii]))
                    ax2 = f3.add_subplot(312)
                    ax2.plot(tx, trE[0].data,linewidth = lw)
                    ax2.plot(txf, trEf[0].data/Ymax2*Ymax,linewidth = lwf)
                    ax2.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax2.text(dPt1/2,Ymax/2,'T')
                    ax2.set_xlim([0, tend])
                    ax2.set_ylim([-Ymax,Ymax])
                    ax2.yaxis.tick_right()
                    ax3 = f3.add_subplot(313)
                    ax3.plot(tx, trZ[0].data,linewidth = lw)
                    ax3.plot(txf, trZf[0].data/Ymax2*Ymax,linewidth = lwf)
                    ax3.plot([dPt1,dPt1],[-Ymax,Ymax])
                    ax3.text(dPt1/2,Ymax/2,'Z')
                    ax3.set_xlabel('SEC.')
                    ax3.set_xlim([0, tend])
                    ax3.set_ylim([-Ymax,Ymax])
                    ax3.yaxis.tick_right()
                    f3.subplots_adjust(left=0.05, right=0.8, top=0.9, bottom=0.1)
                    canvas3 = FigureCanvasTkAgg(f3, self.frame3)
                    canvas3.get_tk_widget().grid(row=0, column=0)
                    canvas3.draw()
                    canvas3.show()
                def zcorseis(self, ii):
                    print(self.Zcor_e.get())


                def update_stn(self):

                    okstn = np.zeros(N_STN) - 1
                    ZcorVal = np.zeros(N_STN)
                    for ii in range(N_STN):
                        if self.button[ii].get():
                            self.frame2.canvas.figure.axes[0].lines[ii+1].set_color(COLOR_ON_STN)
                            okstn[ii] = ii

                            if self.Zcor_e[ii].get():
                                ZcorVal[ii] = self.Zcor_e[ii].get()
                                print(ZcorVal)

                        else:
                            self.frame2.canvas.figure.axes[0].lines[ii+1].set_color(COLOR_OFF_STN)


                    Iok = okstn > -1
                    okstn1 = okstn[Iok]
                    self.Final_STN_list = []
                    self.Final_Dist = []
                    self.Final_baz = []
                    self.Lon_Stn = []
                    self.Lat_Stn = []
                    self.Zecor = []
                    for ii in range(len(okstn1)):
                        jj = int(okstn1[ii])
                        self.Final_STN_list.append(STN_LIST[jj])
                        self.Final_baz.append(baz[jj])
                        self.Final_Dist.append(Dist[jj])
                        self.Lon_Stn.append(Longs_stn[jj])
                        self.Lat_Stn.append(Lats_stn[jj])
                        self.Zecor.append(ZcorVal[jj])

                    self.frame2.canvas.figure.axes[0].osm.draw()
                    self.frame2.canvas.draw()
                    for ii in range(len(self.Zecor)):
                        print('%s %s' % (self.Final_STN_list[ii], self.Zecor[ii]))


                def plotfinalstn(self):
                    print(self.Final_STN_list)
                    ZCORS = np.zeros(len(self.Final_STN_list))
                    D1, D2 = SetDepthRang(Depth0, Min_Depth, Max_Depth, dr0)

                    TDMTRUN(Depth0, Lat0, Long0, GFtype, Green, self.Final_STN_list, self.Final_Dist, self.Final_baz, self.Lat_Stn, self.Lon_Stn, SEIS_DISP_ROT_F, Sampl4inv, smpl, D1, D2, self.vol_p, self.Zecor, faultsT,AUTOMODE, Event, SEIS0, fcutlow, fcuthigh, Taper2)


            app = Application(root)
            root.mainloop()

        else:
            Istn = S2N > thresh_s2n
            Final_STN_list = STN_LIST[Istn]
            if len(Final_STN_list) == 0:
                print('S2N too low, no stations are selected')
                Path2pdf = '-1'
                return Path2pdf
            else:
                Final_Dist = Dist[Istn]
                Final_baz = baz[Istn]
                Final_Lat_Stn = Lats_stn[Istn]
                Final_Lon_Stn = Longs_stn[Istn]
                vol_p = 0
                ZCORS = np.zeros(len(Final_STN_list))
                D1, D2 = SetDepthRang(Depth0,Min_Depth, Max_Depth,dr0)

                print('Start inversion!')
                Path2pdf = TDMTRUN(Depth0, Lat0, Long0, GFtype, Green, Final_STN_list, Final_Dist, Final_baz, Final_Lat_Stn, Final_Lon_Stn, SEIS_DISP_ROT_F, Sampl4inv, dt0, D1, D2, vol_p, ZCORS, faultsT, AUTOMODE, Event, SEIS0, fcutlow, fcuthigh,Taper2)
                return Path2pdf
    else:
        print('No stations are detected')
        Path2pdf = '-1'
        return Path2pdf
if __name__ == '__main__':


    args = parser.parse_args()
    print(args.Origintime)
    # log = logging.getLogger('TDMTW')
    # log.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d | %(name)s | %(levelname)s | %(message)s', datefmt='%Y-%m-%dT%H:%M:%S')
    # filehandler = TimedRotatingFileHandler('TDMTW.log',
    #                                    when='midnight',
    #                                    utc=True)
    # filehandler.setFormatter(formatter)
    # filehandler.setLevel(logging.DEBUG)
    # log.addHandler(filehandler)




    if not args.Origintime:
        print("Testing EVENT !!!!")


        args.Origintime = "2018-10-11T16:55:07" # Island
        args.Lat0 =  30.9572
        args.Long0 = 35.4711
        args.Depth0 = 5
        args.Mw0 = 2.4
        args.Invlength = 20
        args.Client = 'GSI'
        args.STN_LIST0 = 'ALLA'
        args.Fmin = 0.5
        args.Fmax = 1.5
        args.Taper2 = 1
        args.Radious = 35
        args.GreenFunction = 'GREEN_DSB__0.1_0.0-0.0_fullNNN_1.0_99.0/MSEED/'

        # args.Origintime = "2021-07-27T22:12:15" # Island
        # args.Lat0 =  64.53
        # args.Long0 = -17.53
        # args.Depth0 = 10
        # args.Mw0 = 5.0
        # args.Client = 'IRIS'
        # args.STN_LIST0 = 'BORG'
        # args.Fmin = 0.02
        # args.Fmax = 0.05

        # args.Origintime = "2021-07-28T08:27:07" # DSB
        # args.Lat0 =  31.156
        # args.Long0 = 35.292
        # args.Depth0 = 5
        # args.Mw0 = 2.9
        # args.Radious = 35
        # args.Client = 'GSI'


        # args.Origintime = "2019-05-15T16:53:49" # E.Med
        # args.Lat0 =  32.77
        # args.Long0 = 32.824
        # args.Depth0 = 14
        # args.Mw0 = 4.5
        # args.Radious = 250
        # args.Client = 'GSI,GFZ'
        # # args.STN_LIST0 = 'IZRL'
        # args.Fmin = 0.02
        # args.Fmax = 0.05

        # args.Origintime = " 2021-07-04T13:59:13" # south Galilee
        # args.Lat0 =  32.444
        # args.Long0 = 35.127
        # args.Depth0 = 13
        # args.Mw0 = 2.7
        # args.Radious = 50
        # args.Client = 'GSI'
        # args.STN_LIST0 = 'IZRL'
        # args.Fmin = 0.2
        # args.Fmax = 1.5
        # args.Invlength = 50

        # args.Origintime = "2021-06-18T19:48:39" # Explosion Florida
        # args.Lat0 =  29.742
        # args.Long0 = -79.347
        # args.Depth0 = 1
        # args.Mw0 = 4.0
        # args.Client = 'IRIS'
        # args.STN_LIST0 = 'ALLA'
        # args.Fmin = 0.015
        # args.Fmax = 0.2
        # args.Invlength = 180

        # args.Origintime = "2021-06-17T05:37:58" # test by Andrey
        # args.Lat0 =  31.92
        # args.Long0 = 35.53
        # args.Depth0 = 13
        # args.Mw0 = 6.5
        # args.Radious = 350
        # args.Client = 'SYNTH'
        # args.STN_LIST0 = 'ALLA'
        # args.Fmin = 0.015
        # args.Fmax = 0.2
        # args.Invlength = 350


        # args.Origintime = "2021-06-15T10:35:36" # test by Andrey
        # args.Lat0 =  31.94
        # args.Long0 = 35.51
        # args.Depth0 = 14
        # args.Mw0 = 6.5
        # args.Radious = 350
        # args.Client = 'SYNTH'
        # args.STN_LIST0 = 'ALLA'
        # args.Fmin = 0.015
        # args.Fmax = 0.09
        # args.Invlength = 350


        # args.Origintime = "2021-06-08T10:22:59" # test
        # args.Lat0 =  31.934 # 31.934 # SC3 31.9215 # 20/70/0 M0=7.19*10^18
        # args.Long0 = 35.5136 # 35.5136 # SC3 35.6015 #
        # args.Depth0 = 14
        # args.Mw0 = 6.5
        # args.Radious = 350
        # args.Client = 'SYNTH'
        # args.STN_LIST0 = 'LOD'
        # args.Fmin = 0.02
        # args.Fmax = 1.5
        # args.Invlength = 80

        # args.Origintime = "2021-06-08T13:10:41" # test
        # args.Lat0 =  31.8714 # 31.8714 # SC3 31.8717
        # args.Long0 = 35.5136 # 35.5136 # SC3 35.5008
        # args.Depth0 = 14
        # args.Mw0 = 5.3 # 5.03
        # args.Radious = 450
        # args.Client = 'SYNTH'
        # args.STN_LIST0 = 'ALLA'
        # args.Fmin = 0.02
        # args.Fmax = 1.5
        # args.Invlength = 30


        # args.Origintime = "2020-01-14T00:25:24" # DSB
        # args.Lat0 =  31.4864
        # args.Long0 = 35.5224
        # args.Depth0 = 14
        # args.Mw0 = 2.4
        # args.Fmin = 0.2
        # args.Fmax = 1.5
        # args.Radious = 50
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'MSLM'
        # args.Invlength = 20
        # args.GreenFunction = 'GREEN_DSB__0.1_0.0-0.0_fullNNN_1.0_99.0/MSEED/'
        # args.Taper2 = 1

        # args.Origintime = "2021-05-31T06:23:41" # Cyprus
        # args.Lat0 =  35.075
        # args.Long0 = 32.769
        # args.Depth0 = 31
        # args.Mw0 = 4.0
        # args.Fmin = 0.02
        # args.Fmax = 0.05
        # args.Radious = 350
        # args.Client = 'IRIS'
        # args.STN_LIST0 = 'MSLM'
        # args.Invlength = 20
        # args.GreenFunction = 'GREEN_DSB__0.1_0.0-0.0_fullNNN_1.0_99.0/MSEED/'
        # args.Taper2 = 1

        #
        # args.Origintime = "2019-09-04T03:45:00" # DSB
        # args.Lat0 =  31.16
        # args.Long0 = 35.493
        # args.Depth0 = 21
        # args.Mw0 = 2.7
        # args.Fmin = 0.2
        # args.Fmax = 0.5
        # args.Radious = 60
        # args.Client = 'ISN'
        # # args.STN_LIST0 = 'GHAJ'
        # args.Invlength = 35
        # # args.GreenFunction = 'GREEN_DSB__0.1_0.0-0.0_fullNNN_1.0_99.0/MSEED/'
        # args.Taper2 = 1

        # args.Origintime = "2019-10-14T15:55:43" # DSB
        # args.Lat0 =  31.5098
        # args.Long0 = 35.5137
        # args.Depth0 = 13
        # args.Mw0 = 3.0
        # args.Fmin = 0.2
        # args.Fmax = 0.5
        # # args.Radious = 150
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'GHAJ'
        # args.GreenFunction = 'GREEN_DSB__0.1_0.0-0.0_fullNNN_1.0_99.0/MSEED/'
        # args.Taper2 = 1

        # args.Origintime = "2021-02-16T21:57:00" # Cyprus
        # args.Lat0 =  34.539
        # args.Long0 = 33.940
        # args.Depth0 = 30
        # args.Mw0 = 3.2
        # args.Fmin = 0.02
        # args.Fmax = 0.05
        # args.Radious = 150
        # args.Client = 'IRIS'
        # args.STN_LIST0 = 'CY602'

        # args.Origintime = "2021-03-03T10:16:11" # Grees
        # args.Lat0 =  39.77
        # args.Long0 = 22.14
        # args.Depth0 = 10
        # args.Mw0 = 6.3
        # # args.Fmin = 0.02
        # # args.Fmax = 0.05
        # # args.Radious = 50
        # args.Client = 'NOA'
        # args.STN_LIST0 = 'APE'

        # args.Origintime = "2019-04-29T12:07:25" # LRB
        # args.Lat0 =  34.192
        # args.Long0 = 36.152
        # args.Depth0 = 1
        # args.Mw0 = 3.8 # 3.1
        # args.Fmin = 0.05
        # args.Fmax = 0.1
        # args.Radious = 500
        # args.Client = 'ISN,KOERI'
        # args.Invlength = 180

        # args.Origintime = "2020-12-27T06:37:34" # Turkey
        # args.Lat0 =  38.51
        # args.Long0 = 39.32
        # args.Depth0 = 10
        # args.Mw0 = 5.4
        # args.Fmin = 0.02
        # args.Fmax = 0.05
        # args.Radious = 350
        # args.Client = 'KOERI'
        # args.Invlength = 180

        # args.Origintime = "2020-08-04T15:08:18" # Beirut blast
        # args.Lat0 =  33.9
        # args.Long0 = 35.52
        # args.Depth0 = 0
        # args.Mw0 = 5.5 # 3.1
        # # args.STN_LIST0 = 'GHAJ'
        # args.Client = 'ISN'
        # args.dr0 = 2

        # args.Origintime = "2020-06-16T14:30:27"
        # args.Lat0 =  27.4
        # args.Long0 = 34.69
        # args.Depth0 = 5
        # args.Mw0 = 5.2
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'EIL'

        # args.Origintime = "2020-11-30T16:36:29"
        # args.Lat0 =  8.1103
        # args.Long0 = -82.91
        # args.Depth0 = 2
        # args.Mw0 = 4.8
        # args.Client = 'OVSICORI'
        # args.STN_LIST0 = 'PEZE'
        # args.Fmin = 0.02
        # args.Fmax = 0.05
        # args.Radious = 200

        # args.Origintime = "2020-09-15T09:24:02"
        # args.Lat0 =  32.505
        # args.Long0 = 35.199
        # args.Depth0 = 2
        # args.Mw0 = 5.2
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'EIL'

        # args.Origintime = "2020-09-01T07:52:03"
        # args.Lat0 =  31.19
        # args.Long0 = 35.32
        # args.Depth0 = 1
        # args.Mw0 = 2.5
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'ALLA'

        # args.Origintime = "2020-08-11T08:40:25"
        # args.Lat0 =  32.98
        # args.Long0 = 35.49
        # args.Depth0 = 10
        # args.Mw0 = 3.7 # 3.1
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'MRON'


        # args.Origintime = "2020-03-22T05:24:02"
        # args.Lat0 =  45.88
        # args.Long0 = 15.99
        # args.Depth0 = 5
        # args.Mw0 = 5.3 # 3.1

        # args.Origintime = "2018-07-04T19:51:24"
        # args.Lat0 =  32.8406
        # args.Long0 = 35.5809
        # args.Depth0 = 5
        # args.Mw0 = 2.6 # 3.5
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'ALLA'



        # args.Origintime = "2019-11-25T11:39:31"
        # args.Lat0 =  31.3502
        # args.Long0 = 35.5328
        # args.Depth0 = 12
        # args.Mw0 = 3.5 # 3.5
        # args.Client = 'ISN'
        # args.Radious = 70
        # args.ChannelPriority = 'EEEE'
        # args.Invlength = 30

        # args.Origintime = "2019-09-25T23:46:46.1"
        # args.Lat0 =  -3.56
        # args.Long0 = 128.45
        # args.Depth0 = 10
        # args.Mw0 = 6.3
        # args.Client = 'GFZ'
        # args.AUTOMODE = 1


        # args.Origintime = "2019-11-01T05:25:46"
        # args.Lat0 =  40.5
        # args.Long0 = 20.79
        # args.Depth0 = 10
        # args.Mw0 = 5.4
        # # args.Client = 'ISN'
        # args.AUTOMODE = 1

        # args.Origintime = "2019-11-01T18:40:41.32"
        # args.Lat0 = 32.8636
        # args.Long0 = 35.5706
        # args.Depth0 = 8
        # args.Mw0 = 3.2
        # args.Client = 'ISN'
        # args.AUTOMODE = 0

        # args.Origintime = "2018-07-04T01:50:07.39"
        # args.Lat0 = 32.8471
        # args.Long0 = 35.5905
        # args.Depth0 = 7
        # args.Mw0 = 4.2
        # args.Client = 'ISN'
        # args.STN_LIST0 = 'ALLB'
        # args.Radious = 20
        # args.Invlength = 30

        # args.Origintime = "2019-10-14T15:55:43"
        # args.Lat0 = 31.510
        # args.Long0 = 35.514
        # args.Depth0 = 13
        # args.Mw0 = 3.0
        # args.Client = 'ISN'

        # print("Testing EVENT Local!!!!")
        # args.Origintime = "2019-09-16T10:35:31"
        # args.Lat0 = 29.975
        # args.Long0 = 34.403
        # args.Depth0 = 7
        # args.Mw0 = 3.6
        # args.STN_LIST0 = 'EIL'

        # args.Origintime = "2019-08-08T11:25:31.0"
        # args.Lat0 =  37.94
        # args.Long0 = 29.59
        # args.Depth0 = 10
        # args.Mw0 = 5.7
        # args.Client = 'GFZ'
        # args.AUTOMODE = 0
        # args.STN_LIST0 = 'APE'


        # args.Origintime = "1995-11-22T04:15:11"
        # args.Lat0 =  28.5
        # args.Long0 = 34.6
        # args.Depth0 = 15
        # args.Mw0 = 7.2
        # args.Client = 'GFZ'
        # args.AUTOMODE = 0
        # args.GreenFunction = 'GREEN_Git05__1.0_0.02-0.05_fullNN/'

        # args.Origintime = "2019-05-15T16:53:47.4"
        # args.Lat0 =  32.8107
        # args.Long0 = 32.7967
        # args.Depth0 = 10
        # args.Mw0 = 5.7
        # args.Client = 'GFZ'
        # args.AUTOMODE = 0


        # args.Origintime = "2019-09-30T11:53:21"
        # args.Lat0 = -22.23
        # args.Long0 = -68.66
        # args.Depth0 = 110.0
        # args.Mw0 = 5.6
        # args.Client = 'GFZ'
        # args.AUTOMODE = 1


        # args.Origintime = "2019-11-27T07:23:40"
        # args.Lat0 = 35.7
        # args.Long0 = 23.19
        # args.Depth0 = 54.0
        # args.Mw0 = 6.0
        # args.Client = 'GFZ'
        # args.AUTOMODE = 0

    if not args.Client:
        print('Automated selection of FDSN server')
        args.Client = 'AUTO'
    else:
        print('Selected FDSN: %s' % args.Client)


    run_TDMTW([args.Origintime], args.AUTOMODE, args.STN_LIST0, args.Lat0, args.Long0, args.Depth0, args.Mw0, args.Client, args.DepthRang, args.Radious, args.Radious0, args.ChannelPriority, args.Invlength, args.GreenFunction, args.Fmin, args.Fmax, args.Taper, args.Taper2)


