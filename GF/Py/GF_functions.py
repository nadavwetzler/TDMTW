import os
import subprocess, sys
import platform
from obspy import Stream, Trace, read
import numpy as np
import datetime
import pandas as pd
now = datetime.datetime.now()
tt=now.timetuple()
b2s = "b2s.par"

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


def set_path_gf(name, path2folder,  dt0, f_low, f_high, DIST1, DIST2):
    Path2_Green = '%s/GREEN_%s__%s_%s-%s_fullNNN_%s_%s/' % (path2folder, name, dt0, f_low, f_high, DIST1, DIST2)

    if os.path.exists(Path2_Green):
        os.chdir(Path2_Green)

    else:
        os.mkdir(Path2_Green)
        os.chdir(Path2_Green)
    return Path2_Green


def loadVmodel(name, path2model):
    model = pd.read_csv('%s/%s.csv' % (path2model,name))
    MODEL = np.array(model)
    MODEL[0] = MODEL[0] + 0.01
    return MODEL
    
# ----------build b2s.par file --------------------------

def mk_b2s(Path2_Green, npts0, dt0):
    b2s_F_name = Path2_Green +'/'+ b2s

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
    filt_P_name = Path2_Green +'/'+ f_name
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


def mk_gf1(MODEL, MODELN, DEPTH, npts0,dt0, Path2_bin, Path2_Green,MT_F, DIST, ReductionVel):
    dH = 0.05 # Thickness of artificial layer in case Z == layer interface

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
    MT_F_name = Path2_Green +'/'+ MT_F +"_" + str(DEPTH) + ".model.txt"
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
    os.system("%s/FKRPROG<%s" % (Path2_bin, MT_F_name))

    # run wvint9 d
    os.system("%s/wvint9d" % Path2_bin)


    # ---------------Make perl fk file -------------------
    fk_name = 'run_fkrsortiso'
    if os.path.exists(fk_name): os.remove(fk_name)
    fk_P_name = Path2_Green +'/'+ fk_name
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
    fk_P.write('mkHelm format="(6e13.5)" ntr=10 dt=$dt nt=$npts < junk > %s_{$dist[$j]}d{$depth}.disp\n' % MT_F)
    fk_P.write('./run_filtsyniso %s_{$dist[$j]}d{$depth}.disp %s_{$dist[$j]}d{$depth}\n' % (MT_F, MT_F))
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

