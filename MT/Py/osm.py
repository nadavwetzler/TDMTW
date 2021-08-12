#!/usr/bin/env python
# /********************************************************************************
#    Copyright (C) by Ran Novitsky Nof                                            *
#                                                                                 *
#    This file is part of E2ReviewTool                                            *
#                                                                                 *
#    E2ReviewTool is free software: you can redistribute it and/or modify         *
#    it under the terms of the GNU Lesser General Public License as published by  *
#    the Free Software Foundation, either version 3 of the License, or            *
#    (at your option) any later version.                                          *
#                                                                                 *
#    This program is distributed in the hope that it will be useful,              *
#    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
#    GNU Lesser General Public License for more details.                          *
#                                                                                 *
#    You should have received a copy of the GNU Lesser General Public License     *
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
# *********************************************************************************
import sys
#sys.path.insert(0, '/Users/nadavwetzler/Dropbox/Moment-tensor/Py_MT/TDMTW/')
import numpy as np
import scipy.misc
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO
import os
from httplib2 import Http
from PIL import Image
import matplotlib
import logging

# logging
logging.basicConfig(
    format='%(asctime)s | %(levelname)s | %(message)s', datefmt='%Y-%m-%dT%H:%M:%S')
log = logging.getLogger('OSM')
#log.setLevel(logging.DEBUG)
log.setLevel(logging.ERROR)

# util class for Open Street Map

TILESARCHIVE = os.path.split(__file__)[0] + os.sep + 'tiles'

class OSM(object):
    latmax = 90.0
    latmin = -90.0
    lonmax = 180.0
    lonmin = -180.0
    maxlevel = 17
    maxtilecash = 500
    # Server = "http://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/"
    Server = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/"
    def __init__(self, ax, tileurl=str(Server), tilepat="{Z}/{Y}/{X}.png", tilearchive=TILESARCHIVE, mplconnect=True):
        # constants and defaults
        self.cashedtiles = {}
        self.tilesorder = []
        self.currentimages = []
        self.latmax = 90.0
        self.latmin = -90.0
        self.lonmax = 180.0
        self.lonmin = -180.0
        self.maxlevel = 17
        self.maxtilecash = 500
        # constants and defaults
        self.ax = ax  # where to plot tiles
        self.ax.set_aspect('equal') # make sure it is equale aspect
        self.canvas = self.ax.figure.canvas
        if mplconnect:
            self.mplconnect()
        self._button = False
        self.http = Http()
        self.tileurl = tileurl
        self.tilepat = tilepat
        if not os.path.exists(tilearchive):
            log.warning('No tile archive was set.')
            tilearchive = ''  # make sure tile archive directory exists
        self.tilearchive = tilearchive
        # Add copyright text
        self.ax.text(0.99, 0.01, 'Leaflet | Tiles ' + u'\N{COPYRIGHT SIGN}' + 'Esri', color='#555555', ha='right', transform=ax.transAxes, name='Helvetica', fontsize=8, zorder=999)
   
    def mplconnect(self):
        self.canvas.mpl_connect('draw_event', self.draw)
        self.canvas.mpl_connect('button_press_event', self.click)
        self.canvas.mpl_connect('button_release_event', self.unclick)        
        
    def click(self,event):
        self._button = self.canvas._button

    def unclick(self,event):
        self._button = None
        self.draw()

    def tile2image(self, datafile, limits=[latmin, latmax, lonmin, lonmax]):
        'add tile to axes images'
        y0, y1, x0, x1 = limits  # get tile extent
        try:
            # data = scipy.misc.fromimage(Image.open(datafile))  # read tile image and convert to 2D array
            data = np.asarray(Image.open(datafile))
        except IOError:
            log.error("Bad image file: {}. Please remove the file for next time.".format(datafile))
            data = np.zeros((256, 256, 3))  # use black image with red X.
            data[np.where(np.identity(256, dtype=np.uint8))] = [1, 0, 0]
            data[np.where(np.identity(256, dtype=np.uint8)[::-1])] = [1, 0, 0]
        im = matplotlib.image.AxesImage(self.ax)  # create an image object
        im.set_data(data)  # set the data to image
        im._extent = [x0, x1, y0, y1]  # set image extents
        return im

    def relim(self, x0, x1, y0, y1):
        'set axes limits'
        self.ax.set_ylim(y0, y1)
        self.ax.set_xlim(x0, x1)
        self.ax.apply_aspect()

    def relimcorrected(self, x0, x1, y0, y1):
        'fix axes limits to geo limits'
        x0, x1 = self.checklonrange(x0, x1)
        y0, y1 = self.checklatrange(y0, y1)
        self.relim(x0, x1, y0, y1)

    def checklonrange(self, lonmin, lonmax):
        'fix longitude limits'
        # The bounds are chosen such that they give the correct results up
        # to zoom level 30 (zoom levels up to 18 actually make sense):
        # lon2tilex(-180.0, 30) == 0
        # lon2tilex(179.9999999, 30) == 1073741823 == 2**30 - 1
        if lonmin < -180.0:
            lonmin = -180.0
        if lonmin > 179.9999999:
            lonmin = 179.9999999
        if lonmax < -180.0:
            lonmax = -180.0
        if lonmax > 179.9999999:
            lonmax = 179.9999999
        return lonmin, lonmax

    def checklatrange(self, latmin, latmax):
        'fix latitude limits'
        # The bounds are chosen such that they give the correct results up
        # to zoom level 30 (zoom levels up to 18 actually make sense):
        # lat2tiley(85.0511287798, 30) == 0
        # lat2tiley(-85.0511287798, 30) == 1073741823 == 2**30 - 1
        if latmin < -85.0511287798:
            latmin = -85.0511287798
        if latmin > 85.0511287798:
            latmin = 85.0511287798
        if latmax < -85.0511287798:
            latmax = -85.0511287798
        if latmax > 85.0511287798:
            latmax = 85.0511287798
        return latmin, latmax

    def lon2tilex(self, lon, zoom):
        'calculate tile longitude index at zoom level'
        return int(((lon + 180) / 360.0 * 2**zoom))

    def tilex2lims(self, tilex, zoom):
        'calculate tile longitude limits at zoom level'
        return tilex * 360.0 / 2.0**zoom - 180.0, (tilex + 1) * 360.0 / 2.0**zoom - 180.0

    def lat2tiley(self, lat, zoom):
        'calculate tile latitude index at zoom level'
        lata = lat * np.pi / 180.0
        return int(((1 - np.log(np.tan(lata) + 1.0 / np.cos(lata)) / np.pi) / 2.0 * 2**zoom))

    def tiley2lims(self, tiley, zoom):
        'calculate tile latitude limits at zoom level'
        return np.degrees(np.arctan(np.sinh(np.pi * (1 - 2 * (tiley + 1) / 2.0**zoom)))), np.degrees(np.arctan(np.sinh(np.pi * (1 - 2 * tiley / 2.0**zoom))))

    def tiles2lims(self, tilex, tiley, zoom):
        'calculate tile limits at zoom level'
        return np.append(self.tiley2lims(tiley, zoom), self.tilex2lims(tilex, zoom))

    def tile2path(self, tilex, tiley, zoom):
        'get tile path'
        if not os.path.exists(os.sep.join([self.tilearchive, self.tilepat.format(Z=zoom, X=tilex, Y=tiley)])):
            return self.tileurl + self.tilepat.format(Z=zoom, X=tilex, Y=tiley)  # download if not archived
        else:
            return os.sep.join([self.tilearchive, self.tilepat.format(Z=zoom, X=tilex, Y=tiley)])  # get form archive

    def tilepath2zoomxy(self, tile):
        'convert tile path to lat-lon indices and zoom'
        tileorder = os.path.splitext(self.tilepat.replace('{', '').replace('}', ''))[0].split('/')[-3:]
        tilexyz = dict(zip(tileorder, os.path.splitext(tile)[0].split(os.sep)[-3:]))
        return int(tilexyz['X']), int(tilexyz['Y']), int(tilexyz['Z'])

    def gettiles(self, lonmin, lonmax, latmin, latmax):
        'get a list of tiles for limits. make sure at least two tiles cover the range'
        latmin, latmax = self.checklatrange(latmin, latmax)  # correct lat limits to tiles
        lonmin, lonmax = self.checklonrange(lonmin, lonmax)  # correct lon limits to tiles
        w, h = self.ax.figure.canvas.get_width_height()  # get figure dimensions
        for i in range(self.maxlevel):  # go through zoom levels
            xtiles = list(set([self.lon2tilex(lon, i) for lon in [lonmin, lonmax]]))  # get lon indices of tiles
            ytiles = list(set([self.lat2tiley(lat, i) for lat in [latmin, latmax]]))  # get lat indices of tiles
            if len(xtiles) == 2 and len(ytiles) == 2 and (abs(np.diff(xtiles)[0]) > round(w / 256) or
                                                          abs(np.diff(ytiles)[0]) > (h / 256)):  # make sure resolution and details are acceptable
                break  # we have found the desired zoom level
        # get tiles indices limits
        mintilex = self.lon2tilex(lonmin, i)
        maxtilex = self.lon2tilex(lonmax, i)
        mintiley = self.lat2tiley(latmax, i)
        maxtiley = self.lat2tiley(latmin, i)
        tiles = []
        # get the list of tiles at indices limits
        for x in range(mintilex, maxtilex + 1):
            for y in range(mintiley, maxtiley + 1):
                tiles.append(self.tile2path(x, y, i))
        return tiles

    def plottiles(self, tiles):
        # remove old tiles
        for im in self.currentimages:
            try:
                self.ax.images.remove(im)
            except Exception as ex:
                log.warning(ex)
        self.currentimages = []
        # read or download tiles and plot them to axes
        for t in tiles:
            tID = t.replace(self.tileurl, '').replace(self.tilearchive, '').replace('/', os.sep)
            if tID not in self.cashedtiles:  # if tiles are not already in cash
                if t.startswith('http'):  # if we need to download
                    log.debug('Downloading {}'.format(t))
                    try:
                        datafile = StringIO(self.http.request(t)[1])  # download tile image to stream
                    except Exception as ex: # in case of disconnection
                        log.warning(ex)
                        self.http = Http()
                        datafile = StringIO(self.http.request(t)[1])
                    if self.tilearchive:  # archive tile if we set the archive
                        d = datafile.read()  # read the downloaded data from stream
                        datafile.seek(0)  # rewind stream
                        archivedir, fname = os.path.split(os.sep.join([self.tilearchive, tID]))  # get archive path of tile
                        if not os.path.exists(archivedir):  # make sure archive directory exists
                            os.makedirs(archivedir)  # or create it
                        with open(os.sep.join([archivedir, fname]), 'wb') as f:
                            f.write(d)  # save tile to archive
                else:
                    log.debug('Reading tile from {}'.format(t))
                    datafile = open(t, 'rb')  # or get the tile file name
                limits = self.tiles2lims(*self.tilepath2zoomxy(t))
                im = self.tile2image(datafile, limits)
                datafile.close()
                self.cashedtiles[tID] = im
                self.tilesorder.append(tID)
                if len(self.tilesorder) > self.maxtilecash:
                    old = self.tilesorder.pop(0)
                    self.cashedtiles.pop(old)
            else:
                im = self.cashedtiles[tID]
            self.ax._set_artist_props(im)
            self.ax.images.append(im)  # add image to axes
            self.currentimages.append(im)  # add image to current images list
        self.ax.stale = True

    def draw(self,*args,**kwargs):
        if self._button:
                return # don't redraw while panning/zooming
        log.info('Updating map, please hold...')
        self.ax.apply_aspect()  # make sure xlim and ylim are updated to screen size. this is because we use equal aspect and datalim. see matplotlib details on axes set_aspect function
        x0, x1 = self.ax.get_xlim()  # get requested limits of x axis
        y0, y1 = self.ax.get_ylim()  # get requested limits of y axis
        self.relimcorrected(x0, x1, y0, y1)  # make sure limits are not out of map phisical boundaries. see osm module for more details.
        tiles = self.gettiles(x0, x1, y0, y1)  # get the needed tiles from buffer or url. see osm module for more details
        self.plottiles(tiles)  # plot the tiles. see osm module for more details
        log.info('Updating finished')
