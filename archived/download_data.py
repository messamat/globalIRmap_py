"""Developer: Mathis L. Messager
Purpose: Download data to update and enhance RiverATLAS hydro-environmental attributes, including:
   - SoilGrids250m v2
   - MODIS 250 m land and water mask
   - Global Land Analysis & Discovery (GLAD) global inland water dynamics data
   - ALOS digital elevation model (DEM) data
"""

#pip install gsutil --ignore-installed six

import os
import arcpy
from arcpy.sa import *
import sys
import re
import subprocess
from utility_functions import *
import json
from cookielib import CookieJar
from urllib import urlencode
import urllib2
from bs4 import BeautifulSoup
import urlparse
import math
import numpy as np
from collections import defaultdict
import tarfile
from functools import reduce

#Folder structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
datdir = os.path.join(rootdir, '/HydroATLAS', 'data')
resdir = os.path.join(rootdir, '/HydroATLAS' , 'results')
glad_dir = os.path.join(datdir, 'GLAD')
pathcheckcreate(glad_dir)

# ------------------------------- Download SoilGrids250 m data -------------------------------------------------------------
#How to make the soil mask? https://github.com/ISRICWorldSoil/SoilGrids250m/blob/1a04d4214c0efabfe15abd621a6c299c173e402c/grids/GlobCover30/soilmask_250m.R

#Create output directory
sg_outdir = os.path.join(datdir, 'SOILGRIDS250')
pathcheckcreate(sg_outdir)

# Download metadata
dlfile("https://files.isric.org/soilgrids/data/recent/META_GEOTIFF_1B.csv", sg_outdir)

#Parse HTML to get all available layers
sg_https = "https://files.isric.org/soilgrids/data/recent/"
sg_r = urllib2.urlopen(sg_https)
sg_soup = BeautifulSoup(sg_r, features="html.parser")

sg_lyrdict = defaultdict(list)
for link in sg_soup.findAll('a', attrs={'href': re.compile(".*[.]tif$")}):
    sg_lyrdict[re.search("^[a-zA-Z1-9]+(?=_)", link.get('href')).group()].append(
        urlparse.urljoin(sg_https, link.get('href')))

#Download all layers of interest
#for lyrk in sg_lyrdict.keys():
lyrk = "SLGWRB"
if lyrk not in ["TAXNWRB", "TAXOUSDA"]:
    outdir = os.path.join(sg_outdir, lyrk)
    pathcheckcreate(outdir)

    for lyrurl in sg_lyrdict[lyrk]:
        print(lyrurl)
        #dlfile(url=lyrurl, outpath=outdir, outfile=os.path.split(lyrurl)[1], fieldnames=None)

#------------------------------- Download MODIS 250 m land and water mask  -------------------------
# Create output directory
mod44w_outdir = os.path.join(datdir, 'mod44w')
pathcheckcreate(mod44w_outdir)

# Download metadata
dlfile(url='https://lpdaac.usgs.gov/documents/109/MOD44W_User_Guide_ATBD_V6.pdf', outpath=mod44w_outdir)

# Parse HTML to get all available layers
mod44w_https = "https://e4ftl01.cr.usgs.gov/MOLT/MOD44W.006/2015.01.01/"
mod44w_r = urllib2.urlopen(mod44w_https)
mod44w_soup = BeautifulSoup(mod44w_r, features="html.parser")

# The user credentials that will be used to authenticate access to the data
with open("../configs.json") as json_data_file:  # https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)

# Download all layers of interest
for lyrurl in [urlparse.urljoin(mod44w_https, link.get('href')) for link in
               mod44w_soup.findAll('a', attrs={'href': re.compile("MOD44W.*[.]hdf([.]xml)*")})]:

        dlfile(url=lyrurl, outpath=mod44w_outdir, outfile=os.path.split(lyrurl)[1], fieldnames=None,
               loginprompter="https://urs.earthdata.nasa.gov",
               username=authdat['earthdata']['username'], password=authdat['earthdata']['password'])

# ----------------------------- Download GLAD data ----------------------------------------------------------------------
glad_dtype = "class99_19.tif"
gsutil_ls_cmd = "gsutil ls gs://earthenginepartners-hansen/water/*/{}".format(glad_dtype)
glad_cloudout = subprocess.check_output(gsutil_ls_cmd)
glad_cloudlist = glad_cloudout.split('\n')

for tile in glad_cloudlist:
    tileroot = os.path.split(os.path.split(tile)[0])[1]
    out_tilen = os.path.join(glad_dir, '{0}_{1}.tif'.format(
        os.path.splitext(glad_dtype)[0], tileroot))
    if not os.path.isfile(out_tilen):
        print(tile)
        gsutil_cp_cmd = "gsutil cp {0} {1}".format(tile, glad_dir)
        subprocess.check_output(gsutil_cp_cmd)
        os.rename(os.path.join(glad_dir, glad_dtype), out_tilen)

# ----------------------------- Download ALOS DEM data ----------------------------------------------------------------------
#Main page https://www.eorc.jaxa.jp/ALOS/en/aw3d30/
alos_outdir = os.path.join(datdir, 'ALOS')
pathcheckcreate(alos_outdir)

# Import user credentials that will be used to authenticate access to the data
with open("../configs.json") as json_data_file:  # https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)

#Source Javascript to get data on tile click: https://www.eorc.jaxa.jp/ALOS/en/aw3d30/data/html_v2003/js/dsm_dl_select_v2003.js?ver=20200330
alosurl = "https://www.eorc.jaxa.jp/ALOS/en/aw3d30/data/html_v2003/dl/download_v2003.htm?N075E010_N079E011"
def parse_alosurl(alosurl):
    urlsub = alosurl[alosurl.index('?') + 1:]
    downurl = reduce(urlparse.urljoin,
                     ['https://www.eorc.jaxa.jp/ALOS/aw3d30/data/release_v2003/',
                      '{}/'.format(urlsub[0:8]),
                      '{}.zip'.format(urlsub[9:9 + 8])])

#Generate server-based list of urls
ytileset1 = {'N{}'.format(str(x).zfill(3)) for x in xrange(0, 95, 5)} | \
            {'S{}'.format(str(x).zfill(3)) for x in xrange(0, 60, 5)}
xtileset1 = {'W{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)} | \
            {'E{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)}

for x1 in xtileset1:
    for y1 in ytileset1:
        if x1[0] == 'W':
            xtileset2 = {'{0}{1}'.format(x1[0], str(i).zfill(3)) for i in xrange(int(x1[1:]) - 5, int(x1[1:]))}
        else:
            xtileset2 = {'{0}{1}'.format(x1[0], str(i).zfill(3)) for i in xrange(int(x1[1:]), int(x1[1:])+5)}

        if y1[0] == 'S':
            ytileset2 = {'{0}{1}'.format(y1[0], str(i).zfill(3)) for i in xrange(int(y1[1:])-5, int(y1[1:]))}
        else:
            ytileset2 = {'{0}{1}'.format(y1[0], str(i).zfill(3)) for i in xrange(int(y1[1:]), int(y1[1:]) + 5)}

        for x2 in xtileset2:
            for y2 in ytileset2:
                alos_tileurl = reduce(urlparse.urljoin,
                                      ['https://www.eorc.jaxa.jp/ALOS/aw3d30/data/release_v2003/',
                                       '{0}{1}/'.format(y1, x1),
                                       '{0}{1}.zip'.format(y2, x2)])
                alos_outtile = os.path.join(alos_outdir, '{0}{1}.zip'.format(y2, x2))
                print(alos_tileurl)
                try:
                    if not os.path.exists(alos_outtile):
                        dlfile(url=alos_tileurl, outpath=alos_outdir, outfile=os.path.split(alos_outtile)[1],
                               fieldnames=None,
                               ignore_downloadable=False,
                               loginprompter="https://www.eorc.jaxa.jp",
                               username=authdat['alos']['username'], password=authdat['alos']['password'])
                    else:
                        print('{} was already processed...'.format(alos_outtile))
                except:
                    traceback.print_exc()



########################################################################################################################
# ######################################################################################################################
# #############################        Not used in the end      ########################################################
"""
#pip install pyproj==1.9.6 owslib==0.18 - 0.19 dropped python 2.7
from owslib.wcs import WebCoverageService  # OWSlib module to access WMS services from SDAT

#------------------------------- Download HYSOGS250 m data -------------------------------------------------------------
hysogdir = os.path.join(datdir, 'hysog')
pathcheckcreate(hysogdir)
hysogresgdb = os.path.join(resdir, 'hysog.gdb')
pathcheckcreate(hysogresgdb)
hysogmosaic = os.path.join(hysogresgdb, 'hysog_mosaic')

#DOI of dataset is https://doi.org/10.3334/ORNLDAAC/1566
#This last number is the dataset ID to be used for search in WebMapService
sdatwcs = WebCoverageService('https://webmap.ornl.gov/ogcbroker/wcs')
print(str(len(sdatwcs.contents)) + ' layers found from ' + sdatwcs.identification.title)
# filter layers
hysog_wcsid = filter(lambda x: x.startswith('1566_'), sdatwcs.contents)[0]
print(hysog_wcsid)

hysogbblist = divbb(bbox=sdatwcs[hysog_wcsid].boundingBoxWGS84,
                    res=sdatwcs[hysog_wcsid].grid.offsetvectors[0][0],
                    divratio=10)
hysogoutlist = ['{0}_{1}.tif'.format(os.path.join(hysogdir, 'hysog'), i)
                for i in xrange(0, len(hysogbblist))]
if not all([os.path.isfile(i) for i in hysogoutlist]):
    x=0
    for bb in hysogbblist:
        #print(bb)
        outtile = hysogoutlist[x]
        if not os.path.isfile(outtile):
            print(outtile)
            hysog_wc = sdatwcs.getCoverage(identifier=hysog_wcsid,
                                           bbox=bb,
                                           crs='EPSG:4326',
                                           format='Geotiff_BYTE',
                                           interpolation='NEAREST',
                                           resx=sdatwcs[hysog_wcsid].grid.offsetvectors[0][0],
                                           resy=sdatwcs[hysog_wcsid].grid.offsetvectors[0][0])

            with open(outtile, "wb") as local_file:
                local_file.write(hysog_wc.read())
        else:
            print("{} already exists...".format(outtile))

        x+=1

# #Only keep tiles with data - only works with numpy > 1.9.3 but breaks arcpy
# for tilepath in hysogoutlist:
#     print(tilepath)
#     tiledat = gdal_array.LoadFile(tilepath)
#     if tiledat.max() == 0:
#         hysogoutlist.remove(tilepath)

#mosaic them
print('Mosaicking hysogs tiles...')
arcpy.MosaicToNewRaster_management(hysogoutlist, output_location=hysogresgdb,
                                   raster_dataset_name_with_extension= 'hysog_mosaic',
                                   pixel_type= '8_BIT_UNSIGNED',
                                   number_of_bands = 1,
                                   mosaic_colormap_mode = 'FIRST')


# ------------------------------- Download EarthEnv-DEM90 data -------------------------------------------------------------
ee_outdir = os.path.join(datdir, 'earthenv')
pathcheckcreate(ee_outdir)

# Parse HTML to get all available layers
ee_https = "https://www.earthenv.org/DEM.html" #Doesn't work, return 403 with urlibb2

# xrange(0, 85, 5) #xrange(0, 60, 5)
ytileset = {'N{}'.format(str(x).zfill(2)) for x in xrange(0, 90, 5)} | \
           {'S{}'.format(str(x).zfill(2)) for x in xrange(0, 60, 5)}
xtileset = {'W{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)} | \
           {'E{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)}

for x in xtileset:
    for y in ytileset:
        tile = "http://mirrors.iplantcollaborative.org/earthenv_dem_data/EarthEnv-DEM90/" \
               "EarthEnv-DEM90_{0}{1}.tar.gz".format(y, x)
        outtile = os.path.join(ee_outdir, os.path.split(tile)[1])
        try:
            if not (os.path.exists(os.path.splitext(outtile)[0]) or \
                    os.path.exists("{}.bil".format(os.path.splitext(os.path.splitext(outtile)[0])[0]))):
                dlfile(url=tile, outpath=ee_outdir, ignore_downloadable = True, outfile=os.path.split(outtile)[1])
                print('Deleting {}...'.format(outtile))
                os.remove(outtile)
            else:
                print('{} was already processed...'.format(outtile))
        except:
            traceback.print_exc()
        del tile

ee_tarlist = getfilelist(ee_outdir, '.*[.]tar$', gdbf=False, nongdbf=True)
for f in ee_tarlist:
    print(f)
    if not os.path.exists('{}.bil'.format(os.path.splitext(f)[0])):
        with tarfile.open(f) as ftar:
            ftar.extractall(ee_outdir)  # specify which folder to extract to
"""