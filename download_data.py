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

#pip install pyproj==1.9.6 owslib==0.18 - 0.19 dropped python 2.7
from owslib.wcs import WebCoverageService  # OWSlib module to access WMS services from SDAT

#Folder structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
datdir = os.path.join(rootdir, 'HydroATLAS', 'data')
resdir = os.path.join(rootdir, 'HydroATLAS' , 'results')
glad_dir = os.path.join(datdir, 'GLAD')
pathcheckcreate(glad_dir)

# ------------------------------- Download EarthEnv-DEM90 data -------------------------------------------------------------
ee_outdir = os.path.join(datdir, 'earthenv')
pathcheckcreate(ee_outdir)

# Parse HTML to get all available layers
ee_https = "https://www.earthenv.org/DEM.html" #Doesn't work, return 403 with urlibb2

ytileset = {'N{}'.format(x) for x in xrange(0, 85, 5)} | \
           {'S{}'.format(x) for x in xrange(0, 60, 5)}
xtileset = {'W{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)} | \
           {'E{}'.format(str(x).zfill(3)) for x in xrange(0, 185, 5)}

for x in xtileset:
    for y in ytileset:
        tile = "http://mirrors.iplantcollaborative.org/earthenv_dem_data/EarthEnv-DEM90/" \
               "EarthEnv-DEM90_{0}{1}.tar.gz".format(y, x)
        outtile = os.path.join(ee_outdir, os.path.split(tile)[1])
        try:
            if not os.path.exists(os.path.splitext(outtile)[0]):
                dlfile(url=tile, outpath=ee_outdir, ignore_downloadable = True, outfile=os.path.split(outtile)[1])
                print('Deleting {}...'.format(outtile))
                os.remove(outtile)
        except:
            traceback.print_exc()
        del tile

# ------------------------------- Download SoilGrids250 m data -------------------------------------------------------------
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

#------------------------------- Download MODIS 250 m land and water mask to enhance SoilGrids -------------------------
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
with open("configs.json") as json_data_file:  # https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)

# Download all layers of interest
for lyrurl in [urlparse.urljoin(mod44w_https, link.get('href')) for link in
               mod44w_soup.findAll('a', attrs={'href': re.compile("MOD44W.*[.]hdf([.]xml)*")})]:

        dlfile(url=lyrurl, outpath=mod44w_outdir, outfile=os.path.split(lyrurl)[1], fieldnames=None,
               loginprompter="https://urs.earthdata.nasa.gov",
               username=authdat['earthdata']['username'], password=authdat['earthdata']['password'])

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


# ----------------------------- Download GLAD data ----------------------------------------------------------------------
glad_dtype = "class99_18.tif"
gsutil_ls_cmd = "gsutil ls gs://earthenginepartners-hansen/water/*/{}".format(glad_dtype)
glad_cloudout = subprocess.check_output(gsutil_ls_cmd)
glad_cloudlist = glad_cloudout.split('\n')

for tile in glad_cloudlist:
    out_tilen = os.path.join(glad_dir, '{0}_{1}.tif'.format(
        os.path.splitext(glad_dtype)[0], tileroot))
    if not os.path.isfile(out_tilen):
        print(tile)
        tileroot = os.path.split(os.path.split(tile)[0])[1]
        gsutil_cp_cmd = "gsutil cp {0} {1}".format(tile, glad_dir)
        subprocess.check_output(gsutil_cp_cmd)
        os.rename(os.path.join(glad_dir, glad_dtype), out_tilen)












########################################################################################################################
########################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------
#Extra stuff
# The user credentials that will be used to authenticate access to the data
with open("configs.json") as json_data_file: #https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)

# The url of the file we wish to retrieve
urlhysog = "https://daac.ornl.gov/daacdata/global_soil/Global_Hydrologic_Soil_Group/data/HYSOGs250m.tif?_ga=2.3478517.499069434.1587733359-1920916281.1587510472"

# Create a password manager to deal with the 401 reponse that is returned from Earthdata Login
password_manager = urllib2.HTTPPasswordMgrWithDefaultRealm()
password_manager.add_password(None, "https://urs.earthdata.nasa.gov",
                              authdat['earthdata']['username'], authdat['earthdata']['password'])

# Create a cookie jar for storing cookies. This is used to store and return the session cookie given to use by
# the data server (otherwise it will just keep sending us back to Earthdata Login to authenticate).
# Ideally, we should use a file based cookie jar to preserve cookies between runs. This will make it much more efficient.
cookie_jar = CookieJar()

# Install all the handlers.
opener = urllib2.build_opener(
    urllib2.HTTPBasicAuthHandler(password_manager),
    # urllib2.HTTPHandler(debuglevel=1),    # Uncomment these two lines to see
    # urllib2.HTTPSHandler(debuglevel=1),   # details of the requests/responses
    urllib2.HTTPCookieProcessor(cookie_jar))
urllib2.install_opener(opener)

# Create and submit the request. There are a wide range of exceptions that
# can be thrown here, including HTTPError and URLError. These should be
# caught and handled.

request = urllib2.Request(urlhysog)
response = urllib2.urlopen(request)
list(response.info())
response.info()['content-type']

outhysog = os.path.join(datdir, 'hysog', 'hysog')
pathcheckcreate(os.path.split(outhysog)[0])

if out is None:
    return response.read().decode('utf-8')
else:
import shutil
shutil.copyfileobj(response, outhysog)


with open(outhysog, "wb") as local_file:
    local_file.write(response.read())
# Unzip downloaded file
try:
    unzip(response + '.zip')
except:
    z = zipfile.ZipFile(io.BytesIO(response.content))
    if isinstance(z, zipfile.ZipFile):
        z.extractall(os.path.split(outhysog)[0])

# Print out the result (not a good idea with binary data!)

body = response.read()
print body




