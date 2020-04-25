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
import math
import numpy as np
from osgeo import gdal
from osgeo import gdal_array

#pip install pyproj==1.9.6 owslib==0.18 - 0.19 dropped python 2.7
from owslib.wcs import WebCoverageService  # OWSlib module to access WMS services from SDAT

#Folder structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
datdir = os.path.join(rootdir, 'HydroATLAS', 'data')
resdir = os.path.join(rootdir, 'HydroATLAS' , 'results')
glad_dir = os.path.join(datdir, 'GLAD')
pathcheckcreate(glad_dir)

#------------------------------- Download HYSOGS250 m data -------------------------------------------------------------
HYSOGSdir = os.path.join(datdir, 'HYSOGS')
pathcheckcreate(HYSOGSdir)
HYSOGSresdir = os.path.join(resdir, 'HYSOGS')
pathcheckcreate(HYSOGSresdir)
HYSOGmosaic = os.path.join(HYSOGSresdir, 'HYSOGS_mosaic')
HYSOGnodata = os.path.join(HYSOGSresdir, 'HYSOGS_nodata')

#DOI of dataset is https://doi.org/10.3334/ORNLDAAC/1566
#This last number is the dataset ID to be used for search in WebMapService
sdatwcs = WebCoverageService('https://webmap.ornl.gov/ogcbroker/wcs')
print(str(len(sdatwcs.contents)) + ' layers found from ' + sdatwcs.identification.title)
# filter layers
hysog_wcsid = filter(lambda x: x.startswith('1566_'), sdatwcs.contents)[0]
print(hysog_wcsid)

HYSOGbblist = divbb(bbox=sdatwcs[hysog_wcsid].boundingBoxWGS84,
                    res=sdatwcs[hysog_wcsid].grid.offsetvectors[0][0],
                    divratio=10)
HYSOGoutlist = ['{0}_{1}.tif'.format(os.path.join(HYSOGSdir, 'HYSOGS'), i)
                for i in xrange(0, len(HYSOGbblist))]
if not all([os.path.isfile(i) for i in HYSOGoutlist]):
    x=0
    for bb in HYSOGbblist:
        #print(bb)
        outtile = HYSOGoutlist[x]
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
# for tilepath in HYSOGoutlist:
#     print(tilepath)
#     tiledat = gdal_array.LoadFile(tilepath)
#     if tiledat.max() == 0:
#         HYSOGoutlist.remove(tilepath)

#mosaic them
arcpy.MosaicToNewRaster_management(HYSOGoutlist, output_location=HYSOGSresdir,
                                   raster_dataset_name_with_extension= 'HYSOGS_mosaic',
                                   pixel_type= '8_BIT_UNSIGNED',
                                   number_of_bands = 1,
                                   mosaic_colormap_mode = 'FIRST')

Con()- Raster(os.path.join(HYSOGSresdir, 'HYSOGS_mosaic')


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





#-----------------------------------------------------------------------------------------------------------------------
#Extra stuff
# The user credentials that will be used to authenticate access to the data
with open("configs.json") as json_data_file: #https://martin-thoma.com/configuration-files-in-python/
    authdat = json.load(json_data_file)

# The url of the file we wish to retrieve
urlHYSOG = "https://daac.ornl.gov/daacdata/global_soil/Global_Hydrologic_Soil_Group/data/HYSOGs250m.tif?_ga=2.3478517.499069434.1587733359-1920916281.1587510472"

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

request = urllib2.Request(urlHYSOG)
response = urllib2.urlopen(request)
list(response.info())
response.info()['content-type']

outHYSOGS = os.path.join(datdir, 'HYSOGS', 'HYSOGS')
pathcheckcreate(os.path.split(outHYSOGS)[0])

if out is None:
    return response.read().decode('utf-8')
else:
import shutil
shutil.copyfileobj(response, outHYSOGS)


with open(outHYSOGS, "wb") as local_file:
    local_file.write(response.read())
# Unzip downloaded file
try:
    unzip(response + '.zip')
except:
    z = zipfile.ZipFile(io.BytesIO(response.content))
    if isinstance(z, zipfile.ZipFile):
        z.extractall(os.path.split(outHYSOGS)[0])

# Print out the result (not a good idea with binary data!)

body = response.read()
print body




