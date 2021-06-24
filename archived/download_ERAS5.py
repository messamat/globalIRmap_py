import os
#https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows
#https://confluence.ecmwf.int/pages/viewpage.action?pageId=139068264

#Virtual environments
#create from pycharm to make sure that all arcpy and other packages are passed on BUT
#Run virtualenv\Scripts\activate in prompt before using pip

import cdsapi
import json
import xarray as xr
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
#import cdstoolbox as ct

from utility_functions import *

# Loop over years, months, and datasets
def days_of_month(y, m):
    d0 = datetime(y, m, 1)
    d1 = (datetime(y, m, 1) + relativedelta(months=1) - relativedelta(days=1)).day
    out = list()
    while d0 < d1:
        out.append(d0.strftime('%Y-%m-%d'))
        d0 += timedelta(days=1)
    return out

if __name__ == '__main__':
    #Folder structure
    rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
    datdir = os.path.join(rootdir, '\\HydroATLAS', 'data')
    resdir = os.path.join(rootdir, '\\HydroATLAS' , 'results')

    eras_outdir = os.path.join(datdir, 'eras5')
    pathcheckcreate(eras_outdir)

    ptest = os.path.join('C:/globalIRmap/src/globalIRmap_py/test3.nc')
    os.path.exists(ptest)
    sourcef = xr.open_dataarray(ptest)
    #sourcef = xr.open_dataarray(ptest, chunks={'time': 10}).isel(time=0)

    sourcef1 = sourcef.where(sourcef>0 & ~np.isnan(sourcef), drop=True)


    sourcef_ar = sourcef.to_array(dim='tprate')
    sourcef_ar.values
    (sourcef_ar['tprate']==0).values

    #Make sure that API key is accessible to Copernicus API module
    cdsapi_idfile = os.path.join(os.path.expanduser("~"), ".cdsapirc")
    if not os.path.exists(cdsapi_idfile):
        with open("../configs.json") as json_data_file:  # https://martin-thoma.com/configuration-files-in-python/
            authdat = json.load(json_data_file)

        with open(cdsapi_idfile, 'w') as out_file:
            apiid_lines = "url: https://cds.climate.copernicus.eu/api/v2\n" \
                          "key: {0}:{1}\n" \
                          "verify: 0".format(
                authdat['copernicusAPI']['UID'], authdat['copernicusAPI']['APIkey'])
            out_file.write(apiid_lines)

    #Change working directory for downloading eras5-land data
    os.chdir(eras_outdir)

    #Launch cdsapi client
    c = cdsapi.Client()



    varstodl = 'total_precipitation'

    for y in range(1981, 2020):
        for m in range(1, 13):
            for d in range(1, (datetime(y, m+1, 1) - datetime(y, m, 1)).days + 1):
                for var in varstodl:
                    print('Downlading {0} for {1}'.format(var, d))
                    c.retrieve('reanalysis-era5-land',
                               {
                                   'format': 'netcdf',
                                   'variable': 'total_precipitation',
                                   'year': y,
                                   'month': m,
                                   'day': '01',
                                   'time': [
                                       '00:00', '01:00', '02:00',
                                       '03:00', '04:00', '05:00',
                                       '06:00', '07:00', '08:00',
                                       '09:00', '10:00', '11:00',
                                       '12:00', '13:00', '14:00',
                                       '15:00', '16:00', '17:00',
                                       '18:00', '19:00', '20:00',
                                       '21:00', '22:00', '23:00',
                                   ],
                               },
                               "{var}_{day}.nc".format(var=var, day=d)
                               )




#Copernicus Data Service toolbox
#https://cds.climate.copernicus.eu/toolbox/doc/tutorial.html#what-is-the-toolbox
#Run daily mean callculation using