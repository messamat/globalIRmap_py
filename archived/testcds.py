import cdstoolbox as ct
from datetime import datetime, timedelta
import os

#https://stackoverflow.com/questions/40001892/reading-named-command-arguments
#https://docs.python.org/2/library/argparse.html
#Exclude drizzle < 1mm/day - https://journals.ametsoc.org/doi/full/10.1175/JCLI3672.1 and https://journals.ametsoc.org/doi/full/10.1175/2010JCLI3571.1

variables = 'total_precipitation'

# Folder structure
rootdir = 'C:'
datdir = os.path.join(rootdir, 'HydroATLAS', 'data')
resdir = os.path.join(rootdir, 'HydroATLAS', 'results')
eras_outdir = os.path.join(datdir, 'eras5')
drizzle_cutoff = 5 / (3600 * 24 * 1000)  # Only keep days with at least 1 mm/day (more than drizzle)

var = 'total_precipitation'
y = '2018'
m = '12'


# d = [str(i).zfill(2) for i in range(1, (datetime(y, m + 1, 1) - datetime(y, m, 1)).days + 1)]

# Change working directory for downloading eras5-land data
@ct.application(title='Download daily prec')
@ct.output.download()
# Change working directory for downloading eras5-land data
def application(var='total_precipitation',
                y='2018',
                m='12'):
    data = ct.catalogue.retrieve(
        'reanalysis-era5-land',
        {
            'format': 'netcdf',
            'variable': var,
            'year': y,
            'month': m,
            'day': [
                '01'
            ],
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
        }
    )

    dailymean = ct.climate.daily_mean(data)
    dailysub = ct.cube.where(dailymean > drizzle_cutoff, x=dailymean, drop=True)
    return dailysub