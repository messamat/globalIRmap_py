import os
import arcpy
from arcpy.sa import *
import sys
import re
from utility_functions import *
import numpy as np

arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True

#Set up dir structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
datdir = os.path.join(rootdir, 'HydroATLAS', 'data')
resdir = os.path.join(rootdir, 'HydroATLAS', 'results')

#---------------------------------------- Create HydroSHEDS landmask ---------------------------------------------------
hydromaskdir = os.path.join(resdir, 'HydroSHEDS', 'mask')
hydromaskgdb = os.path.join(resdir, 'HydroSHEDS', 'mask.gdb')
pathcheckcreate(hydromaskdir)
if not arcpy.Exists(hydromaskgdb):
    arcpy.CreateFileGDB_management(os.path.split(hydromaskgdb)[0], os.path.split(hydromaskgdb)[1])


hydrodir_list = []
for (dirpath, dirnames, filenames) in \
        arcpy.da.Walk(os.path.join(datdir, 'HydroSHEDS'), topdown=True, datatype="RasterDataset", type="GRID"):
    for filename in filenames:
        if re.search('.*dir_15s$', filename):
            #print(filename)
            hydrodir_list.append(os.path.join(dirpath, filename))


for tilepath in hydrodir_list:
    outmask = os.path.join(hydromaskgdb, re.sub('dir', 'mask', os.path.split(tilepath)[1]))
    #if not arcpy.Exists(outmask):
    print('Processing {}...'.format(outmask))
    arcpy.env.XYResolution = "0.0000000000000001 degrees"
    arcpy.env.extent = arcpy.env.snapRaster = tilepath
    #Because HydroSHEDS cell size is set with 16 digits while arcpy.env assignment by layer only handles 12,
    #explicity use arcpy.Describe
    arcpy.env.cellSize = arcpy.Describe(hydrodir_list[0]).meanCellWidth
    print('%.17f' % float(arcpy.env.cellSize))
    print(arcpy.env.extent)
    print(arcpy.Describe(tilepath).Extent)

    try:
        #Create a 1-NoData mask
        #arcpy.CopyRaster_management(Con(~IsNull(Raster(tilepath)), 1).save(outmask), outmask, pixel_type='8_BIT_UNSIGNED')

        #Check whether everything is the same
        arcpy.RasterCompare_management(outmask, tilepath, compare_type='RASTER_DATASET',
                                       continue_compare='CONTINUE_COMPARE',
                                       out_compare_file=os.path.join(resdir,
                                                                     re.sub('dir_15s', 'compare',
                                                                            os.path.split(tilepath)[1])))
        print(arcpy.Describe(outmask).Extent)
        print(arcpy.Describe(outmask).meanCellWidth)


    except:
        arcpy.ResetEnvironments()

