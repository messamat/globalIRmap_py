import os
import arcpy
from arcpy.sa import *
import sys
import re
import math
from utility_functions import *
import numpy as np
import cProfile

arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True

# Set up dir structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\globalIRmap')[0]
datdir = os.path.join(rootdir, 'HydroATLAS', 'data')
resdir = os.path.join(rootdir, 'HydroATLAS', 'results')
scratchdir = os.path.join(rootdir, 'HydroATLAS', 'scratch')
scratchgdb = os.path.join(scratchdir, 'scratch.gdb')
pathcheckcreate(scratchgdb)
arcpy.env.scratchWorkspace = scratchgdb

hydrotemplate = os.path.join(resdir, 'HydroSHEDS', 'mask.gdb', 'af_mask_15s')  # Grab HydroSHEDS layer for one continent as template

glad_dir = os.path.join(datdir, 'GLAD')
gladresgdb = os.path.join(resdir, 'glad.gdb')
pathcheckcreate(gladresgdb)

#GLAD values mean the following
#0: NoData
#1: Land
#2: Permanent water:
#3: Stable seasonal
#4: Water gain
#5: Water loss
#6: Dry period
#7: Wet period
#8: High frequency
#10: Probable land
#11: Probable water
#12: Sparse data? or land. Exclude

rawtilelist = getfilelist(glad_dir, 'class99_18.*[.]tif$')

#Check aggregation ratio
cellsize_ratio = arcpy.Describe(hydrotemplate).meanCellWidth / arcpy.Describe(rawtilelist[0]).meanCellWidth
print('Aggregating DEM by cell size ratio of {0} would lead to a difference in resolution of {1} mm'.format(
    math.floor(cellsize_ratio),
    11100000 * (arcpy.Describe(hydrotemplate).meanCellWidth - math.floor(cellsize_ratio) * arcpy.Describe(
        rawtilelist[0]).meanCellWidth)
))
# Make sure that the cell size ratio is a multiple of the number of rows and columns in DEM tiles to not have edge effects
float(arcpy.Describe(rawtilelist[0]).height) / math.floor(cellsize_ratio)
float(arcpy.Describe(rawtilelist[0]).width) / math.floor(cellsize_ratio)

hysogagg_dict = {}
for tile in rawtilelist:
    print(tile)
    outaggk = re.sub('(^.*class99_18_)|([.]tif)', '', tile)
    hysogagg_dict[outaggk] = os.path.join(gladresgdb, '{0}_agg'.format(os.path.splitext(os.path.split(tile)[1])[0]))

    # Get unique categorical values
    if not arcpy.Exists(hysogagg_dict[outaggk]):
        if not Raster(tile).hasRAT and arcpy.Describe(tile).bandCount == 1:  # Build attribute table if doesn't exist
            try:
                arcpy.BuildRasterAttributeTable_management(tile)  # Does not work
            except Exception:
                e = sys.exc_info()[1]
                print(e.args[0])
                arcpy.DeleteRasterAttributeTable_management(tile)

        gladvals = {row[0] for row in arcpy.da.SearchCursor(tile, 'Value')}

        if len(gladvals) == 1:  # If only one value across entire tile
            if list(gladvals)[0] == 0:  # And that value is NoData
                print('Tile only has NoData values, deleting...')
                arcpy.Delete_management(tile)  # Delete tile

        else:
            # Divide and aggregate each band
            if compsr(tile, hydrotemplate):  # Make sure that share spatial reference with HydroSHEDS
                print('Divide into {1} bands and aggregate by rounded value of {2}'.format(
                    tile, len(gladvals), math.floor(cellsize_ratio)))
                arcpy.CompositeBands_management(in_rasters=catdivagg_list(inras=tile,
                                                                          vals=gladvals,
                                                                          exclude_list=[0, 12],
                                                                          aggratio=math.floor(cellsize_ratio)),
                                                out_raster=hysogagg_dict[outaggk])

#Differentiate sea water from inland water



#Mosaick tiles


#Resample to snap
hysog_500 = {}
for cont in hydrodir_list:
    hysog_500[cont] = os.path.join(hysogresgdb, 'hysogs_{}_500m'.format(cont))
    if not arcpy.Exists(hysog_500):
        print(hysog_500[cont])
        # Resample with nearest cell assignment
        arcpy.env.extent = arcpy.env.snapRaster = arcpy.env.cellSize = arcpy.env.mask = hydromask_dict[cont]
        arcpy.Resample_management(hysogagg, hysog_500[cont], cell_size=arcpy.env.cellSize, resampling_type='NEAREST')

