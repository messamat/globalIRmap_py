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


# Get unique categorical values
rawtilelist = getfilelist(glad_dir, 'class99_18.*[.]tif')
hysogagg_dict = {}
for tile in rawtilelist:
    outaggk = re.sub('(^.*class99_18_)|([.]tif)', '', tile)
    hysogagg_dict[outaggk] = os.path.join(gladresgdb, '{0}_agg'.format(os.path.splitext(tile)[0]))

    if not arcpy.Exists(hysogagg_dict[outaggk]):
        if not Raster(tile).hasRAT and arcpy.Describe(tile).bandCount == 1:  # Build attribute table if doesn't exist
            try:
                arcpy.BuildRasterAttributeTable_management(tile)  # Does not work
            except Exception:
                e = sys.exc_info()[1]
                print(e.args[0])
                arcpy.DeleteRasterAttributeTable_management(tile)

        gladvals = {row[0] for row in arcpy.da.SearchCursor(tile, 'Value')}

        # Divide and aggregate each band
        if compsr(tile, hydrotemplate):  # Make sure that share spatial reference with HydroSHEDS
            # Check cellsize ratio
            cellsize_ratio = arcpy.Describe(hydrotemplate).meanCellWidth / arcpy.Describe(tile).Children[
                0].meanCellWidth
            print('Divide {0} into {1} bands and aggregate by rounded value of {2}'.format(
                hysogmosaic, len(hysogvals) - 1, cellsize_ratio))
            arcpy.CompositeBands_management(in_rasters=catdivagg_list(inras=hysogmosaic, vals=hysogvals,
                                                                      exclude_list=[0], aggratio=math.floor(cellsize_ratio)),
                                            out_raster=hysogagg)

