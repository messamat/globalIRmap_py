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

#Input dir and lyrs
hysogresgdb = os.path.join(resdir, 'hysogs.gdb')
hysogmosaic = os.path.join(hysogresgdb, 'hysogs_mosaic')


# Aggregate each band by cell size ratio
def catdivagg_list(inras, vals, exclude_list, aggratio):
    return (
        [Aggregate((Con(Raster(inras) == v, 1)),
                   aggratio, aggregation_type='SUM', extent_handling='EXPAND', ignore_nodata='DATA')
         for v in vals if v not in exclude_list])

#---------------------------------------- Create HydroSHEDS landmask ---------------------------------------------------
hydromaskdir = os.path.join(resdir, 'HydroSHEDS', 'mask')
hydromaskgdb = os.path.join(resdir, 'HydroSHEDS', 'mask.gdb')
pathcheckcreate(hydromaskdir)
if not arcpy.Exists(hydromaskgdb):
    arcpy.CreateFileGDB_management(os.path.split(hydromaskgdb)[0], os.path.split(hydromaskgdb)[1])

hydrodir_list = {}
for (dirpath, dirnames, filenames) in \
        arcpy.da.Walk(os.path.join(datdir, 'HydroSHEDS'), topdown=True, datatype="RasterDataset", type="GRID"):
    for filename in filenames:
        if re.search('.*dir_15s$', filename):
            #print(filename)
            hydrodir_list[re.search('^[a-z]*(?=_dir_15s$)', filename).group()] = os.path.join(dirpath, filename)

continents = hydrodir_list.keys()

hydromask_dict = {}
for cont in hydrodir_list.keys():
    tilepath = hydrodir_list[cont]
    hydromask_dict[cont] = os.path.join(hydromaskgdb, re.sub('dir', 'mask', os.path.split(tilepath)[1]))
    tiledesc = arcpy.Describe(tilepath)

    if not arcpy.Exists(hydromask_dict[cont]):
        print('Processing {}...'.format(hydromask_dict[cont]))
        arcpy.env.extent = arcpy.env.snapRaster = tilepath
        #Because HydroSHEDS cell size is set with 16 digits while arcpy.env assignment by layer only handles 12,
        #explicity use arcpy.Describe — still doesn't work unless use file gdb
        arcpy.env.XYResolution = "0.0000000000000001 degrees"
        arcpy.env.cellSize = arcpy.Describe(tilepath).meanCellWidth
        print('%.17f' % float(arcpy.env.cellSize))

        try:
            #Create a 1-NoData mask
            arcpy.CopyRaster_management(Con(~IsNull(Raster(tilepath)), 1),
                                        hydromask_dict[cont], pixel_type='8_BIT_UNSIGNED')

            #Check whether everything is the same
            maskdesc = arcpy.Describe(hydromask_dict[cont])
            print('Equal extents? {}'.format(maskdesc.Extent.JSON == tiledesc.Extent.JSON))
            print('Equal cell size? {}'.format(maskdesc.meanCellWidth == tiledesc.meanCellWidth))
            print('Same Spatial Reference? {}'.format(compsr(tilepath, hydromask_dict[cont])))

        except Exception:
            print("Exception in user code:")
            traceback.print_exc(file=sys.stdout)
            arcpy.ResetEnvironments()

hydrotemplate = hydromask_dict[hydromask_dict.keys()[0]]  # Grab HydroSHEDS layer for one continent as template

#----------------------------------------- Format MODIS 250m water mask ------------------------------------------------
mod44w_outdir = os.path.join(datdir, 'mod44w')
mod44w_resgdb = os.path.join(resdir, 'mod44w.gdb')
mod44w_mosaic = os.path.join(mod44w_resgdb, 'mod44w_mosaic')
mod44w_QAmosaic = os.path.join(mod44w_resgdb, 'mod44w_QAmosaic')
pathcheckcreate(mod44w_resgdb)

for tile in getfilelist(mod44w_outdir, '.*[.]hdf$'):
    #Generate land-water mask
    outgrid = os.path.join(mod44w_resgdb, re.sub('[.]', '_', os.path.splitext(os.path.split(tile)[1])[0]))
    outgrid_wgs = '{}_wgs84'.format(outgrid)
    if not arcpy.Exists(outgrid):
        print(outgrid)
        arcpy.ExtractSubDataset_management(in_raster=tile, out_raster=outgrid, subdataset_index=0)
    if not arcpy.Exists(outgrid_wgs):
        print(outgrid_wgs)
        arcpy.env.snapRaster = hysogmosaic
        arcpy.ProjectRaster_management(outgrid, outgrid_wgs, out_coor_system=hydrotemplate,
                                       resampling_type='NEAREST', cell_size=arcpy.Describe(hysogmosaic).meanCellWidth)

    #Generate QA mask that also identifies water
    outQA = os.path.join(mod44w_resgdb, "QA_{}".format(re.sub('[.]', '_', os.path.splitext(os.path.split(tile)[1])[0])))
    outQA_wgs = '{}_wgs84'.format(outQA)
    if not arcpy.Exists(outQA):
        print(outQA)
        arcpy.ExtractSubDataset_management(in_raster=tile, out_raster=outQA, subdataset_index=1)
    if not arcpy.Exists(outQA_wgs):
        print(outQA_wgs)
        arcpy.env.snapRaster = hysogmosaic
        arcpy.ProjectRaster_management(outQA, outQA_wgs, out_coor_system=hydrotemplate,
                                       resampling_type='NEAREST', cell_size=arcpy.Describe(hysogmosaic).meanCellWidth)
    arcpy.ResetEnvironments()

mod44w_QAmosaic = os.path.join(mod44w_resgdb, 'mod44w_QAmosaic')
arcpy.MosaicToNewRaster_management(
    getfilelist(dir=mod44w_resgdb, repattern='^QA.*_wgs84$', gdbf=True, nongdbf=False),
    output_location=mod44w_resgdb,
    raster_dataset_name_with_extension=os.path.split(mod44w_QAmosaic)[1],
    number_of_bands=1)

#Turn all land pixels to NoData, then run Euclidean Allocation to 


# ---------------------------------------- Format hysogs ---------------------------------------------------
hysogagg = os.path.join(hysogresgdb, 'hysogs_agg')

testdat = os.path.join(hysogresgdb, 'test')

# Identify permanent ice pixels as defined by SoilGrids
#


# Identify inland water vs ocean mask
# Change the former to 0, rest as NoData



# Get unique categorical values
if not arcpy.Exists(hysogagg):
    if not Raster(hysogmosaic).hasRAT and arcpy.Describe(hysogmosaic).bandCount == 1: #Build attribute table if doesn't exist
        try:
            arcpy.BuildRasterAttributeTable_management(hysogmosaic) #Does not work
        except Exception:
            e = sys.exc_info()[1]
            print(e.args[0])
            arcpy.DeleteRasterAttributeTable_management(hysogmosaic)

    hysogvals = {row[0] for row in arcpy.da.SearchCursor(hysogmosaic, 'Value')}

    #Divide and aggregate each band
    if compsr(hysogmosaic, hydrotemplate):  # Make sure that share spatial reference with HydroSHEDS
        # Check cellsize ratio
        cellsize_ratio = arcpy.Describe(hydrotemplate).meanCellWidth / arcpy.Describe(hysogmosaic).Children[0].meanCellWidth
        print('Divide {0} into {1} bands and aggregate by rounded value of {2}'.format(
            hysogmosaic, len(hysogvals)-1, cellsize_ratio))
        arcpy.CompositeBands_management(in_rasters=catdivagg_list(inras=hysogmosaic, vals=hysogvals,
                                                                  exclude_list=[0], aggratio=round(cellsize_ratio)),
                                        out_raster=hysogagg)

hysog_500 = {}
for cont in hydrodir_list:
    hysog_500[cont] = os.path.join(hysogresgdb, 'hysogs_{}_500m'.format(cont))
    if not arcpy.Exists(hysog_500):
        print(hysog_500[cont])
        #Resample with nearest cell assignment
        arcpy.env.extent = arcpy.env.snapRaster = arcpy.env.cellSize = arcpy.env.mask = hydromask_dict[cont]
        arcpy.Resample_management(hysogagg, hysog_500[cont], cell_size=arcpy.env.cellSize, resampling_type='NEAREST')


#Notes for cleaning and formatting
#Check out the generalization toolset
#https://desktop.arcgis.com/en/arcmap/10.7/tools/spatial-analyst-toolbox/an-overview-of-the-generalization-tools.htm
#Think of using Nibble to fill-in NoData regions
#Look at Region Group.