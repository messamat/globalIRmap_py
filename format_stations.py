#Purpose: Format GRDC station data and associated environmental variables extract from HydroATLAS

import arcpy
import os
import re

arcpy.env.overwriteOutput = True

#Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

#Parameters
wgs84 = arcpy.SpatialReference(4326)
arcpy.env.overwriteOutput = 'True'
arcpy.env.qualifiedFieldNames = 'False'

#Input variables
GRDCstations = os.path.join(datdir, 'GRDC_curated', 'high_qual_daily_stations.csv')
riveratlas = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
hydromask = os.path.join(datdir, 'HydroATLAS', 'Masks_20200525','hydrosheds_landmask_15s.gdb', 'hys_land_15s')

#Output variables
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
if not arcpy.Exists(outgdb):
    arcpy.CreateFileGDB_management(os.path.split(outgdb)[0], os.path.split(outgdb)[1])
GRDCp = os.path.join(outgdb, 'GRDCstations')
GRDCpjoin = os.path.join(outgdb, 'GRDCstations_riverjoin')
riveratlas_csv = os.path.join(resdir, 'RiverATLAS_v10tab.csv')
GRDCp_aeqd = os.path.join(outgdb, "GRDCstations_aeqd")
GRDCbuf = os.path.join(outgdb, 'GRDCstations_buf50k')

#Create points for GRDC stations
stations_coords = {row[0]:[row[1], row[2]]
                   for row in arcpy.da.SearchCursor(GRDCstations, ['GRDC_NO', 'LONG_NEW', 'LAT_NEW'])}
arcpy.CreateFeatureclass_management(os.path.split(GRDCp)[0], os.path.split(GRDCp)[1],
                                    geometry_type='POINT', spatial_reference=wgs84)
arcpy.AddField_management(GRDCp, 'GRDC_NO', 'TEXT')

with arcpy.da.InsertCursor(GRDCp, ['GRDC_NO', 'SHAPE@XY']) as cursor:
    for k, v in stations_coords.items():
        cursor.insertRow([k, arcpy.Point(v[0], v[1])])

#Join GRDC stations to nearest river reach in RiverAtlas
arcpy.SpatialJoin_analysis(GRDCp, riveratlas, GRDCpjoin, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                           match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                           distance_field_name='station_river_distance')

#Export attribute table of RiverATLAS with selected
arcpy.CopyRows_management(in_rows = riveratlas, out_table=riveratlas_csv)

#------------------------------- Create grid for prediction error mapping ----------------------------------------------
#Buffer gauging stations
hydrodir_list = {}
for (dirpath, dirnames, filenames) in \
        arcpy.da.Walk(os.path.join(datdir, 'HydroATLAS', 'Flow_directions_20200525', 'flow_dir_15s_by_continent.gdb'),
                      topdown=True, datatype="Any"):
    for filename in filenames:
        if re.search('.*dir_15s$', filename):
            #print(filename)
            hydrodir_list[re.search('^[a-z]*(?=_dir_15s$)', filename).group()] = os.path.join(dirpath, filename)

continents = hydrodir_list.keys()

bufmask_dict = {}
for cont in hydrodir_list.keys():
    tilepath = hydrodir_list[cont]
    bufmask_dict[cont] = os.path.join(outgdb, re.sub('dir', 'gaugemask', os.path.split(tilepath)[1]))
    tiledesc = arcpy.Describe(tilepath)

    if not arcpy.Exists(bufmask_dict[cont]):
        print('Processing {}...'.format(bufmask_dict[cont]))
        arcpy.MakeFeatureLayer_management(GRDCp, 'grdcp')
        arcpy.SelectLayerByLocation_management('grdcp', overlap_type='INTERSECT', select_features=tiledesc.extent.polygon)
        arcpy.CopyFeatures_management('grdcp', bufmask_dict[cont])

    arcpy.Buffer_analysis(bufmask_dict[cont], GRDCbuf, buffer_distance_or_field= '50 kilometers',
                          dissolve_option='ALL', method='GEODESIC')
    #Tile by continent
    arcpy.env.snapRaster = hydromask
    arcpy.env.cellSize = hydromask
    outbufras = os.path.join(resdir, 'GRDCstations_bufras{}.tif'.format(cont))
    print('Processing {}...'.format(outbufras))
    arcpy.PolygonToRaster_conversion(GRDCbuf, value_field='OBJECTID', out_rasterdataset=outbufras,
                                     cellsize=arcpy.Describe(hydromask).meanCellWidth)

    outbufrasproj = os.path.join(resdir, 'GRDCstations_bufrasproj{}.tif'.format(cont))
    print('Processing {}...'.format(outbufrasproj))
    arcpy.ProjectRaster_management(outbufras, outbufrasproj,
                                   out_coor_system= arcpy.SpatialReference(54032), resampling_type='NEAREST',
                                   cell_size=500)

#Export hydrosheds mask

arcpy.CopyRaster_management(hydromask, out_rasterdataset=os.path.join(resdir, 'hys_land_15s.tif'))

