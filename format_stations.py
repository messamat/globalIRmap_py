#Purpose: Format station data and associated environmental variables extract from HydroATLAS

import arcpy
import os
import re
from utility_functions import *

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
riveratlas = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
riveratlasv11 = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v11.gdb', 'RiverATLAS_v11')
hydromask = os.path.join(datdir, 'HydroATLAS', 'Masks_20200525','hydrosheds_landmask_15s.gdb', 'hys_land_15s')

#Output variables
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
if not arcpy.Exists(outgdb):
    arcpy.CreateFileGDB_management(os.path.split(outgdb)[0], os.path.split(outgdb)[1])
riveratlas_csv = os.path.join(resdir, 'RiverATLAS_v10tab.csv')
riveratlasv11_csv = os.path.join(resdir, 'RiverATLAS_v11tab.csv')

#grdc stations
grdcstations = os.path.join(datdir, 'grdc_curated', 'high_qual_daily_stations.csv')
grdcp = os.path.join(outgdb, 'grdcstations')
grdcpjoin = os.path.join(outgdb, 'grdcstations_riverjoin')
grdcp_aeqd = os.path.join(outgdb, "grdcstations_aeqd")
grdcbuf = os.path.join(outgdb, 'grdcstations_buf50k')
gaugebufdiss = os.path.join(outgdb, grdcbuf + 'diss')

#GSIM stations
gsimdatdir = os.path.join(datdir, 'GSIM')
gsimmeta = os.path.join(gsimdatdir, 'GSIM_metadata', 'GSIM_catalog', 'GSIM_metadata.csv')
gsimsuspect = os.path.join(gsimdatdir, 'GSIM_metadata', 'GSIM_catalog', 'GSIM_suspect_coordinates_stations.csv')
gsiminddir = os.path.join(gsimdatdir, 'GSIM_indices', 'TIMESERIES', 'monthly')

gsimresdir = os.path.join(resdir, 'GSIM')
gsimresgdb = os.path.join(gsimresdir, 'GSIM.gdb')
pathcheckcreate(gsimresgdb)
gsimmeta_format = os.path.join(gsimresgdb, 'GSIM_metadata_format')
gsimpraw = os.path.join(gsimresgdb, 'GSIMstations_raw')
gsimpmeta = os.path.join(gsimresgdb, 'GSIMstations_metadata')


######################################## FORMAT GRDC STATIONS ##########################################################
#Create points for grdc stations
if not arcpy.Exists(grdcp):
    print('Create points for grdc stations')
    stations_coords = {row[0]:[row[1], row[2]]
                       for row in arcpy.da.SearchCursor(grdcstations, ['GRDC_NO', 'LONG_NEW', 'LAT_NEW'])}
    arcpy.CreateFeatureclass_management(os.path.split(grdcp)[0], os.path.split(grdcp)[1],
                                        geometry_type='POINT', spatial_reference=wgs84)
    arcpy.AddField_management(grdcp, 'GRDC_NO', 'TEXT')

    with arcpy.da.InsertCursor(grdcp, ['GRDC_NO', 'SHAPE@XY']) as cursor:
        for k, v in stations_coords.items():
            cursor.insertRow([k, arcpy.Point(v[0], v[1])])

#Join grdc stations to nearest river reach in RiverAtlas
if not arcpy.Exists(grdcpjoin):
    print('Join grdc stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(grdcp, riveratlas, grdcpjoin, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                               distance_field_name='station_river_distance')

#Export attribute table of RiverATLAS with selected
if not arcpy.Exists(riveratlas_csv):
    arcpy.CopyRows_management(in_rows = riveratlas, out_table=riveratlas_csv)
if not arcpy.Exists(riveratlasv11_csv):
    arcpy.CopyRows_management(in_rows=riveratlasv11, out_table=riveratlasv11_csv)

######################################## FORMAT GSIM STATIONS ##########################################################
#Replace points with underscores in metadata table and import into geodatabase with correct field types
if not arcpy.Exists(gsimmeta_format):
    metarawtab = pd.read_csv(gsimmeta)
    metarawtab.columns = [re.sub('[.]', '_', s) for s in metarawtab.columns]
    x = metarawtab.reset_index()
    z = np.rec.fromrecords(x, names=x.columns.tolist())
    arcpy.da.NumPyArrayToTable(z, gsimmeta_format)

#Create points for GSIM stations
if not arcpy.Exists(gsimpmeta):
    arcpy.MakeXYEventLayer_management(gsimmeta_format, in_x_field='longitude', in_y_field='latitude',
                                      out_layer='gsimpmetatest', spatial_reference=wgs84)
    arcpy.SaveToLayerFile_management('gsimpmetatest', os.path.join(gsimresgdb, 'test'))
    arcpy.CopyFeatures_management(os.path.join(gsimresgdb, 'test.lyr'), gsimpmeta)

#Only keep those with drainage area information and at least 10 years of data (with max 20 years of missing data per year)
gsimpsub = os.path.join(gsimresgdb, 'GSIMSstations_sub')
if not arcpy.Exists(gsimpsub):
    arcpy.MakeFeatureLayer_management(gsimpmeta, 'gsimpmetalyr',  where_clause= '(area > 0) AND (number_available_days > 3450)')
    arcpy.CopyFeatures_management('gsimpmetalyr', gsimpsub)

    #Delete those that are already in the curated GRDC database or with suspect coordinates
    grdcids = [row[0] for row in arcpy.da.SearchCursor(grdcp, 'GRDC_NO')]
    suspectids= pd.read_csv(gsimsuspect)['reference.no'].to_list()
    x=0
    with arcpy.da.UpdateCursor(gsimpsub, ['grdb_no', 'reference_no']) as cursor:
        for row in cursor:
            if row[0] is not None:
                if str(int(row[0])) in grdcids:
                    cursor.deleteRow()
                    x+=1
                    continue
            if row[1] in suspectids:
                cursor.deleteRow()
                x += 1
    print('Deleted {} gauges in GSIM that were already in the curated GRDC datasets or had suspect coordinates'.format(x))
    del x
    del grdcids

#Merge all monthly discharge indices time series for gauges in the initial selection
gsimsubno = [row[0] for row in arcpy.da.SearchCursor(gsimpsub, 'gsim_no')]
def format_gsimmontab(tab, filterlist):
    tabno = os.path.splitext(os.path.split(tab)[1])[0]
    print(tabno)
    if tabno in filterlist:
        tsdf = pd.read_csv(tab, sep=',\\t', header=21)
        tsdf.columns = [re.sub('"', '', s) for s in tsdf.columns]
        tsdf['gsim_no'] = tabno
        return(tsdf)

gsimsub_modf = pd.concat([format_gsimmontab(file, filterlist=gsimsubno) for file in getfilelist(gsiminddir)],
                         axis=0) \
    .sort_index()

gsimsub_modf['year'] = gsimsub_modf['date']
gsimsub_modf['month'] =
gsimsub_modf['day'] =

#Only keep those that have at least 10 years of data, only counting years with less than 20 days of missing data,
#and that have average discharge of at least 0.01 m3/s OR a reported drainage area >= 9 km2




#Spatial Join to HydroATLAS


#------------------------------- Create grid for prediction error mapping ----------------------------------------------
#Buffer gauging stations
arcpy.Buffer_analysis(grdcp, grdcbuf, buffer_distance_or_field='50 kilometers',
                      dissolve_option='ALL', method='GEODESIC')
arcpy.Dissolve_management(grdcbuf, gaugebufdiss, multi_part='FALSE')

#Split
bufgdb = os.path.join(resdir, 'gaugebuf.gdb')
arcpy.CreateFileGDB_management(os.path.split(bufgdb)[0], os.path.split(bufgdb)[1])
arcpy.AddField_management(gaugebufdiss, 'UID', 'SHORT')
with arcpy.da.UpdateCursor(gaugebufdiss , 'UID') as cursor:
    x = 0
    for row in cursor:
        print(x)
        row[0] = x
        x+=1
        cursor.updateRow(row)

arcpy.SplitByAttributes_analysis(gaugebufdiss, bufgdb, 'UID')

bufrasdir = os.path.join(resdir, 'bufrasdir', 'bufras_{}proj')
os.mkdir(bufrasdir)
arcpy.env.workspace = bufgdb
for buf in arcpy.ListFeatureClasses('T*', feature_type='Polygon'):
    arcpy.env.snapRaster = hydromask
    arcpy.env.cellSize = hydromask
    outbufras = os.path.join(bufgdb, 'bufras_{}'.format(buf))
    print('Processing {}...'.format(outbufras))
    if not arcpy.Exists(outbufras):
        arcpy.PolygonToRaster_conversion(os.path.join(bufgdb, buf), value_field='OBJECTID', out_rasterdataset=outbufras,
                                         cellsize=arcpy.Describe(hydromask).meanCellWidth)

    outbufrasproj = os.path.join(bufrasdir, 'bufras_{}proj.tif'.format(buf))
    print('Processing {}...'.format(outbufrasproj))
    if not arcpy.Exists(outbufrasproj):
        arcpy.ResetEnvironments()
        arcpy.ProjectRaster_management(in_raster=outbufras, out_raster=outbufrasproj,
                                       out_coor_system=arcpy.SpatialReference(54032), resampling_type='NEAREST',
                                       cell_size=500)

########################################################################################################################
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
        arcpy.MakeFeatureLayer_management(grdcp, 'grdcp')
        arcpy.SelectLayerByLocation_management('grdcp', overlap_type='INTERSECT', select_features=tiledesc.extent.polygon)
        arcpy.CopyFeatures_management('grdcp', bufmask_dict[cont])

    arcpy.Buffer_analysis(bufmask_dict[cont], grdcbuf, buffer_distance_or_field= '50 kilometers',
                          dissolve_option='ALL', method='GEODESIC')
    #Tile by continent
    arcpy.env.snapRaster = hydromask
    arcpy.env.cellSize = hydromask
    outbufras = os.path.join(resdir, 'grdcstations_bufras{}.tif'.format(cont))
    print('Processing {}...'.format(outbufras))
    arcpy.PolygonToRaster_conversion(grdcbuf, value_field='OBJECTID', out_rasterdataset=outbufras,
                                     cellsize=arcpy.Describe(hydromask).meanCellWidth)

    outbufrasproj = os.path.join(resdir, 'grdcstations_bufrasproj{}.tif'.format(cont))
    print('Processing {}...'.format(outbufrasproj))
    arcpy.ProjectRaster_management(outbufras, outbufrasproj,
                                   out_coor_system= arcpy.SpatialReference(54032), resampling_type='NEAREST',
                                   cell_size=500)

#Export hydrosheds mask

arcpy.CopyRaster_management(hydromask, out_rasterdataset=os.path.join(resdir, 'hys_land_15s.tif'))

