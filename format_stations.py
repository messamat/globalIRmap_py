#Purpose: Format station data and associated environmental variables extract from HydroATLAS

from utility_functions import *

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

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
seamaskbuf3k = os.path.join(datdir, 'GLAD', 'class99_19_rsp9_buf3k1.tif')
basinatlasl05 = os.path.join(datdir, 'HydroATLAS', 'BasinATLAS_v10.gdb', 'BasinATLAS_v10_lev05')

#Output variables
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
if not arcpy.Exists(outgdb):
    arcpy.CreateFileGDB_management(os.path.split(outgdb)[0], os.path.split(outgdb)[1])

#grdc stations
grdcstations = os.path.join(datdir, 'grdc_curated', 'high_qual_daily_stations.csv')
grdc_disdatdir = os.path.join(datdir, 'GRDCdat_day')
grdcp = os.path.join(outgdb, 'grdcstations')
grdcpjoin = os.path.join(outgdb, 'grdcstations_riverjoin')
grdcpclean = os.path.join(outgdb, 'grdcstations_clean')
grdcpcleanjoin = os.path.join(outgdb, 'grdcstations_cleanjoin')
basin5grdcpjoin = os.path.join(outgdb, 'BasinATLAS_v10_lev05_GRDCstations_join')

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

gsimpsub = os.path.join(gsimresgdb, 'GSIMSstations_sub')
gsimpsub2 = os.path.join(gsimresgdb, 'GSIMSstations_sub2')
gsimpsub2join = os.path.join(gsimresgdb, 'GSIMSstations_sub2_basinjoin')
gsimpsub3= os.path.join(gsimresgdb, 'GSIMSstations_sub3')
gsimpsub3join = os.path.join(gsimresgdb, 'GSIMSstations_sub3_riverjoin')
gsimpsub3b= os.path.join(gsimresgdb, 'GSIMSstations_sub3b')
gsimpsub3joinb = os.path.join(gsimresgdb, 'GSIMSstations_sub3_riverjoinb')

grdcp_10ymin = os.path.join(gsimresgdb, 'grdcstations_10ymin')
grdcpjoinedit = os.path.join(outgdb, 'grdcstations_riverjoinedit')
gsimpsnap = os.path.join(gsimresgdb, 'GSIMSstations_riversnap')
gsimpsnapedit = os.path.join(gsimresgdb, 'GSIMSstations_riversnapedit')
gsimpsnapb= os.path.join(gsimresgdb, 'GSIMSstations_riversnapb')
gsimpsnapeditb = os.path.join(gsimresgdb, 'GSIMSstations_riversnapeditb')

gsimpsnapedit_merge = os.path.join(gsimresgdb, 'GSIMstations_riversnapedit_merge')
gsimpsnapclean = os.path.join(gsimresgdb, 'GSIMstations_riversnapclean')
gsimpcleanjoin = os.path.join(gsimresgdb, 'GSIMstations_cleanriverjoin')

######################################## FORMAT GRDC STATIONS ##########################################################
#Create points for grdc stations
if not arcpy.Exists(grdcp):
    print('Create points for grdc stations')
    stations_coords = {row[0]:[row[1], row[2], row[3]]
                       for row in arcpy.da.SearchCursor(grdcstations, ['GRDC_NO', 'LONG_NEW', 'LAT_NEW', 'AREA'])}
    arcpy.CreateFeatureclass_management(os.path.split(grdcp)[0], os.path.split(grdcp)[1],
                                        geometry_type='POINT', spatial_reference=wgs84)
    arcpy.AddField_management(grdcp, 'GRDC_NO', 'TEXT')
    arcpy.AddField_management(grdcp, 'GRDC_AREA', 'DOUBLE')

    with arcpy.da.InsertCursor(grdcp, ['GRDC_NO', 'GRDC_AREA', 'SHAPE@XY']) as cursor:
        for k, v in stations_coords.items():
            cursor.insertRow([k, v[2], arcpy.Point(v[0], v[1])])

#Join grdc stations to nearest river reach in RiverAtlas
if not arcpy.Exists(grdcpjoin):
    print('Join grdc stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(grdcp, riveratlas, grdcpjoin, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                               distance_field_name='station_river_distance')

    arcpy.AddField_management(grdcpjoin, field_name='DApercdiff', field_type='FLOAT')
    arcpy.CalculateField_management(grdcpjoin, field='DApercdiff',
                                    expression='(!GRDC_AREA!-!UPLAND_SKM!)/!UPLAND_SKM!',
                                    expression_type='PYTHON')

#Check and correct all those that are more than 50 meters OR whose |DApercdiff| > 0.10
if not arcpy.Exists(grdcpjoinedit):
    arcpy.CopyFeatures_management(grdcpjoin, grdcpjoinedit)
    arcpy.AddField_management(grdcpjoinedit, 'manualsnap_mathis', 'SHORT')
    arcpy.AddField_management(grdcpjoinedit, 'snap_comment_mathis', 'TEXT')

#Manual editing

#Delete those that were marked with -1
if not arcpy.Exists(grdcpclean):
    arcpy.CopyFeatures_management(grdcpjoinedit, grdcpclean)
    with arcpy.da.UpdateCursor(grdcpclean, ['manualsnap_mathis']) as cursor:
        for row in cursor:
            if row[0] == -1:
                cursor.deleteRow()

    #Delete useless fields and re-spatial join
    for f1 in arcpy.ListFields(grdcpclean):
        if f1.name not in [f2.name for f2 in arcpy.ListFields(grdcp)]+['manualsnap_mathis', 'snap_comment_mathis']:
            arcpy.DeleteField_management(grdcpclean, f1.name)

arcpy.SpatialJoin_analysis(grdcpclean, riveratlas, grdcpcleanjoin, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                           match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                           distance_field_name='station_river_distance')

#Extract GLAD sea mask value to check for tidal reversal as cause of intermittency
ExtractMultiValuesToPoints(in_point_features=grdcpcleanjoin, in_rasters=seamaskbuf3k,
                           bilinear_interpolate_values='NONE')

#Add coordinates
arcpy.AddGeometryAttributes_management(grdcpcleanjoin, Geometry_Properties='POINT_X_Y_Z_M')

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

    #Correct reported drainage area for US (convert from sq.miles to sq.km)
    arcpy.AddField_management(gsimpmeta, 'area_correct', 'DOUBLE')
    with arcpy.da.UpdateCursor(gsimpmeta, ['reference_db', 'area', 'area_correct']) as cursor:
        for row in cursor:
            if (row[0] == 'usgs') and (row[1] is not None):
                row[2] = row[1]*2.5899881103464
            else:
                row[2] = row[1]
            cursor.updateRow(row)

#Only keep those with drainage area information and at least 10 years of data (with max 20 years of missing data per year)
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

gsimsub_modf['year'] = gsimsub_modf['date'].str.slice(start=0, stop=4)
gsimsub_modf['month'] = gsimsub_modf['date'].str.slice(start=5, stop=7)
gsimsub_modf['day'] = gsimsub_modf['date'].str.slice(start=8, stop=10)

#Only keep those that have at least 10 years of data, only counting years with less than 20 days of missing data,
gsimsub_modf_ymiss = pd.merge(gsimsub_modf,
                              gsimsub_modf.groupby(['gsim_no','year'], as_index=False)['n.available'].sum(),
                              on=['gsim_no', 'year'], suffixes=('', '_yearsum'))
gsimsub_modf_u20miss = gsimsub_modf_ymiss[gsimsub_modf_ymiss['n.available_yearsum'] >= 345]


#Compute average discharge for years with < 20 missing days
gsimsub_annualdiss = gsimsub_modf_u20miss.groupby('gsim_no')['MEAN'].mean().reset_index()

#Get reported drainage area
gsimsub_area = pd.DataFrame.from_dict({row[0]:row[1] for row in
                                       arcpy.da.SearchCursor(gsimpsub, ['gsim_no', 'area_correct'])},
                                      orient='index').reset_index()
gsimsub_area.columns = ['gsim_no', 'area_correct']

#Compute proportion of months (within years with < 20 missing days) that have MIN==0
gsimsub_modf_u20miss_interprop = gsimsub_modf_u20miss.groupby('gsim_no')['MIN'].agg(lambda x: x.eq(0).sum()/float(x.count()))

#Merge number of years with < 20 missing days, drainage area, mean annual discharge, and intermittency proportion
gsimsub_selstats = pd.merge(gsimsub_modf_u20miss_interprop,
                            pd.merge(
                                pd.merge(gsimsub_annualdiss,
                                         gsimsub_area, on='gsim_no'),
                                gsimsub_modf_u20miss.groupby('gsim_no')['year'].nunique().reset_index(),
                                on='gsim_no', suffixes=['', '_keptcount']),
                            on='gsim_no')
gsimsub_selstats.columns = ['gsim_no', 'ir_moprop', 'myrdiss', 'DA', 'years_kept']

#Only keep those with mean dis > 0.01
gsimsub2_df = gsimsub_selstats[((gsimsub_selstats['myrdiss'] > 0.01) | (gsimsub_selstats['DA'] > 5)) &
                               (gsimsub_selstats['years_kept'] >= 10)]

#Further subset the point dataset to only keep those remaining stations
if not arcpy.Exists(gsimpsub2):
    arcpy.CopyFeatures_management(gsimpsub, gsimpsub2)
    arcpy.AddField_management(gsimpsub2, field_name='ir_moprop', field_type='FLOAT')

    with arcpy.da.UpdateCursor(gsimpsub2, ['gsim_no', 'ir_moprop']) as cursor:
        for row in cursor:
            if (gsimsub2_df['gsim_no'] == row[0]).any():
                row[1] = float(gsimsub2_df.loc[gsimsub2_df['gsim_no'] == row[0]]['ir_moprop'])
                cursor.updateRow(row)
            else:
                print('Deleting {}...'.format(row[0]))
                cursor.deleteRow()

if not arcpy.Exists(gsimpsub2join):
    print('Join gsim stations to HydroBASINS level 5')
    arcpy.SpatialJoin_analysis(target_features=gsimpsub2,
                               join_features=basinatlasl05,
                               out_feature_class=gsimpsub2join,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type="KEEP_COMMON",
                               match_option='INTERSECT'
                               )

#Remove all stations that are not intermittent or in basins level 05 (give average area) where we did not yet have any observations
# Check basins where we already have GRDC stations with a record that spans at least 10 years
# Get number of years of data in GRDC
if not arcpy.Exists(basin5grdcpjoin):
    grdcdistablist = getfilelist(dir=grdc_disdatdir)
    print('Compute number of years of data in GRDC for {} station...'.format(len(grdcdistablist)))
    grdc_u20missy = {}
    x=0
    for file in grdcdistablist:
        if (x % 100)==0 :
            print('Processed {} stations...'.format(x))
        disdf = pd.read_csv(file, sep=';', parse_dates=[0], skiprows=40)
        disdf['year'] = disdf['YYYY-MM-DD'].dt.strftime('%Y')
        disdf_ydatdays = disdf[(disdf[' Original'] != -999) & (disdf[' Original'] != None)].groupby('year')[' Original'].count()
        grdc_u20missy[os.path.splitext(os.path.split(file)[1])[0]] = sum(disdf_ydatdays >= 345)
        x+=1

    arcpy.CopyFeatures_management(grdcp, grdcp_10ymin)
    with arcpy.da.UpdateCursor(grdcp_10ymin, 'GRDC_NO') as cursor:
        for row in cursor:
            if grdc_u20missy[row[0]] < 10:
                cursor.deleteRow()

    arcpy.SpatialJoin_analysis(target_features=basinatlasl05,
                               join_features=grdcp_10ymin,
                               out_feature_class=basin5grdcpjoin,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type="KEEP_COMMON",
                               match_option='INTERSECT')

basgrdcpdict = {row[0] : [row[1], row[2]] for row in
                arcpy.da.SearchCursor(basin5grdcpjoin, ['HYBAS_ID', 'Join_Count', 'SUB_AREA'])} #Get basins where there are GRDC stations
if not arcpy.Exists(gsimpsub3):
    arcpy.CopyFeatures_management(gsimpsub2join, gsimpsub3)
    with arcpy.da.UpdateCursor(gsimpsub3, ['ir_moprop', 'HYBAS_ID']) as cursor:
        for row in cursor:
            if (row[0] < 1/12.0) and (row[1] in basgrdcpdict.keys()):
                cursor.deleteRow()

if not arcpy.Exists(gsimpsub3b):
    arcpy.CopyFeatures_management(gsimpsub2join, gsimpsub3b)
    with arcpy.da.UpdateCursor(gsimpsub3b, ['ir_moprop', 'HYBAS_ID', 'area_correct']) as cursor:
        for row in cursor:
            if ((row[0] >= 1/12.0) or (row[1] not in basgrdcpdict.keys())) or (row[2] >= 100):
                cursor.deleteRow()

#Remove extraneous fields
for f1 in arcpy.ListFields(gsimpsub3):
    if f1.name not in [f2.name for f2 in arcpy.ListFields(gsimpsub)]+['ir_moprop']:
        arcpy.DeleteField_management(gsimpsub3, f1.name)
        arcpy.DeleteField_management(gsimpsub3b, f1.name)

if not arcpy.Exists(gsimpsub3join):
    print('Join gsim stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(gsimpsub3,
                               riveratlas, gsimpsub3join, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.1,
                               distance_field_name='station_river_distance')

    arcpy.AddField_management(gsimpsub3join, field_name='DApercdiff', field_type='FLOAT')
    arcpy.CalculateField_management(gsimpsub3join, field='DApercdiff',
                                    expression='(!area_correct!-!UPLAND_SKM!)/!UPLAND_SKM!',
                                    expression_type='PYTHON')

if not arcpy.Exists(gsimpsub3joinb):
    print('Join gsim stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(gsimpsub3b,
                               riveratlas, gsimpsub3joinb, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.1,
                               distance_field_name='station_river_distance')

    arcpy.AddField_management(gsimpsub3joinb, field_name='DApercdiff', field_type='FLOAT')
    arcpy.CalculateField_management(gsimpsub3joinb, field='DApercdiff',
                                    expression='(!area_correct!-!UPLAND_SKM!)/!UPLAND_SKM!',
                                    expression_type='PYTHON')

#Snap to river network when within 200 m
arcpy.CopyFeatures_management(gsimpsub3join, gsimpsnap)
snapenv = [[riveratlas, 'EDGE', '200 meters']]
arcpy.Snap_edit(gsimpsnap, snapenv)

arcpy.CopyFeatures_management(gsimpsub3joinb, gsimpsnapb)
snapenv = [[riveratlas, 'EDGE', '200 meters']]
arcpy.Snap_edit(gsimpsnapb, snapenv)

#Flag those to check and edit them: beyond 200 m from HydroSHEDS reach OR DApercdiff > 5%
arcpy.CopyFeatures_management(gsimpsnap, gsimpsnapedit)

arcpy.AddField_management(gsimpsnapedit, 'manualsnap_mathis', 'SHORT')
arcpy.AddField_management(gsimpsnapedit, 'snap_commentmathis', 'TEXT')

with arcpy.da.UpdateCursor(gsimpsnapedit, ['DApercdiff', 'station_river_distance', 'manualsnap_mathis']) as cursor:
    for row in cursor:
        if (abs(row[0]) > 0.05) or (row[1] > 200):
            row[2] = 2
            cursor.updateRow(row)

arcpy.CopyFeatures_management(gsimpsnapb, gsimpsnapeditb)

arcpy.AddField_management(gsimpsnapeditb, 'manualsnap_mathis', 'SHORT')
arcpy.AddField_management(gsimpsnapeditb, 'snap_commentmathis', 'TEXT')

with arcpy.da.UpdateCursor(gsimpsnapeditb, ['DApercdiff', 'station_river_distance', 'manualsnap_mathis']) as cursor:
    for row in cursor:
        if (abs(row[0]) > 0.05) or (row[1] > 200):
            row[2] = 2
            cursor.updateRow(row)

#Process to snapping, manual editing and deleting

#Merge the two edited layers
arcpy.Merge_management([gsimpsnapedit, gsimpsnapeditb], gsimpsnapedit_merge)


#Delete those that were marked with -1
if not arcpy.Exists(gsimpsnapclean):
    arcpy.CopyFeatures_management(gsimpsnapedit_merge, gsimpsnapclean)
    with arcpy.da.UpdateCursor(gsimpsnapclean, ['manualsnap_mathis']) as cursor:
        for row in cursor:
            if row[0] == -1:
                cursor.deleteRow()

    #Re-snap points
    snapenv2 = [[riveratlas, 'EDGE', '10 meters']]
    arcpy.Snap_edit(gsimpsnapclean, snapenv2)

    #Delete unneeded columns
    keepcols_gsim = [f.name for f in arcpy.ListFields(gsimpsub3)] + ['manualsnap_mathis', 'DApercdiff', 'snap_commentmathis']
    for f in arcpy.ListFields(gsimpsnapclean):
        if f.name not in keepcols_gsim:
            print('Deleting {}'.format(f.name))
            arcpy.DeleteField_management(gsimpsnapclean, f.name)

#Join stations to nearest river reach in RiverAtlas
if not arcpy.Exists(gsimpcleanjoin):
    print('Join gsim stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(gsimpsnapclean, riveratlas, gsimpcleanjoin,
                               join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                               distance_field_name='station_river_distance')

#Extract GLAD sea mask value to check
ExtractMultiValuesToPoints(in_point_features=gsimpcleanjoin, in_rasters=seamaskbuf3k,
                           bilinear_interpolate_values='NONE')

#Get coordinates
arcpy.AddGeometryAttributes_management(gsimpcleanjoin, Geometry_Properties='POINT_X_Y_Z_M')


###########################################################################################
'''Not used in the end - formatting for error assessment and kriging of residuals
hydromask = os.path.join(datdir, 'HydroATLAS', 'Masks_20200525','hydrosheds_landmask_15s.gdb', 'hys_land_15s')
grdcbuf = os.path.join(outgdb, 'grdcstations_buf50k')
gaugebufdiss = os.path.join(outgdb, grdcbuf + 'diss')

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
'''