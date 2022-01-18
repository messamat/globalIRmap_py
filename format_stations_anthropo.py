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

riveratlas = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')

outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
grdcp = os.path.join(outgdb, 'grdcstations')
grdcpjoinedit = os.path.join(outgdb, 'grdcstations_riverjoinedit')
grdcpjoinedit_anthropo = os.path.join(outgdb, 'grdcstations_riverjoinedit_anthropo')
grdcpclean_anthropo = os.path.join(outgdb, 'grdcstations_clean_anthropo')
grdcpcleanjoin_anthropo = os.path.join(outgdb, 'grdcstations_cleanjoin_anthropo')

gsimresdir = os.path.join(resdir, 'GSIM')
gsimresgdb = os.path.join(gsimresdir, 'GSIM.gdb')
gsimpsub3= os.path.join(gsimresgdb, 'GSIMSstations_sub3')
gsimpsnapedit = os.path.join(gsimresgdb, 'GSIMSstations_riversnapedit')
gsimpsnapeditb = os.path.join(gsimresgdb, 'GSIMSstations_riversnapeditb')
gsimpsnapedit_anthropo = os.path.join(gsimresgdb, 'GSIMSstations_riversnapedit_anthropo')
gsimpsnapeditb_anthropo = os.path.join(gsimresgdb, 'GSIMSstations_riversnapeditb_anthropo')
gsimpsnapedit_anthropo_merge = os.path.join(gsimresgdb, 'GSIMSstations_riversnapedit_anthropo_merge')
gsimpsnapclean_anthropo = os.path.join(gsimresgdb, 'GSIMstations_riversnapclean_anthropo')
gsimpcleanjoin_anthropo = os.path.join(gsimresgdb, 'GSIMstations_cleanriverjoin_anthropo')

######################################## FORMAT GRDC STATIONS ##########################################################
if not arcpy.Exists(grdcpjoinedit_anthropo):
    arcpy.CopyFeatures_management(grdcpjoinedit, grdcpjoinedit_anthropo)

############ Manual edits to re-introduce those  that were included because downstream of reservoirs or diversions

#Delete those that were marked with -1
if not arcpy.Exists(grdcpclean_anthropo):
    arcpy.CopyFeatures_management(grdcpjoinedit_anthropo, grdcpclean_anthropo)
    with arcpy.da.UpdateCursor(grdcpclean_anthropo, ['manualsnap_mathis']) as cursor:
        for row in cursor:
            if row[0] == -1:
                cursor.deleteRow()

    #Delete useless fields and re-spatial join
    for f1 in arcpy.ListFields(grdcpclean_anthropo):
        if f1.name not in [f2.name for f2 in arcpy.ListFields(grdcp)]+['manualsnap_mathis', 'snap_comment_mathis']:
            arcpy.DeleteField_management(grdcpclean_anthropo, f1.name)

arcpy.SpatialJoin_analysis(grdcpclean_anthropo, riveratlas, grdcpcleanjoin_anthropo, join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                           match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                           distance_field_name='station_river_distance')

#Extract GLAD sea mask value to check for tidal reversal as cause of intermittency
# ExtractMultiValuesToPoints(in_point_features=grdcpcleanjoin_anthropo, in_rasters=seamaskbuf3k,
#                            bilinear_interpolate_values='NONE')

#Add coordinates
arcpy.AddGeometryAttributes_management(grdcpcleanjoin_anthropo, Geometry_Properties='POINT_X_Y_Z_M')

######################################## FORMAT GSIM STATIONS ##########################################################
if not arcpy.Exists(gsimpsnapedit_anthropo):
    arcpy.CopyFeatures_management(gsimpsnapedit, gsimpsnapedit_anthropo)

if not arcpy.Exists(gsimpsnapeditb_anthropo):
    arcpy.CopyFeatures_management(gsimpsnapeditb, gsimpsnapeditb_anthropo)

############ Manual edits to re-introduce those  that were included because downstream of reservoirs or diversions

arcpy.Merge_management([gsimpsnapedit_anthropo, gsimpsnapeditb_anthropo], gsimpsnapedit_anthropo_merge)

#Delete those that were marked with -1
if not arcpy.Exists(gsimpsnapclean_anthropo):
    arcpy.CopyFeatures_management(gsimpsnapedit_anthropo_merge, gsimpsnapclean_anthropo)
    with arcpy.da.UpdateCursor(gsimpsnapclean_anthropo, ['manualsnap_mathis']) as cursor:
        for row in cursor:
            if row[0] == -1:
                cursor.deleteRow()

    #Re-snap points
    snapenv2 = [[riveratlas, 'EDGE', '10 meters']]
    arcpy.Snap_edit(gsimpsnapclean_anthropo, snapenv2)

    #Delete unneeded columns
    keepcols_gsim = [f.name for f in arcpy.ListFields(gsimpsub3)] + ['manualsnap_mathis', 'DApercdiff', 'snap_commentmathis']
    for f in arcpy.ListFields(gsimpsnapclean_anthropo):
        if f.name not in keepcols_gsim:
            print('Deleting {}'.format(f.name))
            arcpy.DeleteField_management(gsimpsnapclean_anthropo, f.name)

#Join stations to nearest river reach in RiverAtlas
if not arcpy.Exists(gsimpcleanjoin_anthropo):
    print('Join gsim stations to nearest river reach in RiverAtlas')
    arcpy.SpatialJoin_analysis(gsimpsnapclean_anthropo, riveratlas, gsimpcleanjoin_anthropo,
                               join_operation='JOIN_ONE_TO_ONE', join_type="KEEP_COMMON",
                               match_option='CLOSEST_GEODESIC', search_radius=0.0005,
                               distance_field_name='station_river_distance')

#Extract GLAD sea mask value to check
# ExtractMultiValuesToPoints(in_point_features=gsimpcleanjoin_anthropo, in_rasters=seamaskbuf3k,
#                            bilinear_interpolate_values='NONE')

#Get coordinates
arcpy.AddGeometryAttributes_management(gsimpcleanjoin_anthropo, Geometry_Properties='POINT_X_Y_Z_M')


