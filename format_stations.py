#Purpose: Format GRDC station data and associated environmental variables extract from HydroATLAS

import arcpy
import os

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

#Output variables
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
if not arcpy.Exists(outgdb):
    arcpy.CreateFileGDB_management(os.path.split(outgdb)[0], os.path.split(outgdb)[1])
GRDCp = os.path.join(outgdb, 'GRDCstations')
GRDCpjoin = os.path.join(outgdb, 'GRDCstations_riverjoin')
riveratlas_csv = os.path.join(resdir, 'RiverATLAS_v10tab.csv')

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