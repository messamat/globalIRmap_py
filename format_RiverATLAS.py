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
riveratlasv11 = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v11.gdb', 'RiverATLAS_v11')
hydrolakes = os.path.join(datdir, 'hydrolakes', 'HydroLAKES_polys_v10.gdb', 'HydroLAKES_polys_v10')
basinatlasl03 = os.path.join(datdir, 'HydroATLAS', 'BasinATLAS_v10.gdb', 'BasinATLAS_v10_lev03')
basinatlasl05 = os.path.join(datdir, 'HydroATLAS', 'BasinATLAS_v10.gdb', 'BasinATLAS_v10_lev05')
hydromask = os.path.join(datdir, 'HydroATLAS', 'Masks_20200525','hydrosheds_landmask_15s.gdb', 'hys_land_15s')
seamaskbuf3k = os.path.join(datdir, 'GLAD', 'class99_19_rsp9_buf3k1.tif')

#Output variables
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
if not arcpy.Exists(outgdb):
    arcpy.CreateFileGDB_management(os.path.split(outgdb)[0], os.path.split(outgdb)[1])
riveratlas_csv = os.path.join(resdir, 'RiverATLAS_v10tab.csv')
riveratlasv11_csv = os.path.join(resdir, 'RiverATLAS_v11tab.csv')
riveratlas_b03 = os.path.join(outgdb, 'RiverATLASbas3join')
riveratlas_b05 = os.path.join(outgdb, 'RiverATLASbas5join')

# ---------- Intersect river reaches with lakes ------
rivlakeinters = os.path.join(outgdb, 'riveratlas_hydrolakes_inters')
if not arcpy.Exists(rivlakeinters):
    arcpy.Intersect_analysis(in_features = [riveratlas, hydrolakes], out_feature_class=rivlakeinters,
                             join_attributes= 'ONLY_FID')
    arcpy.AddGeometryAttributes_management(rivlakeinters, Geometry_Properties='LENGTH_GEODESIC', Length_Unit='kilometers')

if not 'INLAKEPERC' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.JoinField_management(in_data=riveratlas, in_field=arcpy.Describe(riveratlas).OIDFieldName,
                               join_table=rivlakeinters, join_field='FID_RiverATLAS_v10', fields='LENGTH_GEO')
    arcpy.AlterField_management(riveratlas, 'LENGTH_GEO', new_field_name='LENGTH_lake',
                                new_field_alias='LENGTH_lake')
    arcpy.AddField_management(riveratlas, 'INLAKEPERC', field_type='float')
    with arcpy.da.UpdateCursor(riveratlas, ['LENGTH_KM', 'LENGTH_lake', 'INLAKEPERC']) as cursor:
        for row in cursor:
            if row[1] is not None:
                row[2] = row[1]/row[0]
            else:
                row[2] = 0
            cursor.updateRow(row)

# ---------- Flag reaches with mean annual discharge == 0 (dis_m3_pyr) AND
# ---------- (ORD_STRA == 1 OR downstream of another reach with discharge == 0)
if not 'NOFLOW' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'NOFLOW', 'SHORT')

nextdownset = set()
x = 0
with arcpy.da.UpdateCursor(riveratlas, ['ORD_STRA', 'dis_m3_pyr', 'NOFLOW', 'NEXT_DOWN']) as cursor:
    for row in cursor:
        if x % 100000 == 0:
            print(x)
        if row[0] == 1:
            if row[1] == 0: #If Strahler order ==1 and mean annual natural discharge == 0
                row[2] = 1 #Set NOFLOW to 0
                nextdownset.add(row[3]) #Add reach ID to set
            else:
                row[2] = 0
        cursor.updateRow(row)
        x += 1

while len(nextdownset) > 0:
    print('Processing {0} reaches to determine inclusion status in analysis'.format(len(nextdownset)))
    with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'dis_m3_pyr', 'NOFLOW', 'NEXT_DOWN'],
                               where_clause='NOFLOW IS NULL') as cursor:
        for row in cursor:
            if (row[0] in nextdownset):
                if (row[1] == 0):  # If downstream of a noflow reach and mean annual natural discharge == 0
                    row[2] = 1  # Set NOFLOW to 0
                    nextdownset.add(row[3])  # Add reach ID to set
                else:
                    row[2] = 0
                nextdownset.remove(row[0]) #Remove HYRIV_ID from nextdownset
                cursor.updateRow(row)

# ---------- Associate reaches with HydroBASIN level 05 ------
arcpy.SpatialJoin_analysis(target_features=riveratlas, join_features=basinatlasl05,
                           out_feature_class=riveratlas_b05, join_operation="JOIN_ONE_TO_ONE",
                           join_type = 'KEEP_ALL', match_option="HAVE_THEIR_CENTER_IN")
basdict = {row[0]:row[1] for row in arcpy.da.SearchCursor(riveratlas_b05, ['HYRIV_ID', 'PFAF_ID'])}
arcpy.AddField_management(riveratlas, field_name='PFAF_ID05', field_type='LONG')
with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'PFAF_ID05']) as cursor:
    for row in cursor:
        row[1] = basdict[row[0]]
        cursor.updateRow(row)

# ---------- Compute coordinates ------
arcpy.AddGeometryAttributes_management(riveratlas, 'LINE_START_MID_END')

# ---------- Export attribute table of RiverATLAS with selected ------------------
if not arcpy.Exists(riveratlas_csv):
    print('Exporting CSV table of RiverATLAS v1.0 attributes')
    arcpy.CopyRows_management(in_rows = riveratlas, out_table=riveratlas_csv)
if not arcpy.Exists(riveratlasv11_csv):
    print('Exporting CSV table of RiverATLAS v1.1 attributes')
    arcpy.CopyRows_management(in_rows=riveratlasv11, out_table=riveratlasv11_csv)
