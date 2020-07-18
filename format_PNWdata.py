from setup_localIRformatting import *

#PNW
datdir_pnw = os.path.join(insitudatdir, 'PNW')
pathcheckcreate(datdir_pnw)
resdir_pnw = os.path.join(insituresdir, 'pnw.gdb')
pathcheckcreate(resdir_pnw)

facraw_pnw = os.path.join(datdir_pnw, 'fac_taudem_17all_int.tif')
obsraw_pnw = os.path.join(datdir_pnw, 'StreamflowPermObs.shp')
obsub_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_sub')
obsnodupli_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_nodupli')
sparse_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_sparse')
sparseattri_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_sparseattri')
obsjoin_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_nhdjoin')
obsjoin_atlas = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlas')
obsjoin_atlassub = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlasjoinsub')
obsjoin_atlasnap = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlassnap')
obsjoin_atlasnapedit = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlassnapedit')
obsjoin_atlasnapclean = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlassnapclean')
obsjoin_atlasnapcleanedit = os.path.join(resdir_pnw, 'StreamflowPermObs_RiverAtlassnapclean_edit')
obsfinal_pnw = os.path.join(resdir_pnw, 'StreamflowPermObs_final')

nhdpnw = os.path.join(resdir_pnw, 'NHDpnw')
nhdvaapnw = os.path.join(resdir_pnw, 'NHDvaapnw')


#Only keep points for which use is:
arcpy.MakeFeatureLayer_management(obsraw_pnw, out_layer='pnwlyr',
                                  where_clause='(Use IN {}) AND (Month > 6)'.format(str(("No; Outside 2004-2016", "Yes"))))
if not arcpy.Exists(obsub_pnw):
    arcpy.CopyFeatures_management('pnwlyr', obsub_pnw)

#Flag groups of duplicates in obsub_pnw (utility_functions)
group_duplishape(in_features=obsub_pnw, deletedupli=True, out_featuresnodupli=obsnodupli_pnw)

#Remove all points within 1 km of each other
#if not arcpy.Exists(sparseattri_pnw):
#    arcpy.CreateRandomPoints_management(out_path=os.path.split(sparse_pnw)[0],
#                                        out_name=os.path.split(sparse_pnw)[1],
#                                        constraining_feature_class=obsnodupli_pnw,
#                                        number_of_points_or_field=nunique,
#                                        minimum_allowed_distance='1000 meters')
#    arcpy.MakeFeatureLayer_management(sparse_pnw, 'sparselyr')
#    arcpy.AddJoin_management('sparselyr', in_field='CID',
#                             join_table=obsnodupli_pnw, join_field=arcpy.Describe(obsnodupli_pnw).OIDFieldName)
#    arcpy.CopyFeatures_management('sparselyr', sparseattri_pnw)

#Get NHDplus for all HUC4 basins extracted from HUC_8 and merge them
huc4list = {row[0][0:4] for row in arcpy.da.SearchCursor(obsnodupli_pnw, 'HUC_8')}
if not arcpy.Exists(nhdpnw):
    nhdlist = [os.path.join(datdir_us, "NHDPLUS_H_{0}_HU4_GDB.gdb".format(huc4), "Hydrography", "NHDFlowLine") for huc4 in huc4list]
    arcpy.Merge_management(nhdlist, nhdpnw)

#Spatial join points to NHDplus high res
if not arcpy.Exists(obsjoin_pnw):
    print('Spatial join, processing {}...'.format(obsjoin_pnw))
    arcpy.SpatialJoin_analysis(obsnodupli_pnw, nhdpnw, out_feature_class=obsjoin_pnw, join_operation="JOIN_ONE_TO_ONE",
                               join_type='KEEP_ALL', match_option='CLOSEST_GEODESIC', distance_field_name='distnhd')

#Get flow accumulation and convert it to km2 (rather than pixels)
ExtractMultiValuesToPoints(obsjoin_pnw, in_rasters=[[facraw_pnw, 'facnhd']])
arcpy.AddField_management(obsjoin_pnw, 'facnhd_km2', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_pnw, ['facnhd', 'facnhd_km2']) as cursor:
    for row in cursor:
        if row[0] is not None:
            row[1] = row[0]*900/float(10**6)
            cursor.updateRow(row)

#Spatial join to RiverATLAS and compute ratio of drainage areas bqsed on the line
obsjoin_pnwproj = "{}_proj".format(obsjoin_pnw)
arcpy.Project_management(obsjoin_pnw, obsjoin_pnwproj, out_coor_system=riveratlas)
if not arcpy.Exists(obsjoin_pnw):
    arcpy.SpatialJoin_analysis(obsjoin_pnwproj, riveratlas,
                               out_feature_class=obsjoin_atlas, join_operation="JOIN_ONE_TO_ONE",
                               join_type='KEEP_ALL', match_option='CLOSEST_GEODESIC', distance_field_name='distatlas')

arcpy.AddField_management(obsjoin_atlas, 'facratio_line', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_atlas, ['facnhd_km2', 'UPLAND_SKM', 'facratio']) as cursor:
    for row in cursor:
        if row[0] is not None:
            row[2] = row[1]/float(row[0])
            cursor.updateRow(row)

#Remove those that are over 500 m away and NHDplus upstream area < 10 km2 OR > 1000 m
arcpy.MakeFeatureLayer_management(obsjoin_atlas, 'atlasjoinlyr',
                                  where_clause='NOT ((((distatlas > 500) OR (facratio_line > 3)) AND (facnhd_km2 < 10)) '
                                               'OR (facnhd_km2 Is Null)) OR (')
arcpy.CopyFeatures_management('atlasjoinlyr', obsjoin_atlassub)

#Snap all those that are within 1000 m from a RiverATLAS segment
arcpy.CopyFeatures_management(obsjoin_atlassub, obsjoin_atlasnap)
snapenv = [[riveratlas, 'EDGE', '1000 meters']]
arcpy.Snap_edit(obsjoin_atlasnap, snapenv)

#Extract flow accumulation directly from HydroSHEDS flow accumulation grid
ExtractMultiValuesToPoints(obsjoin_atlasnap, in_rasters=hydroacc, bilinear_interpolate_values='NONE')

#Compute a second flow accumulation ratio for checking
arcpy.AddField_management(obsjoin_atlasnap, 'facratio_ras', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_atlasnap, ['facratio_ras', 'facnhd_km2', 'up_area_skm_15s']) as cursor:
    for row in cursor:
        if row[2] > 0:
            row[0] = row[1]/float(row[2])
        else:
            row[0] = -9999.0
        cursor.updateRow(row)

#When multiple points are associated with a single HydroATLAS reach, keep the one with the lowest facratio
#duplidict = defaultdict(float)
#with arcpy.da.SearchCursor(obsjoin_atlassub, ['HYRIV_ID', 'facratio_line']) as cursor:
#    for row in cursor:
#        if (row[0] not in duplidict) or (abs(row[1] - 1) < duplidict[row[0]]):
#            duplidict[row[0]] = (abs(row[1] - 1))
#with arcpy.da.UpdateCursor(obsjoin_atlassub, ['HYRIV_ID', 'facratio_line']) as cursor:
#    for row in cursor:
#        if (abs(row[1] - 1)) > duplidict[row[0]]:
#            print('Delete {}'.format(row[0]))
#            cursor.deleteRow()

#Copy for manual editing
arcpy.CopyFeatures_management(obsjoin_atlasnap, obsjoin_atlasnapedit)

#Create fields for manual inspection and editing
arcpy.AddField_management(obsjoin_atlasnapedit, 'manualsnap', 'SHORT')
arcpy.AddField_management(obsjoin_atlasnapedit, 'snap_comment', 'TEXT')

######## INSPECT AND EDIT MANUALLY ##############
#Inspect all gauges > 200 m away from river atlas
#Inspect all gauges with 0.70 < facratio < 1.30
#manualsap: 0 when inspected and not moved, 1 when moved, -1 when to delete
#When possible, delete points that are on the side channel of a bifurcation
#Delete points that are too far upstream from a first order HydroSHEDS reach and that are separated from it by a relatively large tributary confluence
#Delete points that are in areas that are messy in the NHD so that the topology can't be compared to HydroSHEDS

#Copy for deleting and resnapping
arcpy.CopyFeatures_management(obsjoin_atlasnapedit, obsjoin_atlasnapclean)

#Delete points with -1
with arcpy.da.UpdateCursor(obsjoin_atlasnapclean, ['manualsnap']) as cursor:
    for row in cursor:
        if row[0] == -1:
            cursor.deleteRow()

#Re-snap points
snapenv = [[riveratlas, 'EDGE', '200 meters']]
arcpy.Snap_edit(obsjoin_atlasnapclean, snapenv)

#Delete unneeded columns
keepcols = [arcpy.Describe(obsjoin_atlasnapclean).OIDFieldName, 'Shape',
            'OBJECTID', 'OBJECTID_2', 'OBJECTID_3', 'Source', 'Date', 'Category', 'Use', 'Edit', 'Year', 'Month',
            'Permanent_', 'FDate', 'Resolution', 'LengthKM', 'ReachCode', 'FType', 'FCode', 'NHDPlusID',
            'fac_km2', 'manualsnap', 'snap_comment', 'facratio']
for f in arcpy.ListFields(obsjoin_atlasnapclean):
    if f.name not in keepcols:
        print('Deleting {}'.format(f.name))
        arcpy.DeleteField_management(obsjoin_atlasnapclean, f.name)

#Copy for deleting and resnapping
arcpy.CopyFeatures_management(obsjoin_atlasnapclean, obsjoin_atlasnapcleanedit)
arcpy.AddField_management(obsjoin_atlasnapcleanedit, 'manualsnap2', 'SHORT')

#Copy to final layer
arcpy.MakeFeatureLayer_management(obsjoin_atlasnapcleanedit, 'cleanlyr',
                                  where_clause='(NOT manualsnap2 = -1) OR (manualsnap2 IS NULL)')
if not arcpy.Exists(obsfinal_pnw):
    arcpy.SpatialJoin_analysis('cleanlyr', riveratlas,
                               out_feature_class=obsfinal_pnw, join_operation="JOIN_ONE_TO_ONE",
                               join_type='KEEP_ALL', match_option='CLOSEST_GEODESIC', distance_field_name='distatlas')