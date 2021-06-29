from setup_localIRformatting import *

#Onde Eau dirs
datdir = os.path.join(insitudatdir, 'OndeEau')
pathcheckcreate(datdir)
resdir = os.path.join(insituresdir, 'OndeEau.gdb')
pathcheckcreate(resdir)

cartnet_raw =os.path.join(datdir,'TronconHydrograElt_FXX.shp')
cartnet = os.path.join(resdir,'carthage_network')
obsraw= os.path.join(resdir,'obsraw')
obsnodupli = os.path.join(resdir, 'obsnodupli')
ondeau_mergecsv = os.path.join(datdir, 'onde_france_merge.csv')
obscartjoin = os.path.join(resdir, 'obs_cartnet_spatialjoin')
obscartjoinedit = os.path.join(resdir, 'obs_cartnet_spatialjoinedit')
obscartjoinclean = os.path.join(resdir, 'obs_cartnet_spatialjoinclean')
cartobsjoin = os.path.join(resdir, 'cartnet_obsclean_spatialjoin')
cartobsvert = os.path.join(resdir, 'cartnet_obsclean_spatialjoinvert')

atlassub = os.path.join(resdir, 'RiverAtlas_sub')
atlassubproj = os.path.join(resdir, 'RiverAtlas_subproj')
atlassub1000 = os.path.join(resdir, 'RiverAtlascartnet_sub1000m')
atlasvert = os.path.join(resdir, 'RiverAtlascartnet_sub1000m_vertices')
atlascartobs_near = os.path.join(resdir, 'RiverAtlASvertices_cartnet_obsclean_near')

obsatlas_near = os.path.join(resdir, 'obsclean_RiverAtlas_near')

cartobsjoin_edit = os.path.join(resdir, 'cartnet_obsclean_spatialjoinedit')
cartobsjoin_editFINAL = os.path.join(resdir, 'cartnet_obsclean_spatialjoinedit20200719')
cartobsjoin_clean = os.path.join(resdir, 'cartnet_obsclean_spatialjoinclean')

obsfinal = os.path.join(resdir, 'obs_final')
atlassub1000route = os.path.join(resdir, 'RiverAtlascartnet_sub1000mroute')
obsfinal_wgs = os.path.join(resdir, 'obs_finalwgs')

################################################### ANALYSIS ###########################################################
# <LbSiteHydro>Name of site
# <CdSiteHydro> Unique code for site
# <Annee>Year of observation
# <TypeCampObservations>:
#       Usuelle: La campagne usuelle (réseau ONDE) vise à acquérir de la donnée pour la constitution d'un réseau
#       stable de connaissance. Elle est commune à l'ensemble des départements, sa fréquence d'observation est mensuelle,
#       au plus près du 25 de chaque mois à plus ou moins 2 jours, sur la période de mai à septembre.
#   `   Crise:La campagne de crise (réseau ONDE) vise à acquérir de la donnée destinée aux services de l'État en charge
#       de la gestion dela crise en période de sécheresse.
# <DtRealObservation>': Datee of observation
# <LbRsObservationDpt>:Observation result label (Departmental typology)
# <RsObservationDpt>: Observation result(Department typology)
# <LbRsObservationNat> Observation result label (National typology)
# <RsObservationNat> Observation result (National typology)
# <NomEntiteHydrographique>
# <CdTronconHydrographique>
# <CoordXSiteHydro>
# <CoordYSiteHydro>
# <ProjCoordSiteHydro>
# FLG: End of line--- not useful - NOT FLAG

# Copy Carthage network to gdb
arcpy.CopyFeatures_management(cartnet_raw, cartnet)

# Merge OndeEau csv
if not arcpy.Exists(ondeau_mergecsv):
    mergedelcsv(dir=datdir,
                repattern="onde_france_[0-9]{4}[.]csv",
                outfile=ondeau_mergecsv,
                returndf=False,
                delete=False,
                verbose=True)

# Convert to points (SpatialReference from ProjCoordSiteHydro referenced to http://mdm.sandre.eaufrance.fr/node/297134 - code 26 (247 observations also have code 2154 which must be a mistake and refers to the ESPG code)
if not arcpy.Exists(obsraw):
    arcpy.MakeXYEventLayer_management(table=ondeau_mergecsv,
                                      in_x_field="<CoordXSiteHydro>",
                                      in_y_field="<CoordYSiteHydro>",
                                      out_layer='ondelyr',
                                      spatial_reference=arcpy.SpatialReference(
                                          2154))  # Lambert 93 for mainland France and Corsica
    arcpy.CopyFeatures_management('ondelyr', obsraw)

#Delete duplicated locations
if not arcpy.Exists(obsnodupli):
    group_duplishape(in_features=obsraw, deletedupli=True, out_featuresnodupli=obsnodupli)

#Join to Carthage network by Spatial Join
if not arcpy.Exists(obscartjoin):
    arcpy.SpatialJoin_analysis(target_features=obsnodupli,
                               join_features=cartnet,
                               out_feature_class=obscartjoin,
                               join_operation='JOIN_ONE_TO_ONE',
                               join_type='KEEP_ALL',
                               match_option='CLOSEST',
                               distance_field_name= 'distcart')

#Check all observations > 10 m from a reach (Troncon) and all those whose non-null 'CdTronconHydrographique' does not match 'CdTronconH'
#If same name and same ID, and within 25 m. Either move to reach if > 10 m or mark 0.
#If Onde Eau has NULL for ID but within 5 m. Mark as 0. > 5 m, check satellite imagery. If not sign of other stream, mark 0. If > 10 m, move to match.
#If different name and/or ID, always check the ones around. If > 10 m from stream in particular, consider moving if better match.
#Always check for close confluence
#Even if there is an ID, always prioritize site name match to river name
#

#NULL - not checked
#-1 - checked, delete location. No corresponding river reach (or too upstream by > 100 m)
#0 - checked, location not adjusted
#1 - checked, location adjusted
#2 - checked, no corresponding river reach in Carthage, but river reach in Pella et al. 2012 theoretic river network
arcpy.CopyFeatures_management(obscartjoin, obscartjoinedit)
arcpy.AddField_management(obscartjoinedit, 'checkjoincart', 'SHORT')
arcpy.AddField_management(obscartjoinedit, 'manualsnapcart', 'SHORT')

with arcpy.da.UpdateCursor(obscartjoinedit,
                           ['distcart', 'F_CdTronconHydrographique_', 'CdTronconH', 'checkjoincart'])  as cursor:
    for row in cursor:
        row[3] = 0

        if row[0]>10: #If distance > 10 m
            row[3] = 1
        elif (row[2] is not None) and (row[2] not in ['', ' ']):
            if row[1] is not None:
                if not re.search(re.compile(re.sub('-', '.', row[1])), row[2]): #OndeEau Codes have erroneous hyphens so use fuzzy matching where there are hyphens
                    row[3] = 1
            else:
                row[3] = 1

        cursor.updateRow(row)

[f.name for f in arcpy.ListFields(obscartjoin)]

#Remove points to be deleted
arcpy.CopyFeatures_management(obscartjoinedit, obscartjoinclean)
with arcpy.da.UpdateCursor(obscartjoinclean, ['manualsnapcart'])  as cursor:
    for row in cursor:
        if row[0] == 1:
            cursor.deleteRow()

#Snap points to Carthage network (max dist 10 m)
arcpy.Snap_edit(in_features=obscartjoinclean, snap_environment=[[cartnet, 'EDGE', '10 meters']])

#Delete additional fields from points and spatially join Carthage network to points
for f in arcpy.ListFields(obscartjoinclean):
    if f not in [f.name for f in arcpy.ListFields(obsraw)]+['checkjoincart', 'manualsnapcart', 'distcart']:
        arcpy.DeleteField_management(obscartjoinclean, f.name)

arcpy.SpatialJoin_analysis(target_features=cartnet,
                           join_features=obscartjoinclean,
                           out_feature_class=cartobsjoin,
                           join_operation="JOIN_ONE_TO_MANY",
                           join_type="KEEP_COMMON",
                           match_option="INTERSECT",
                           search_radius="0.1 meters")

#Convert Carthage segments to vertices (max 100 m apart)
arcpy.Densify_edit(in_features=cartobsjoin, densification_method='DISTANCE',
                   distance='100 meters')
arcpy.FeatureVerticesToPoints_management(in_features=cartobsjoin,
                                         out_feature_class=cartobsvert,
                                         point_location='ALL')

########## JOIN CARTHAGE NETWORK WITH OBSERVATIONS TO RIVERATLAS ##############
#Select RiverATLAS segments within one kilometer from any Carthage reach of interest
#Re-project Carthage observations to WGS84
arcpy.CopyFeatures_management(arcpy.Describe(obscartjoinclean).extent.polygon,
                              os.path.join(resdir, 'obspoly'))
arcpy.Project_management(os.path.join(resdir, 'obspoly'),
                         os.path.join(resdir, 'obspolyproj'),
                         riveratlas)

#Select all RiverATLAS reaches within  50 kilometers from extent
arcpy.MakeFeatureLayer_management(riveratlas, 'atlaslyr')
arcpy.SelectLayerByLocation_management(in_layer='atlaslyr', #First select all RiverAtlas reaches intersecting the extent of the observations
                                       overlap_type="INTERSECT",
                                       select_features=os.path.join(resdir, 'obspolyproj'),
                                       selection_type='NEW_SELECTION',
                                       search_distance='50 kilometers')
arcpy.CopyFeatures_management('atlaslyr', atlassub)

#Project to Carthage projection
arcpy.Project_management(atlassub, atlassubproj, out_coor_system=cartnet)
arcpy.Delete_management('atlaslyr')

#select all reaches within 1000 meters from a Carthage-Onde Eau line observation
arcpy.MakeFeatureLayer_management(atlassubproj, 'atlaslyr')
arcpy.SelectLayerByLocation_management(in_layer='atlaslyr',
                                       overlap_type="WITHIN_A_DISTANCE_GEODESIC",
                                       select_features=cartobsjoin,
                                       search_distance='1000 meters',
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('atlaslyr', atlassub1000)

#Convert HydroSHEDS segments to vertices
arcpy.Densify_edit(in_features=atlassub1000, densification_method='DISTANCE',
                   distance='100 meters')
arcpy.FeatureVerticesToPoints_management(in_features=atlassub1000,
                                         out_feature_class=atlasvert,
                                         point_location='ALL')
vertdict = {row[0]: row[1] for row in arcpy.da.SearchCursor(atlasvert, ['OID@', 'HYRIV_ID'])}
vertndict = defaultdict(int)
for k,v in vertdict.iteritems():
    vertndict[v] +=1

#Associate each Carthage-Onde observation line with all RiverATLAS lines within 1000 m

#Compute distance from Onde point observations to RiverATLAS lines
arcpy.GenerateNearTable_analysis(in_features=obscartjoinclean,
                                 near_features=atlassub1000,
                                 out_table= obsatlas_near,
                                 search_radius='1000 meters',
                                 location='location',
                                 closest='ALL',
                                 method='PLANAR')

#Compute distance from RiverATLAS vertices to nearest point on paired Cathage-Onde line
arcpy.GenerateNearTable_analysis(in_features=atlasvert,
                                 near_features=cartobsjoin,
                                 out_table=atlascartobs_near,
                                 search_radius='1000 meters',
                                 location='location',
                                 closest = 'ALL',
                                 method = 'PLANAR')

#Create dictionary with
#{Hyriv_ID : {Carthage_ID : sum(distance of HydroSHEDS vertices to Carthage line), # of vertices < 5000 m from Carthage line} }
neardict = {}
with arcpy.da.SearchCursor(atlascartobs_near, ['IN_FID', 'NEAR_FID', 'NEAR_DIST']) as cursor:
    for row in cursor:
        hyriv_id = vertdict[row[0]] #ID of in_layer's vertice's line (e.g. HydroSHEDS)
        if hyriv_id not in neardict:
            neardict[hyriv_id] = defaultdict(dict)

        if row[1] not in neardict[hyriv_id]: #If join line (e.g. Carthage line) not already associated with the in_layer (e.g. HydroSHEDS) line
            neardict[hyriv_id][row[1]][0] = row[2] #Distance between (e.g. HydroSHED) vertex to (e.g. Carthage) join line
            neardict[hyriv_id][row[1]][1] = 1 #Number of vertices within 1000 m of line
        else:
            neardict[hyriv_id][row[1]][0] += row[2]
            neardict[hyriv_id][row[1]][1] += 1

#Convert to dataframe and write out if needs manual inspection later
cartoid = arcpy.Describe(cartobsjoin).OIDFieldName
neardf = pd.concat({
    k: pd.DataFrame.from_dict(v, 'index') for k, v in neardict.items()},
    axis=0).reset_index()
neardf.columns = ['HYRIV_ID', cartoid, 'distsum', 'nvertu1k']
neardf['CARTtoATLAS_dist'] = neardf['distsum']/neardf['nvertu1k'].astype(float)

### Join distance from points to HydroSHEDS
# def getiddf(in_feature, IDfield, valuefields):
#     iddict = {row[0]: row[1] for row in arcpy.da.SearchCursor(in_feature, [IDfield, valuefields])}
#     iddf = pd.DataFrame.from_dict(iddict, 'index').reset_index()
#     iddf.columns = [arcpy.Describe(in_feature).OIDFieldName, IDfield]
#     return(iddf)

obsiddf = get_arctab_df(in_tab=obscartjoinclean, IDfield='OID@', valuefields='F_CdSiteHydro_', resetid=True)
atlasiddf = get_arctab_df(in_tab=atlassub1000, IDfield='OID@', valuefields='HYRIV_ID', resetid=True)
obsatlasneardf = get_arctab_df(in_tab=obsatlas_near, IDfield='OID@', valuefields=['IN_FID', 'NEAR_FID', 'NEAR_DIST'])

nearptdf = obsatlasneardf.merge(obsiddf, left_on='IN_FID', right_on=arcpy.Describe(obscartjoinclean).OIDFieldName).\
    merge(atlasiddf, left_on='NEAR_FID', right_on=arcpy.Describe(atlassub1000).OIDFieldName)
nearptdf = nearptdf.drop([f for f in ['OBJECTID', 'OBJECTID_x', 'OBJECTID_y', 'IN_FID', 'NEAR_FID']
                          if f in nearptdf.columns], axis=1)
nearptdf.columns = ['OBStoATLAS_dist'] + list(nearptdf.columns[1:])

### Compute bearing
hydrobear = pd.DataFrame.from_dict(
    linevert_wbearing(in_vert=atlasvert, in_IDfield='HYRIV_ID'),
    'index').reset_index()
hydrobear.columns = ['HYRIV_ID', 'hydrobearing']

cartbear = pd.DataFrame.from_dict(
    linevert_wbearing(in_vert = cartobsvert, in_IDfield = 'IdTronconH'),
    'index').reset_index()
cartbear.columns = ['IdTronconH', 'cartbearing']

cart_objectidtroncon_df = get_arctab_df(in_tab=cartobsjoin, IDfield='OID@',
                                        valuefields=['IdTronconH', 'F_CdSiteHydro_'], resetid=True)

#Join distance df with bearing dfs (getting rid of Carthage-Onde eau observations that do not have any RiverATLAS lines within 1000 meters)
nearbeardf = cart_objectidtroncon_df.merge(neardf, how='left', on=cartoid). \
    dropna(axis='index', subset=['HYRIV_ID']). \
    merge(hydrobear, on='HYRIV_ID'). \
    merge(cartbear, on='IdTronconH')

#Join df with point near distance
nearbeardfmerge = nearbeardf.merge(nearptdf, on=['F_CdSiteHydro_', 'HYRIV_ID'])

mindistdf = nearbeardfmerge.groupby("IdTronconH").agg({'CARTtoATLAS_dist':'min', 'OBStoATLAS_dist':'min'})
nearbeardfmerge = nearbeardfmerge.merge(mindistdf, on= "IdTronconH", suffixes=('', 'min'))

maxvertdf = pd.DataFrame(nearbeardfmerge.groupby("IdTronconH")['nvertu1k'].max())
nearbeardfmerge = nearbeardfmerge.merge(maxvertdf, on="IdTronconH", suffixes=('', 'max'))

#Compute absolute difference in bearing between HydroSHEDS and Carthage lines
nearbeardfmerge['bearingdiff'] = abs(nearbeardfmerge['hydrobearing'] - nearbeardfmerge['cartbearing'])
nearbeardfmerge['CARTtoATLAS_distindex'] = nearbeardfmerge['CARTtoATLAS_distmin']/nearbeardfmerge['CARTtoATLAS_dist']
nearbeardfmerge['OBStoATLAS_distindex'] = nearbeardfmerge['OBStoATLAS_distmin']/nearbeardfmerge['OBStoATLAS_dist']
nearbeardfmerge['nvertindex'] = nearbeardfmerge['nvertu1k']/nearbeardfmerge['nvertu1kmax']

nearbeardfmerge['matchindex'] = (2*nearbeardfmerge['CARTtoATLAS_distindex'] +
                                 2*nearbeardfmerge['OBStoATLAS_distindex'] +
                                 2*(1-(nearbeardfmerge['bearingdiff']/180.0)) +
                                 2*nearbeardfmerge['nvertindex'])/8

#Export to CSV
nearbeardfmerge.sort_values(["IdTronconH", "matchindex"]).to_csv(os.path.join(insituresdir, 'netjoindistbear.csv'))

#Join to Carthage subset
arcpy.AddField_management(cartobsjoin, 'HYRIV_ID', 'LONG')
arcpy.AddField_management(cartobsjoin, 'nvertu1k', 'LONG')
arcpy.AddField_management(cartobsjoin, 'CARTtoATLAS_dist', 'FLOAT')
arcpy.AddField_management(cartobsjoin, 'OBStoATLAS_dist', 'FLOAT')
arcpy.AddField_management(cartobsjoin, 'bearingdiff', 'FLOAT')
arcpy.AddField_management(cartobsjoin, 'matchindex', 'FLOAT')
nearbeardf_cartsub = nearbeardfmerge.sort_values(["IdTronconH", "matchindex"], ascending=False).\
    drop_duplicates(subset='IdTronconH')
with arcpy.da.UpdateCursor(cartobsjoin, ['IdTronconH', 'HYRIV_ID',
                                              'nvertu1k', 'CARTtoATLAS_dist', 'OBStoATLAS_dist',
                                              'bearingdiff', 'matchindex']) as cursor:
    for row in cursor:
        dfrow = nearbeardf_cartsub.loc[nearbeardf_cartsub['IdTronconH'] == row[0]]
        if len(dfrow['HYRIV_ID']) != 0:
            row[1] = dfrow['HYRIV_ID'].values[0]
            row[2] = dfrow['nvertu1k'].values[0]
            row[3] = dfrow['CARTtoATLAS_dist'].values[0]
            row[4] = dfrow['OBStoATLAS_dist'].values[0]
            row[5] = dfrow['bearingdiff'].values[0]
            row[6] = dfrow['matchindex'].values[0]
            cursor.updateRow(row)

#Inspect and edit
arcpy.CopyFeatures_management(cartobsjoin, cartobsjoin_edit)
arcpy.AddField_management(cartobsjoin_edit, 'manualnetmatch', 'TEXT')

#Modification in the future:
#  - compute distance indices based on 1000 m (max dist) rather than based on min.
#  - check bearing function

#Check all points
#Overlay each point, associated Carthage reach and HYDROSHEDS network, highlighting the HydroSHEDS reach associated with
# each observation. Generalled stayed between 1:25,000 and 1:65,000 zoom levels depending on network configuration
#(overlay cartobsjoin_edit with atlassub1000, and the original onde observations: see onde_checknetjoin.mxd)
#Exclude all those that are too unclear, those >= 500 m from upstream limit of corresponding reach or
# that don't have any corresponding segment in HydroSHEDS.
# Add 0 for all those that were checked, -1 for those to delete, and the HYRIV_ID of the correct HydroSHEDS reach when a change was needed

#Inspect and edit
arcpy.CopyFeatures_management(cartobsjoin_editFINAL, cartobsjoin_clean)

#Remove -1, create new field to add edited HYRIV_ID
if not 'HYRIV_IDjoinedit' in [f for f in arcpy.ListFields(cartobsjoin_clean)]:
    print('Create HYRIV_IDjoinedit field...')
    arcpy.AddField_management(cartobsjoin_clean, 'HYRIV_IDjoinedit', 'LONG')

joindict = {}
with arcpy.da.UpdateCursor(cartobsjoin_clean, ['manualnetmatch', 'HYRIV_IDjoinedit', 'HYRIV_ID', 'F_CdSiteHydro_']) \
        as cursor:
    for row in cursor:
        if str(row[0]) == "-1":
            cursor.deleteRow()
            continue
        elif str(row[0]) == "0":
            row[1] = row[2]
        else:
            row[1] = int(row[0])
        cursor.updateRow(row)

        joindict[row[3]] = row[1]

#Add joining field to point observations
arcpy.CopyFeatures_management(obscartjoinclean, obsfinal)
if not 'HYRIV_IDjoinedit' in [f for f in arcpy.ListFields(obsfinal)]:
    print('Create HYRIV_IDjoinedit field...')
    arcpy.AddField_management(obsfinal, 'HYRIV_IDjoinedit', 'LONG')

with arcpy.da.UpdateCursor(obsfinal, ['HYRIV_IDjoinedit', 'F_CdSiteHydro_']) as cursor:
    for row in cursor:
        if row[1] in joindict:
            row[0] = joindict[row[1]]
            cursor.updateRow(row)
        else:
            cursor.deleteRow()
            continue

# Snap to specific feature (create fast function)
indivsnap(in_features = obsfinal,
          in_field = 'HYRIV_IDjoinedit',
          join_field = 'HYRIV_ID',
          indivsnap_env = [[atlassub1000, 'EDGE', '2000 meters']])

#Compute how far down the line the site is as length and percentage of line's length
if not all(v in [f for f in arcpy.ListFields(atlassub1000)] for v in ['fromM', 'toM']):
    arcpy.AddField_management(atlassub1000, 'fromM', 'DOUBLE')
    arcpy.AddField_management(atlassub1000, 'toM', 'DOUBLE')

    with arcpy.da.UpdateCursor(atlassub1000, ['fromM', 'toM', 'SHAPE@LENGTH']) as cursor:
        for row in cursor:
            row[0] = 0
            row[1] = row[2]
            cursor.updateRow(row)

if not arcpy.Exists(atlassub1000route):
    arcpy.CreateRoutes_lr(in_line_features=atlassub1000,
                          route_id_field='HYRIV_ID',
                          out_feature_class=atlassub1000route,
                          measure_source='TWO_FIELDS',
                          from_measure_field='fromM',
                          to_measure_field='toM')

obsnet_locatetab = os.path.join(resdir, 'obsfinal_atlas_locate')
if not arcpy.Exists(obsnet_locatetab):
    arcpy.LocateFeaturesAlongRoutes_lr(in_features=obsfinal,
                                       in_routes=atlassub1000route,
                                       route_id_field='HYRIV_ID',
                                       radius_or_tolerance='0.1',
                                       out_table=obsnet_locatetab,
                                       out_event_properties= "HYRIV_ID POINT fromM toM",
                                       route_locations='FIRST',
                                       in_fields='NO_FIELDS')

if not all(v in [f for f in arcpy.ListFields(obsfinal)] for v in ['HYDROSHEDSdis', 'HYDROSHEDSDA']):
    arcpy.JoinField_management(in_data=obsfinal, in_field=arcpy.Describe(obsfinal).OIDFieldName,
                               join_table=obsnet_locatetab, join_field='INPUTOID', fields='fromM')
    arcpy.JoinField_management(in_data=obsfinal, in_field='HYRIV_IDjoinedit',
                               join_table=atlassub1000, join_field='HYRIV_ID',
                               fields=['Shape_Length', 'dis_m3_pyr', 'UPLAND_SKM'])
    arcpy.AlterField_management(obsfinal, 'dis_m3_pyr', new_field_name='HYDROSHEDSdis', new_field_alias='HYDROSHEDSdis')
    arcpy.AlterField_management(obsfinal, 'UPLAND_SKM', new_field_name='HYDROSHEDSDA', new_field_alias='HYDROSHEDSDA')

# Extract upstream area and discharge for the point and ratio of point upstream area to line's pourpoint upstream area
if not arcpy.Exists(obsfinal_wgs):
    arcpy.Project_management(obsfinal, obsfinal_wgs, out_coor_system=arcpy.Describe(hydrolink).SpatialReference)
    arcpy.RepairGeometry_management(obsfinal_wgs)

extract_obsdisacc(in_obs=obsfinal_wgs,
                  in_net=atlassub,
                  obsin_field='HYRIV_IDjoinedit',
                  netjoin_field='HYRIV_ID',
                  in_diss=hydrodisras,
                  in_acc=hydroacc,
                  resdir=resdir,
                  delintermediate=False)

#Run a second time as sometimes a couple of obs don't compute
extract_obsdisacc(in_obs=obsfinal_wgs,
                  in_net=atlassub,
                  obsin_field='HYRIV_IDjoinedit',
                  netjoin_field='HYRIV_ID',
                  in_diss=hydrodisras,
                  in_acc=hydroacc,
                  resdir=resdir,
                  delintermediate=True)