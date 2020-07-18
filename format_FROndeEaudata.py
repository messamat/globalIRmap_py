from setup_localIRformatting import *

#Onde Eau dirs
datdir_onde = os.path.join(insitudatdir, 'OndeEau')
pathcheckcreate(datdir_onde)
resdir_onde = os.path.join(insituresdir, 'OndeEau.gdb')
pathcheckcreate(resdir_onde)


cartnet_raw =os.path.join(datdir_onde,'TronconHydrograElt_FXX.shp')
cartnet = os.path.join(resdir_onde,'carthage_network')
obsraw_onde= os.path.join(resdir_onde,'obsraw')
obsnodupli_onde = os.path.join(resdir_onde, 'obsnodupli')
ondeau_mergecsv = os.path.join(datdir_onde, 'onde_france_merge.csv')
obscartjoin_onde = os.path.join(resdir_onde, 'obs_cartnet_spatialjoin')
obscartjoinedit_onde = os.path.join(resdir_onde, 'obs_cartnet_spatialjoinedit')
obscartjoinclean_onde = os.path.join(resdir_onde, 'obs_cartnet_spatialjoinclean')
cartobsjoin_onde = os.path.join(resdir_onde, 'cartnet_obsclean_spatialjoin')
cartobsvert_onde = os.path.join(resdir_onde, 'cartnet_obsclean_spatialjoinvert')
cartobsjoin_ondeedit = os.path.join(resdir_onde, 'cartnet_obsclean_spatialjoinedit')
obsatlas_near = os.path.join(resdir_onde, 'obsclean_RiverAtlas_near')

atlassub_onde = os.path.join(resdir_onde, 'RiverAtlas_sub')
atlassubproj_onde = os.path.join(resdir_onde, 'RiverAtlas_subproj')
atlassub1000_onde = os.path.join(resdir_onde, 'RiverAtlascartnet_sub1000m')
atlasvert_onde = os.path.join(resdir_onde, 'RiverAtlascartnet_sub1000m_vertices')
atlascartobs_near = os.path.join(resdir_onde, 'RiverAtlASvertices_cartnet_obsclean_near')

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
    mergedel(dir=datdir_onde,
             repattern="onde_france_[0-9]{4}[.]csv",
             outfile=ondeau_mergecsv,
             returndf=False,
             delete=False,
             verbose=True)

# Convert to points (SpatialReference from ProjCoordSiteHydro referenced to http://mdm.sandre.eaufrance.fr/node/297134 - code 26 (247 observations also have code 2154 which must be a mistake and refers to the ESPG code)
if not arcpy.Exists(obsraw_onde):
    arcpy.MakeXYEventLayer_management(table=ondeau_mergecsv,
                                      in_x_field="<CoordXSiteHydro>",
                                      in_y_field="<CoordYSiteHydro>",
                                      out_layer='ondelyr',
                                      spatial_reference=arcpy.SpatialReference(
                                          2154))  # Lambert 93 for mainland France and Corsica
    arcpy.CopyFeatures_management('ondelyr', obsraw_onde)

#Delete duplicated locations
if not arcpy.Exists(obsnodupli_onde):
    group_duplishape(in_features=obsraw_onde, deletedupli=True, out_featuresnodupli=obsnodupli_onde)

#Join to Carthage network by Spatial Join
if not arcpy.Exists(obscartjoin_onde):
    arcpy.SpatialJoin_analysis(target_features=obsnodupli_onde,
                               join_features=cartnet,
                               out_feature_class=obscartjoin_onde,
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
arcpy.CopyFeatures_management(obscartjoin_onde, obscartjoinedit_onde)
arcpy.AddField_management(obscartjoinedit_onde, 'checkjoincart', 'SHORT')
arcpy.AddField_management(obscartjoinedit_onde, 'manualsnapcart', 'SHORT')

with arcpy.da.UpdateCursor(obscartjoinedit_onde,
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

[f.name for f in arcpy.ListFields(obscartjoin_onde)]

#Remove points to be deleted
arcpy.CopyFeatures_management(obscartjoinedit_onde, obscartjoinclean_onde)
with arcpy.da.UpdateCursor(obscartjoinclean_onde, ['manualsnapcart'])  as cursor:
    for row in cursor:
        if row[0] == 1:
            cursor.deleteRow()

#Snap points to Carthage network (max dist 10 m)
arcpy.Snap_edit(in_features=obscartjoinclean_onde, snap_environment=[[cartnet, 'EDGE', '10 meters']])

#Delete additional fields from points and spatially join Carthage network to points
for f in arcpy.ListFields(obscartjoinclean_onde):
    if f not in [f.name for f in arcpy.ListFields(obsraw_onde)]+['checkjoincart', 'manualsnapcart', 'distcart']:
        arcpy.DeleteField_management(obscartjoinclean_onde, f.name)

arcpy.SpatialJoin_analysis(target_features=cartnet,
                           join_features=obscartjoinclean_onde,
                           out_feature_class=cartobsjoin_onde,
                           join_operation="JOIN_ONE_TO_MANY",
                           join_type="KEEP_COMMON",
                           match_option="INTERSECT",
                           search_radius="0.1 meters")

#Convert Carthage segments to vertices (max 100 m apart)
arcpy.Densify_edit(in_features=cartobsjoin_onde, densification_method='DISTANCE',
                   distance='100 meters')
arcpy.FeatureVerticesToPoints_management(in_features=cartobsjoin_onde,
                                         out_feature_class=cartobsvert_onde,
                                         point_location='ALL')

########## JOIN CARTHAGE NETWORK WITH OBSERVATIONS TO RIVERATLAS ##############
#Select RiverATLAS segments within one kilometer from any Carthage reach of interest
#Re-project Carthage observations to WGS84
arcpy.CopyFeatures_management(arcpy.Describe(obscartjoinclean_onde).extent.polygon,
                              os.path.join(resdir_onde, 'obspoly'))
arcpy.Project_management(os.path.join(resdir_onde, 'obspoly'),
                         os.path.join(resdir_onde, 'obspolyproj'),
                         riveratlas)

#Select all RiverATLAS reaches within  50 kilometers from extent
arcpy.MakeFeatureLayer_management(riveratlas, 'atlaslyr')
arcpy.SelectLayerByLocation_management(in_layer='atlaslyr', #First select all RiverAtlas reaches intersecting the extent of the observations
                                       overlap_type="INTERSECT",
                                       select_features=os.path.join(resdir_onde, 'obspolyproj'),
                                       selection_type='NEW_SELECTION',
                                       search_distance='50 kilometers')
arcpy.CopyFeatures_management('atlaslyr', atlassub_onde)

#Project to Carthage projection
arcpy.Project_management(atlassub_onde, atlassubproj_onde, out_coor_system=cartnet)
arcpy.Delete_management('atlaslyr')

#select all reaches within 1000 meters from a Carthage-Onde Eau line observation
arcpy.MakeFeatureLayer_management(atlassubproj_onde, 'atlaslyr')
arcpy.SelectLayerByLocation_management(in_layer='atlaslyr',
                                       overlap_type="WITHIN_A_DISTANCE_GEODESIC",
                                       select_features=cartobsjoin_onde,
                                       search_distance='1000 meters',
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('atlaslyr', atlassub1000_onde)

#Convert HydroSHEDS segments to vertices
arcpy.Densify_edit(in_features=atlassub1000_onde, densification_method='DISTANCE',
                   distance='100 meters')
arcpy.FeatureVerticesToPoints_management(in_features=atlassub1000_onde,
                                         out_feature_class=atlasvert_onde,
                                         point_location='ALL')
vertdict = {row[0]: row[1] for row in arcpy.da.SearchCursor(atlasvert_onde, ['OID@', 'HYRIV_ID'])}
vertndict = defaultdict(int)
for k,v in vertdict.iteritems():
    vertndict[v] +=1

#Associate each Carthage-Onde observation line with all RiverATLAS lines within 1000 m

#Compute distance from Onde point observations to RiverATLAS lines
arcpy.GenerateNearTable_analysis(in_features=obscartjoinclean_onde,
                                 near_features=atlassub1000_onde,
                                 out_table= obsatlas_near,
                                 search_radius='1000 meters',
                                 location='location',
                                 closest='ALL',
                                 method='PLANAR')

#Compute distance from RiverATLAS vertices to nearest point on paired Cathage-Onde line
arcpy.GenerateNearTable_analysis(in_features=atlasvert_onde,
                                 near_features=cartobsjoin_onde,
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
cartoid = arcpy.Describe(cartobsjoin_onde).OIDFieldName
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

def get_arctab_df(in_tab, IDfield='OID@', valuefields='*', resetid = False):
    if IDfield == 'OID@':
        IDfield = arcpy.Describe(in_tab).OIDFieldName
    if valuefields == '*':
        valuefields = [f.name for f in arcpy.ListFields(in_tab)]
        valuefields.remove(IDfield)

    if isinstance(valuefields, list):
        fl = valuefields
        fl.insert(0, IDfield)
    else:
        fl = [IDfield, valuefields]

    ddict = {row[0]: row[1:] for row in arcpy.da.SearchCursor(in_tab, fl)}

    df = pd.DataFrame.from_dict(ddict, 'index')
    if resetid:
        df.reset_index(inplace=True)
    else:
        fl.remove(IDfield)

    df.columns = fl

    return(df)

obsiddf = get_arctab_df(in_tab=obscartjoinclean_onde, IDfield='OID@', valuefields='F_CdSiteHydro_', resetid=True)
atlasiddf = get_arctab_df(in_tab=atlassub1000_onde, IDfield='OID@', valuefields='HYRIV_ID', resetid=True)
obsatlasneardf = get_arctab_df(in_tab=obsatlas_near, IDfield='OID@', valuefields=['IN_FID', 'NEAR_FID', 'NEAR_DIST'])

nearptdf = obsatlasneardf.merge(obsiddf, left_on='IN_FID', right_on=arcpy.Describe(obscartjoinclean_onde).OIDFieldName).\
    merge(atlasiddf, left_on='NEAR_FID', right_on=arcpy.Describe(atlassub1000_onde).OIDFieldName)
nearptdf = nearptdf.drop([f for f in ['OBJECTID', 'OBJECTID_x', 'OBJECTID_y', 'IN_FID', 'NEAR_FID']
                          if f in nearptdf.columns], axis=1)
nearptdf.columns = ['OBStoATLAS_dist'] + list(nearptdf.columns[1:])

### Compute bearing
#Compute average bearing for each line
def projdistbearpt(x1, y1, x2, y2):
    # Distance between the two points
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    # Angle (counterclockwise)
    if y2-y1 > 0:
        angle = ((math.atan((x2 - x1) / (y2 - y1)) * 180 / math.pi) - 90)
    else:
        angle = -90

    # Revert direction (clockwise) to get bearing from north
    if angle > 0:
        bearing = 360 - angle
    else:
        bearing = -angle

    return ([dist, bearing])

#Append XY position of all vertices for each line
def linevert_wbearing(in_vert, in_IDfield):
    vertazdict = defaultdict(list)
    with arcpy.da.SearchCursor(in_vert, [in_IDfield, 'SHAPE@XY']) as cursor:
        for row in cursor:
            vertazdict[row[0]].append(row[1])

    #Replace vertazdict values with average bearing (out of 360)
    for k, v in vertazdict.iteritems():
        wanglesum = 0
        distsum = 0
        for v1, v2 in zip(v[:-1], v[1:]):
            distbear = projdistbearpt(v1[0], v1[1], v2[0], v2[1])
            wanglesum += distbear[0] * distbear[1]
            distsum += distbear[0]
        vertazdict[k] = wanglesum/float(distsum)

    return(vertazdict)

hydrobear = pd.DataFrame.from_dict(
    linevert_wbearing(in_vert=atlasvert_onde, in_IDfield='HYRIV_ID'),
    'index').reset_index()
hydrobear.columns = ['HYRIV_ID', 'hydrobearing']

cartbear = pd.DataFrame.from_dict(
    linevert_wbearing(in_vert = cartobsvert_onde, in_IDfield = 'IdTronconH'),
    'index').reset_index()
cartbear.columns = ['IdTronconH', 'cartbearing']

cart_objectidtroncon_df = get_arctab_df(in_tab=cartobsjoin_onde, IDfield='OID@',
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
nearbeardfmerge.sort_values(["IdTronconH", "matchindex"]).to_csv(os.path.join(insituresdir, 'netjoindistbear_onde.csv'))

#Join to Carthage subset
arcpy.AddField_management(cartobsjoin_onde, 'HYRIV_ID', 'LONG')
arcpy.AddField_management(cartobsjoin_onde, 'nvertu1k', 'LONG')
arcpy.AddField_management(cartobsjoin_onde, 'CARTtoATLAS_dist', 'FLOAT')
arcpy.AddField_management(cartobsjoin_onde, 'OBStoATLAS_dist', 'FLOAT')
arcpy.AddField_management(cartobsjoin_onde, 'bearingdiff', 'FLOAT')
arcpy.AddField_management(cartobsjoin_onde, 'matchindex', 'FLOAT')
nearbeardf_cartsub = nearbeardfmerge.sort_values(["IdTronconH", "matchindex"], ascending=False).\
    drop_duplicates(subset='IdTronconH')
with arcpy.da.UpdateCursor(cartobsjoin_onde, ['IdTronconH', 'HYRIV_ID',
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

#Inspect and edit (overlay cartobsjoin_ondeedit with atlassub1000_onde, and the original onde observations: see onde_checknetjoin.mxd
arcpy.CopyFeatures_management(cartobsjoin_onde, cartobsjoin_ondeedit)
arcpy.AddField_management(cartobsjoin_ondeedit, 'manualnetmatch', 'TEXT')

#Modification in the future:
#  - compute distance indices based on 1000 m (max dist) rather than based on min.
#  - check bearing function
# matchindex < 0.90
# OBStoATLASdist > 250 m
# CARTtoATLASdist > 350 m
# Total checked:

#Exclude all those that are too unclear, those >= 500 m from upstream limit of corresponding reach

#Copy to clean
#Remove -1, create new field to add edited HYRIV_ID, snap to specific feature (create fast function)
# compute how far down the line the site is as length and percentage of line's length
# Extract upstream area for the point and ratio of point upstream area to line's pourpoint upstream area


