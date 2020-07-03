"""
Download national river networks, model outputs, and on-the-ground observation of river intermittence
"""

from utility_functions import *
import math

rootdir = 'D://PhD//HydroATLAS'
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

arcpy.env.overwriteOutput = True
arcpy.env.qualifiedFieldNames = 'False'
arcpy.CheckOutExtension('Spatial')

compdatdir = os.path.join(datdir, 'Comparison_databases')
compresdir = os.path.join(resdir, 'Comparison_databases')
pathcheckcreate(compdatdir)
pathcheckcreate(compresdir)

insitudatdir = os.path.join(datdir, 'Insitu_databases')
insituresdir = os.path.join(resdir, 'Insitu_databases')
pathcheckcreate(compdatdir)
pathcheckcreate(compresdir)

hydrobasin12 = os.path.join(datdir, 'Bernhard\\HydroATLAS\\HydroATLAS_v10_final_data\\'
                                    'BasinATLAS_v10_shp\\BasinATLAS_v10_lev12.shp')
riveratlas = os.path.join(datdir, 'Bernhard\\HydroATLAS\\HydroATLAS_v10_final_data\\',
                          'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
hydroacc = os.path.join(datdir, 'Bernhard\HydroATLAS\HydroATLAS_Geometry\Accu_area_grids\upstream_area_skm_15s.gdb',
                        'up_area_skm_15s')

#US dirs
datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

inbas_us = os.path.join(datdir_us, 'WBD_National_GDB.gdb', 'WBD', 'WBDHU8')
outnet_us = os.path.join(resdir_us, 'network')
bassel_us = os.path.join(resdir_us, 'hydrobasins12')

#France dirs
datdir_fr = os.path.join(compdatdir, 'France')
pathcheckcreate(datdir_fr)
resdir_fr = os.path.join(compresdir, 'france.gdb')
pathcheckcreate(resdir_fr)
outnet_fr = os.path.join(resdir_fr, 'network')
bassel_fr = os.path.join(resdir_fr, 'hydrobasins12')

#Onde Eau dirs
datdir_onde = os.path.join(insitudatdir, 'OndeEau')
pathcheckcreate(datdir_onde)
resdir_onde = os.path.join(insituresdir, 'OndeEau.gdb')
pathcheckcreate(resdir_onde)

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

#OndeEau data processing
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

atlassub_onde = os.path.join(resdir_onde, 'RiverAtlas_sub')
atlassubproj_onde = os.path.join(resdir_onde, 'RiverAtlas_subproj')
atlassub1000_onde = os.path.join(resdir_onde, 'RiverAtlascartnet_sub1000m')
atlasvert_onde = os.path.join(resdir_onde, 'RiverAtlascartnet_sub1000m_vertices')
atlascartobs_near = os.path.join(resdir_onde, 'RiverAtlASvertices_cartnet_obsclean_near')


# -------------------- FRANCE - ONDE EAU -------------------------------------------------------------------------------
for yr in range(2012, 2021):
    in_url = "https://onde.eaufrance.fr/sites/default/files/fichiers-telechargeables/onde_france_{}.zip".format(yr)
    dlfile(in_url,
           outpath=datdir_onde,
           outfile=os.path.split(in_url)[1],
           ignore_downloadable=True)

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

# Download Carthage Hydrographic network
dlfile(
    "http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2017/France_metropole_entiere/SHP/TronconHydrograElt_FXX-shp.zip",
    outpath=datdir_onde,
    outfile="TronconHydrograElt_FXX-shp.zip",
    ignore_downloadable=True)

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

#select all reaches within 1000 meters from a Carthage-Onde Eau observation
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

# #Compute one-to-many join table by computing average distance of vertices to line (assuming 5000 m dist for those > 5000 m)
# for hyriv_id, hyriv_dict in neardict.iteritems():
#     print(hyriv_id)
#     for cart_id, cart_dict in hyriv_dict.iteritems():
#         nvert_u1k = neardict[hyriv_id][cart_id][1]
#         neardict[hyriv_id][cart_id] = ((neardict[hyriv_id][cart_id][0] + 5000 * (vertndict[hyriv_id] - nvert_u5k)) /
#                                        vertndict[hyriv_id])

#Convert to dataframe and write out if needs manual inspection later
cartoid = arcpy.Describe(cartobsjoin_onde).OIDFieldName
neardf = pd.concat({
    k: pd.DataFrame.from_dict(v, 'index') for k, v in neardict.items()},
    axis=0).reset_index()
neardf.columns = ['HYRIV_ID', cartoid, 'distsum', 'nvertu1k']
neardf['dist'] = neardf['distsum']/neardf['nvertu1k'].astype(float)

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

cart_objectidtroncon_df = pd.DataFrame.from_dict(
    {row[0]: int(row[1]) for row in  arcpy.da.SearchCursor(cartobsjoin_onde, ['OID@', 'IdTronconH'])},
    'index').\
    reset_index()
cart_objectidtroncon_df.columns = [cartoid, 'IdTronconH']

#Join distance df with bearing dfs
nearbeardf = neardf.set_index(cartoid). \
    join(cart_objectidtroncon_df.set_index(cartoid)). \
    set_index('HYRIV_ID'). \
    join(hydrobear.set_index('HYRIV_ID')). \
    reset_index(). \
    set_index('IdTronconH'). \
    join(cartbear.set_index('IdTronconH'))


mindistdf = pd.DataFrame(nearbeardf.groupby("IdTronconH")['dist'].min())
nearbeardf = nearbeardf.join(mindistdf, rsuffix='min')

maxvertdf = pd.DataFrame(nearbeardf.groupby("IdTronconH")['nvertu1k'].max())
nearbeardf = nearbeardf.join(maxvertdf, rsuffix='max').reset_index()

#Compute absolute difference in bearing between HydroSHEDS and Carthage lines
nearbeardf['bearingdiff'] = abs(nearbeardf['hydrobearing'] - nearbeardf['cartbearing'])
nearbeardf['distindex'] = nearbeardf['distmin']/nearbeardf['dist']
nearbeardf['nvertindex'] = nearbeardf['nvertu1k']/nearbeardf['nvertu1kmax']
nearbeardf['matchindex'] = (3*nearbeardf['distindex'] +
                            1-(nearbeardf['bearingdiff']/180.0) +
                            2*nearbeardf['nvertindex'])/6

#Export to CSV
nearbeardf.sort_values(["IdTronconH", "matchindex"]).to_csv(os.path.join(insituresdir, 'netjoindistbear_onde.csv'))

#Join to Carthage subset
arcpy.AddField_management(cartobsjoin_onde, 'HYRIV_ID', 'LONG')
arcpy.AddField_management(cartobsjoin_onde, 'nvertu1k', 'LONG')
arcpy.AddField_management(cartobsjoin_onde, 'dist', 'FLOAT')
arcpy.AddField_management(cartobsjoin_onde, 'bearingdiff', 'FLOAT')
arcpy.AddField_management(cartobsjoin_onde, 'matchindex', 'FLOAT')
nearbeardf_cartsub = nearbeardf.sort_values(["IdTronconH", "matchindex"], ascending=False).\
    drop_duplicates(subset='IdTronconH')
with arcpy.da.UpdateCursor(cartobsjoin_onde, ['IdTronconH', 'HYRIV_ID',
                                              'nvertu1k', 'dist', 'bearingdiff', 'matchindex']) as cursor:
    for row in cursor:
        dfrow = nearbeardf_cartsub.loc[nearbeardf_cartsub['IdTronconH'] == row[0]]
        if len(dfrow['HYRIV_ID']) != 0:
            row[1] = dfrow['HYRIV_ID'].values[0]
            row[2] = dfrow['nvertu1k'].values[0]
            row[3] = dfrow['dist'].values[0]
            row[4] = dfrow['bearingdiff'].values[0]
            row[5] = dfrow['matchindex'].values[0]
            cursor.updateRow(row)

#Inspect and edit
arcpy.CopyFeatures(cartobsjoin_onde, )








#Join to HydroSHEDS subset
arcpy.AddField_management(atlassub1000_onde, 'IdTronconH', 'LONG')
arcpy.AddField_management(atlassub1000_onde, 'dist', 'FLOAT')
arcpy.AddField_management(atlassub1000_onde, 'bearingdiff', 'FLOAT')
arcpy.AddField_management(atlassub1000_onde, 'matchindex', 'SHORT')

#Thoughts for river network matching:
#
# 	* Digitize to points: compute near distance to all lines within a given distance. Average distance to each nearby line.
# 	*
# For every subsegment, compute azimuth, average azimuth as weighted average by length
# 	*
# If can, get river order for each network
# 	*
# Slope from DEM in native resolution
#
#
#
# Maybe have a set of training and test sets of matched segments to then develop a weighting scheme through optimization

####################### REGIONAL RIVER NETWORK DATASETS WITH INTERMITTENCE DATA ########################################
#-------------------- USA ----------------------------------------------------------------------------------------------
#Extra links
#"https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/National/HighResolution/GDB/NHD_H_National_GDB.zip"
#"https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHDPlusHR/Beta/GDB/"

#Download NHDplus data
# nhdplus_linklist = ["https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlusHR/Beta/" \
#                     "GDB/NHDPLUS_H_{}_HU4_GDB.zip".format(str(i).zfill(4)) for i in xrange(1807, 2103)]
# for url in nhdplus_linklist:
#     outgdb = os.path.join(datdir_us, "{}.gdb".format(os.path.splitext(os.path.split(url)[1])[0]))
#     if not arcpy.Exists(outgdb):
#         print(url)
#         dlfile(url=url, outpath=datdir_us, ignore_downloadable=True)

#Download HUC4 polygons
# if not os.path.exists('WBD_National_GDB.zip'):
#     dlfile(url='https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip',
#            outpath=datdir_us,
#            outfile='WBD_National_GDB.zip',
#            ignore_downloadable=True)

#Merge all flow lines
FLlist = [os.path.join(i, 'NHDFlowLine') for i in 
          getfilelist(dir=datdir_us, repattern='Hydrography', gdbf=True, nongdbf=False)]

#Add attributes
#Drainage area
fldl_us = defaultdict(list)
for tab in FLlist:
    print(tab)
    with arcpy.da.SearchCursor(tab, ['NHDPlusID', 'FCode', 'LengthKM']) as cursor:
        for row in cursor:
            fldl_us[row[0]] = [row[1], row[2]]

VAlist = getfilelist(dir=datdir_us, repattern='NHDPlusFlowlineVAA', gdbf=True, nongdbf=False)
for tab in VAlist:
    print(tab)
    with arcpy.da.SearchCursor(tab, ['NHDPlusID', 'ReachCode', 'StreamOrde', 'TotDASqKm']) as cursor:
        for row in cursor:
            fldl_us[row[0]].extend([row[1], row[2], row[3]])

#Discharge
EROMQAlist = getfilelist(dir=datdir_us, repattern='NHDPlusEROMMA', gdbf=True, nongdbf=False)
for tab in EROMQAlist:
    print(tab)
    with arcpy.da.SearchCursor(tab, ['NHDPlusID', 'QAMA']) as cursor:
        for row in cursor:
            fldl_us[row[0]].append(row[1])

usnhd_df = pd.DataFrame.from_dict(fldl_us, orient='index')
usnhd_df.columns = ['FCode', 'LengthKM', 'ReachCode', 'StreamOrde', 'TotDASqKm', 'QAMA']
usnhd_df.to_csv(os.path.join(datdir_us, 'NHD_attris.csv'))

#Create a subselection of HydroSHEDS river basins that overlap with the NHD dataset
arcpy.Project_management(inbas_us, out_dataset=os.path.join(resdir_us, 'WBDHU4_wgs84'),
                         out_coor_system=hydrobasin12)

arcpy.MakeFeatureLayer_management(hydrobasin12, 'hydrous')
arcpy.SelectLayerByLocation_management('hydrous', overlap_type='INTERSECT',
                                       select_features=os.path.join(resdir_us, 'WBDHU4_wgs84'),
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('hydrous', bassel_us)

#Then create spatial join whereby each HydroBASINS level 12 is associated with the HUC8 that it overlaps with the most
basinters_us = os.path.join(resdir_us, 'bassinters')
arcpy.Intersect_analysis(in_features=[bassel_us, inbas_us], out_feature_class=basinters_us, join_attributes='ALL')
arcpy.AddGeometryAttributes_management(basinters_us, Geometry_Properties='AREA_GEODESIC', Area_Unit='SQUARE_KILOMETERS')
intersdict_us = {}
with arcpy.da.SearchCursor(basinters_us, ['HYBAS_ID','HUC8','AREA_GEO']) as cursor:
    for row in cursor:
        if row[0] in intersdict_us:
            if row[2] > intersdict_us[row[0]][1]:
                intersdict_us[row[0]] = [row[1], row[2]]
        else:
            intersdict_us[row[0]] =[row[1], row[2]]

if "HUC8" not in [f.name for f in arcpy.ListFields(bassel_us)]:
    arcpy.AddField_management(bassel_us, 'HUC8', 'TEXT')

with arcpy.da.UpdateCursor(bassel_us, ['HUC8', 'HYBAS_ID']) as cursor:
    for row in cursor:
        row[0] = intersdict_us[row[1]][0]
        cursor.updateRow(row)

#-------------------- France ----------------------------------------------------------------------------------------------
#Data are personal communications from Snelder et al. 2013

#Join original French river network data with Snelder et al.' predictions
arcpy.MakeFeatureLayer_management(os.path.join(datdir_fr, 'rhtvs2_all_phi_qclass.shp'),
                                  out_layer='frlyr')
arcpy.AddJoin_management(in_layer_or_view='frlyr', in_field='ID_DRAIN',
                         join_table=os.path.join(datdir_fr, 'INT_RF.txt'), join_field='AllPred$ID_DRAIN')
arcpy.CopyFeatures_management(in_features='frlyr', out_feature_class=outnet_fr)

arcpy.DefineProjection_management(in_dataset=outnet_fr, coor_system=arcpy.SpatialReference(2192)) #ED_1950_France_EuroLambert

#Create a subselection of HydroSHEDS river sections that overlap with the French dataset
arcpy.MakeFeatureLayer_management(hydrobasin12, 'hydrofr')
arcpy.SelectLayerByLocation_management('hydrofr', overlap_type='CONTAINS', select_features=outnet_fr,
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('hydrofr', bassel_fr)

#-------------------- Brazil ----------------------------------------------------------------------------------------------
#ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/shapefile/
#ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/informacoes_tecnicas/
"ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/shapefile/bcim_2016_shapefiles_21-11-2018.zip"
"ftp://geoftp.ibge.gov.br/cartas_e_mapas/mapa_indice_digital_4ed/produto_mapa_indice_digital/" #For topographic map footprints
#Select streams based on IMPRESSA_250000 in Shap

#Argentina
"https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG"

#Espana
#https://centrodedescargas.cnig.es/CentroDescargas/buscadorCatalogo.do?codFamilia=HIDRO#

#Canada
"http://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/fgdb/Hydro/"
"http://ftp.maps.canada.ca/pub/nrcan_rncan/vector/canvec/fgdb/Hydro/canvec_50K_AB_Hydro_fgdb.zip"
"http://ftp.maps.canada.ca/pub/nrcan_rncan/vector/index/nts_snrc.zip" #Tile system

#Mexico
#https://www.inegi.org.mx/temas/hidrografia/default.html#Descargas - REFERENCE
mexico_linklist = ["http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/" \
                   "bvinegi/productos/geografia/hidrogeolo/region_hidrografica/70282500{}_s.zip".format(i) for i in xrange(6976, 7013)]
for url in mexico_linklist:
    print(url)
    dlfile(url=url, outpath=os.path.join(compdatdir, 'Mexico'), ignore_downloadable=True)

"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007012_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007011_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007010_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007009_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007008_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007007_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007006_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007005_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825007004_s.zip"

"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825006977_s.zip"
"http://internet.contenidos.inegi.org.mx/contenidos/Productos/prod_serv/contenidos/espanol/bvinegi/productos/geografia/hidrogeolo/region_hidrografica/702825006976_s.zip"

#France
http://onde.eaufrance.fr/
https://www.sciencebase.gov/catalog/item/5a0f338de4b09af898d099b9

#Digital Chart of the World
#https://guides.libraries.psu.edu/c.php?g=376207&p=5296088
#https://psu.app.box.com/v/dcw

#South Africa

####################### HIGH QUALITY IN SITU DATA ######################################################################
#-------------------- PNW-PROSPER ----------------------------------------------------------------------------------------------
#Non-perennial includes classifications such as dry channel, no flow, intermittent, and ephemeral

#https://www.sciencebase.gov/catalog/item/5a0f338de4b09af898d099b9 - REFERENCE PAGE
if not os.path.exists(os.path.join(datdir_pnw, 'StreamflowPermObs.zip')):
    pnw_r = urllib2.urlopen("https://www.sciencebase.gov/catalog/item/5a0f338de4b09af898d099b9")
    pnw_soup = BeautifulSoup(pnw_r, features="html.parser")

    for link in pnw_soup.findAll("span", {"class": "sb-file-get sb-download-link"}):
        ltext = link.text
        if ltext == 'StreamflowPermObs.zip':
            dlfile(url=urlparse.urljoin(pnw_r.url, link.get('data-url')),
                   outpath = datdir_pnw,
                   outfile = ltext,
                   ignore_downloadable= True)

if not os.path.exists(os.path.join(datdir_pnw, 'fac_taudem_17all_int.tif')):
    pnw_r = urllib2.urlopen("https://www.sciencebase.gov/catalog/item/5a789cc3e4b00f54eb1e8390")
    pnw_soup = BeautifulSoup(pnw_r, features="html.parser")
    for link in pnw_soup.findAll("span", {"class": "sb-file-get sb-download-link"}):
        ltext = link.text
        if ltext in ['ContributingArea_CPG_17all.xml', 'fac_taudem_17all_int.tif']:
            dlfile(url=urlparse.urljoin(pnw_r.url, link.get('data-url')),
                   outpath = datdir_pnw,
                   outfile = ltext,
                   ignore_downloadable= True)


#Only keep points for which use is:
arcpy.MakeFeatureLayer_management(obsraw_pnw, out_layer='pnwlyr',
                                  where_clause='(Use IN {}) AND (Month > 6)'.format(str(("No; Outside 2004-2016", "Yes"))))
if not arcpy.Exists(obsub_pnw):
    arcpy.CopyFeatures_management('pnwlyr', obsub_pnw)

#Flag groups of duplicates in obsub_pnw (utility_functions)
group_duplishape(in_features=obsub_pnw, deletedupli=True, out_featuresnodupli=obsnodupli_pnw)

#Remove all points within 1 km of each other
if not arcpy.Exists(sparseattri_pnw):
    arcpy.CreateRandomPoints_management(out_path=os.path.split(sparse_pnw)[0],
                                        out_name=os.path.split(sparse_pnw)[1],
                                        constraining_feature_class=obsnodupli_pnw,
                                        number_of_points_or_field=nunique,
                                        minimum_allowed_distance='1000 meters')
    arcpy.MakeFeatureLayer_management(sparse_pnw, 'sparselyr')
    arcpy.AddJoin_management('sparselyr', in_field='CID',
                             join_table=obsnodupli_pnw, join_field=arcpy.Describe(obsnodupli_pnw).OIDFieldName)
    arcpy.CopyFeatures_management('sparselyr', sparseattri_pnw)

#Get NHDplus for all HUC4 basins extracted from HUC_8 and merge them
huc4list = {row[0][0:4] for row in arcpy.da.SearchCursor(sparseattri_pnw, 'HUC_8')}
if not arcpy.Exists(nhdpnw):
    nhdlist = [os.path.join(datdir_us, "NHDPLUS_H_{0}_HU4_GDB.gdb".format(huc4), "Hydrography", "NHDFlowLine") for huc4 in huc4list]
    arcpy.Merge_management(nhdlist, nhdpnw)

#Spatial join points to NHDplus high res
if not arcpy.Exists(obsjoin_pnw):
    arcpy.SpatialJoin_analysis(sparseattri_pnw, nhdpnw, out_feature_class=obsjoin_pnw, join_operation="JOIN_ONE_TO_ONE",
                               join_type='KEEP_ALL', match_option='CLOSEST_GEODESIC', distance_field_name='distnhd')

#Get flow accumulation and convert it to km2 (rather than pixels)
ExtractMultiValuesToPoints(obsjoin_pnw, in_rasters=[[facraw_pnw, 'fac']])
arcpy.AddField_management(obsjoin_pnw, 'fac_km2', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_pnw, ['fac', 'fac_km2']) as cursor:
    for row in cursor:
        if row[0] is not None:
            row[1] = row[0]*900/float(10**6)
            cursor.updateRow(row)

#Spatial join to RiverATLAS and compute ratio of drainage areas
obsjoin_pnwproj = "{}_proj".format(obsjoin_pnw)
arcpy.Project_management(obsjoin_pnw, obsjoin_pnwproj, out_coor_system=riveratlas)
if not arcpy.Exists(obsjoin_pnw):
    arcpy.SpatialJoin_analysis(obsjoin_pnwproj, riveratlas,
                               out_feature_class=obsjoin_atlas, join_operation="JOIN_ONE_TO_ONE",
                               join_type='KEEP_ALL', match_option='CLOSEST_GEODESIC', distance_field_name='distatlas')

arcpy.AddField_management(obsjoin_atlas, 'facratio', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_atlas, ['fac_km2', 'UPLAND_SKM', 'facratio']) as cursor:
    for row in cursor:
        if row[0] is not None:
            row[2] = row[1]/float(row[0])
            cursor.updateRow(row)

#Remove those that are over 250 m away and NHDplus upstream area < 10 km2
arcpy.MakeFeatureLayer_management(obsjoin_atlas, 'atlasjoinlyr',
                                  where_clause='NOT ((((distatlas > 500) OR (facratio > 3)) AND (fac_km2 < 10)) OR (fac_km2 Is Null))')
arcpy.CopyFeatures_management('atlasjoinlyr', obsjoin_atlassub)

#When multiple points are associated with a single HydroATLAS reach, keep the one with the lowest facratio
duplidict = defaultdict(float)
with arcpy.da.SearchCursor(obsjoin_atlassub, ['HYRIV_ID', 'facratio']) as cursor:
    for row in cursor:
        if (row[0] not in duplidict) or (abs(row[1] - 1) < duplidict[row[0]]):
            duplidict[row[0]] = (abs(row[1] - 1))
with arcpy.da.UpdateCursor(obsjoin_atlassub, ['HYRIV_ID', 'facratio']) as cursor:
    for row in cursor:
        if (abs(row[1] - 1)) > duplidict[row[0]]:
            print('Delete {}'.format(row[0]))
            cursor.deleteRow()

#Snap all those that are within 200 m from a RiverATLAS segment
arcpy.CopyFeatures_management(obsjoin_atlassub, obsjoin_atlasnap)
snapenv = [[riveratlas, 'EDGE', '200 meters']]
arcpy.Snap_edit(obsjoin_atlasnap, snapenv)

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

#Extract flow accumulation directly from HydroSHEDS flow accumulation grid
ExtractMultiValuesToPoints(obsjoin_atlasnapclean, in_rasters=hydroacc, bilinear_interpolate_values='NONE')

#Compute a second flow accumulation ratio for checking
arcpy.AddField_management(obsjoin_atlasnapclean, 'facratio_2', 'FLOAT')
with arcpy.da.UpdateCursor(obsjoin_atlasnapclean, ['facratio_2', 'fac_km2', 'up_area_skm_15s']) as cursor:
    for row in cursor:
        if row[2] > 0:
            row[0] = row[1]/float(row[2])
        else:
            row[0] = -9999.0
        cursor.updateRow(row)

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

#######################################################################################################################
# http://www.freshwaterplatform.eu/
# http://project.freshwaterbiodiversity.eu/
# http://www.lifetrivers.eu/
# http://irbas.inrae.fr/people/related-publications
# http://www.ub.edu/fem/index.php/en/inici-riunet-en
# https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e