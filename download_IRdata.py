"""
Download national river networks, model outputs, and on-the-ground observation of river intermittence
"""

from utility_functions import *


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

nhdpnw = os.path.join(resdir_pnw, 'NHDpnw')
nhdvaapnw = os.path.join(resdir_pnw, 'NHDvaapnw')

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

arcpy.DefineProjection_management(in_dataset=outnet_fr, coor_system=arcpy.SpatialReference(2192))

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
if not arcpy.Exists(obsnodupli_pnw):
    arcpy.CopyFeatures_management('pnwlyr', obsnodupli_pnw)

#Remove duplicate locations for spatial processing
duplidel = []
nunique = 0
with arcpy.da.UpdateCursor(obsnodupli_pnw, ['SHAPE@XY']) as cursor:
    for row in cursor:
        if row[0] not in duplidel:
            nunique += 1
            duplidel.append(row[0])
        else:
            print('Delete duplicate {}...'.format(row[0]))
            cursor.deleteRow()

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

#Delete extraneous records





http://www.freshwaterplatform.eu/
http://project.freshwaterbiodiversity.eu/
http://www.lifetrivers.eu/
http://irbas.inrae.fr/people/related-publications
http://www.ub.edu/fem/index.php/en/inici-riunet-en
https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e