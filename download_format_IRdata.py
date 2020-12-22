"""
Download national river networks, model outputs, and on-the-ground observation of river intermittence
"""

from setup_localIRformatting import *
import math
arcpy.CheckOutExtension('Network')

#US dirs
datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

inbas_us = os.path.join(datdir_us, 'WBD_National_GDB.gdb', 'WBD', 'WBDHU8')
outnet_us = os.path.join(resdir_us, 'network')
bassel_us = os.path.join(resdir_us, 'hydrobasins12')
netrefsub_us = os.path.join(resdir_us, 'netref_o10')

#France dirs
datdir_fr = os.path.join(compdatdir, 'France')
datdir_snelder = os.path.join(datdir_fr, 'Snelder')
pathcheckcreate(datdir_snelder)
resdir_snelder = os.path.join(compresdir, 'france_snelder.gdb')
pathcheckcreate(resdir_snelder)

datdir_snelder = os.path.join(datdir_fr, 'Snelder')
pathcheckcreate(datdir_snelder)
resdir_snelder = os.path.join(compresdir, 'france_snelder.gdb')
pathcheckcreate(resdir_snelder)
outnet_snelder = os.path.join(resdir_snelder, 'network')
netrefsub_snelder = os.path.join(resdir_snelder, 'netref_o10')
bassel_fr = os.path.join(resdir_snelder, 'hydrobasins12')

#Brazil dirs
datdir_bz = os.path.join(compdatdir, 'Brazil')
resdir_bz= os.path.join(compresdir, 'brazil.gdb')
pathcheckcreate(datdir_bz)
pathcheckcreate(resdir_bz)

#Australia
datdir_au = os.path.join(compdatdir, 'Australia')
resdir_au= os.path.join(compresdir, 'australia.gdb')
pathcheckcreate(datdir_au)
pathcheckcreate(resdir_au)
netref_au = os.path.join(datdir_au, 'SH_Network_GDB' ,'SH_Network.gdb', 'AHGFNetworkStream')
netproj_au = os.path.join(resdir_au, 'AHGFNetworkStream_WGS84')

#OndeEau dir
datdir_onde = os.path.join(insitudatdir, 'OndeEau')
pathcheckcreate(datdir_onde)

#PNW
datdir_pnw = os.path.join(insitudatdir, 'PNW')
pathcheckcreate(datdir_pnw)

####################### REGIONAL RIVER NETWORK DATASETS WITH INTERMITTENCE DATA ########################################
#-------------------- USA ----------------------------------------------------------------------------------------------
#Make sure to only keep data for the conterminous US
#Subselect HydroSHEDS river basins that overlap with the conterminous US
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

#------NHDPlus High Resolution ----------------------
#Extra links
#"https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/National/HighResolution/GDB/NHD_H_National_GDB.zip"
#"https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHDPlusHR/Beta/GDB/"

#Download NHDplus data
# nhdplus_linklist = ["https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlusHR/Beta/" \
#                     "GDB/NHDPLUS_H_{}_HU4_GDB.zip".format(str(i).zfill(4)) for i in xrange(1807, 2103)]
# for url in nhdplus_linklist:
#     outgdb = os.path.join(datdir_us, "NHDPlus_hr", "{}.gdb".format(os.path.splitext(os.path.split(url)[1])[0]))
#     if not arcpy.Exists(outgdb):
#         print(url)
#         dlfile(url=url, outpath=os.path.join(datdir_us, "NHDPlus_hr"), ignore_downloadable=True)

#Download HUC4 polygons
# if not os.path.exists('WBD_National_GDB.zip'):
#     dlfile(url='https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip',
#            outpath=datdir_us,
#            outfile='WBD_National_GDB.zip',
#            ignore_downloadable=True)

#Merge all flow lines
FLlist = [os.path.join(i, 'NHDFlowLine') for i in 
          getfilelist(dir=os.path.join(datdir_us, "NHDPlus_hr"), repattern='Hydrography', gdbf=True, nongdbf=False)]

#Add attributes
#Drainage area
fldl_us = defaultdict(list)
for fltab in FLlist:
    print(fltab)
    with arcpy.da.SearchCursor(fltab, ['NHDPlusID', 'FCode', 'LengthKM']) as cursor:
        for row in cursor:
            fldl_us[row[0]] = [row[1], row[2]]

VAlist = getfilelist(dir=os.path.join(datdir_us, "NHDPlus_hr"), repattern='NHDPlusFlowlineVAA', gdbf=True, nongdbf=False)
flsubl = []
for tab in VAlist:
    print(tab)
    with arcpy.da.SearchCursor(tab, ['NHDPlusID', 'ReachCode', 'StreamOrde', 'TotDASqKm']) as cursor:
        for row in cursor:
            fldl_us[row[0]].extend([row[1], row[2], row[3]])

    #Join data to line dataset
    fltab = os.path.join(os.path.split(tab)[0], "Hydrography", "NHDFlowLine")

    HUCID = re.search("(?<={}[\\\])NHDPLUS_H_[0-9]*(?=_HU4_GDB[.]gdb[\\\]Hydrography[\\\]NHDFlowLine)".format(
        os.path.split(os.path.join(datdir_us, "NHDPlus_hr"))[1]), fltab).group()
    outflsub = os.path.join(resdir_us, "{}_FlowLinesub".format(HUCID))
    flsubl.append(outflsub)

    if not arcpy.Exists(outflsub):
        arcpy.MakeFeatureLayer_management(fltab, 'fllyr',
                                          where_clause= "(FType IN (460, 558)) AND (InNetwork = 1)")
        arcpy.CopyFeatures_management('fllyr', outflsub)

    fltab_fnames = [f.name for f in arcpy.ListFields(outflsub)]
    if not 'StreamOrde' in fltab_fnames:
        arcpy.AddField_management(outflsub,'StreamOrde', 'SHORT')
    if not 'TotDASqKm' in fltab_fnames:
        arcpy.AddField_management(outflsub, 'TotDASqKm', 'DOUBLE')

    with arcpy.da.UpdateCursor(outflsub, ['NHDPlusID', 'StreamOrde', 'TotDASqKm']) as cursor:
        for row in cursor:
            if len(fldl_us[row[0]]) > 4:
                if fldl_us[row[0]][4] < 10:
                    cursor.deleteRow()
                else:
                    row[1] = fldl_us[row[0]][3]
                    row[2] = fldl_us[row[0]][4]
                    cursor.updateRow(row)

#Discharge
EROMQAlist = getfilelist(dir=datdir_us, repattern='NHDPlusEROMMA', gdbf=True, nongdbf=False)
for tab in EROMQAlist:
    print(tab)
    with arcpy.da.SearchCursor(tab, ['NHDPlusID', 'QAMA']) as cursor:
        for row in cursor:
            fldl_us[row[0]].append(row[1])

usnhd_df = pd.DataFrame.from_dict(fldl_us, orient='index')
usnhd_df.columns = ['FCode', 'LengthKM', 'ReachCode', 'StreamOrde', 'TotDASqKm', 'QAMA']
usnhd_df.to_csv(os.path.join(datdir_us, 'NHDhr_attris.csv'))

#Subsetlect all NHD segments with TotDASqKm > 10km2 and merge them to visually assess predictions
arcpy.Merge_management(inputs=flsubl, output=netrefsub_us)

#Divide by river order
for ord in [[1,3], [4,5], [6,7], [8,9], [10, 11]]:
    subriver_out = os.path.join(resdir_us,
                                '{0}_ORD{1}_{2}'.format(os.path.split(netrefsub_us)[1], ord[0], ord[1]))
    if not arcpy.Exists(subriver_out):
        sqlexp = 'StreamOrde >= {0} AND StreamOrde <= {1}'.format(ord[0], ord[1])
        print(sqlexp)
        arcpy.MakeFeatureLayer_management(netrefsub_us, out_layer='subriver', where_clause=sqlexp)
        arcpy.CopyFeatures_management('subriver', subriver_out)
        arcpy.Delete_management('subriver')


#------ NHDPlus Medium Resolution ------------------
#https://www.epa.gov/waterdata/nhdplus-national-data
#https://s3.amazonaws.com/edap-nhdplus/NHDPlusV21/Data/NationalData/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z
datdir_NHDmr = os.path.join(datdir_us, "NHDPlus_mr")
NHDmr = os.path.join(datdir_NHDmr,
                     "NHDPlusNationalData\NHDPlusV21_National_Seamless_Flattened_Lower48.gdb\NHDSnapshot",
                     "NHDFlowline_Network")

#Export attributes for R analysis
arcpy.MakeFeatureLayer_management(NHDmr, 'NHDmrlyr',
                                  where_clause="FType IN ('ArtificialPath', 'StreamRiver')")
arcpy.CopyRows_management('NHDmrlyr', os.path.join(datdir_us, 'NHDmr_attris.csv'))

#Export subset > 10 km2 natural rivers for mapping
arcpy.SelectLayerByAttribute_management('NHDmrlyr', selection_type='SUBSET_SELECTION',
                                  where_clause="TotDASqKm >= 10")
outflsubmr = os.path.join(resdir_us, "NHDPlusmr_FlowLinesub")
if not arcpy.Exists(outflsubmr):
    arcpy.CopyFeatures_management('NHDmrlyr', outflsubmr)

#Add discharge field in m3/s
arcpy.AddField_management(outflsubmr, 'QE_MAm3', field_type='DOUBLE')
with arcpy.da.UpdateCursor(outflsubmr, ['QE_MAm3', 'QE_MA']) as cursor:
    for row in cursor:
        row[0] = float(row[1]) * 0.028316847
        cursor.updateRow(row)


for dis in [[0,0.1], [0.1,1], [1,10], [10,100], [100,1000], [1000, 10000], [10000, 1000000]]:
    subriver_out = os.path.join(resdir_us,
                                '{0}_DIS{1}_{2}'.format(os.path.split(outflsubmr)[1],
                                                        re.sub('[.]', '', str(dis[0])),
                                                        re.sub('[.]', '', str(dis[1]))))
    sqlexp = 'QE_MAm3 >= {0} AND QE_MAm3 <= {1}'.format(dis[0], dis[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(outflsubmr, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

#-------------------- France - SNELDER ----------------------------------------------------------------------------------------------
#Data are personal communications from Snelder et al. 2013

#Join original French river network data with Snelder et al.' predictions
arcpy.MakeFeatureLayer_management(os.path.join(datdir_snelder, 'rhtvs2_all_phi_qclass.shp'),
                                  out_layer='frlyr')
arcpy.AddJoin_management(in_layer_or_view='frlyr', in_field='ID_DRAIN',
                         join_table=os.path.join(datdir_snelder, 'INT_RF.txt'), join_field='AllPred$ID_DRAIN')
arcpy.CopyFeatures_management(in_features='frlyr', out_feature_class=outnet_snelder)

arcpy.DefineProjection_management(in_dataset=outnet_snelder, coor_system=arcpy.SpatialReference(2192)) #ED_1950_France_EuroLambert

#Subsetlect all segments with drainage area > 10 km2
arcpy.MakeFeatureLayer_management(outnet_snelder,
                                  out_layer='frlyr', where_clause='rhtvs2_all_phi_qclass_SURF_BV >= 10')
arcpy.CopyFeatures_management(in_features='frlyr', out_feature_class=netrefsub_snelder)

#Divide french network by river order
for ord in [[1,3], [4,5], [6,7], [8,9]]:
    subriver_out = os.path.join(resdir_snelder,
                                '{0}_ORD{1}_{2}'.format(os.path.split(netrefsub_snelder)[1], ord[0], ord[1]))
    sqlexp = 'rhtvs2_all_phi_qclass_STRAHLER >= {0} AND rhtvs2_all_phi_qclass_STRAHLER <= {1}'.format(ord[0], ord[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(netrefsub_snelder, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

for dis in [[0,0.1], [0.1,1], [1,10], [10,100], [100,1000], [1000, 10000]]:
    subriver_out = os.path.join(resdir_snelder,
                                '{0}_DIS{1}_{2}'.format(os.path.split(netrefsub_snelder)[1],
                                                        re.sub('[.]', '', str(dis[0])),
                                                        re.sub('[.]', '', str(dis[1]))))
    sqlexp = 'rhtvs2_all_phi_qclass_MODULE >= {0} AND rhtvs2_all_phi_qclass_MODULE <= {1}'.format(dis[0], dis[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(netrefsub_snelder, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

#Create a subselection of HydroSHEDS river sections that overlap with the French dataset
arcpy.MakeFeatureLayer_management(hydrobasin12, 'hydrofr')
arcpy.SelectLayerByLocation_management('hydrofr', overlap_type='CONTAINS', select_features=outnet_snelder,
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('hydrofr', bassel_snelder)


#-------------------- Brazil ----------------------------------------------------------------------------------------------
#ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/shapefile/
#ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/informacoes_tecnicas/
"ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/shapefile/bcim_2016_shapefiles_21-11-2018.zip"
"ftp://geoftp.ibge.gov.br/cartas_e_mapas/mapa_indice_digital_4ed/produto_mapa_indice_digital/" #For topographic map footprints
#Select streams based on IMPRESSA_250000 in Shap

#Download and unzip database
url_bz = "ftp://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bcim/versao2016/shapefile/bcim_2016_shapefiles_21-11-2018.zip"
outzip_bz = os.path.join(datdir_bz, 'bcim_2016_shapefiles_21-11-2018.zip')
if not arcpy.Exists(outzip_bz ):
    print(url_bz)
    dlfile(url=url_bz, outpath=outzip_bz, ignore_downloadable=True)
    unzip(outzip_bz)

netref_bz = os.path.join(datdir_bz, 'hid_trecho_drenagem_l.shp')
netproj_bz = os.path.join(resdir_bz, 'netrawproj')
danglep_bz = os.path.join(resdir_bz, 'danglep')
netmain_bz = os.path.join(resdir_bz, 'netmain')
netdangle_bz = os.path.join(resdir_bz, 'netdangle')

#Unfortunately, network is poorly formatted. Lines are not digitized in correct direction and there are scores
#of topological issues (many lines that should be connected are > 100 m apart)
arcpy.Project_management(netref_bz, netproj_bz, out_coor_system='5880')
arcpy.FeatureVerticesToPoints_management(netproj_bz, danglep_bz, point_location='DANGLE')
danglelist = {row[0] for row in arcpy.da.UpdateCursor(danglep_bz, ['ORIG_FID'])}
danglef = 'dangle'
if danglef not in [f.name for f in arcpy.ListFields(netproj_bz)]:
    arcpy.AddField_management(netproj_bz, danglef, field_type='SHORT')

with arcpy.da.UpdateCursor(netproj_bz, ['OBJECTID', danglef]) as cursor:
    for row in cursor:
        if row[0] in danglelist:
            if row[1] is None:
                row[1] = 1
            elif row[1] > 0:
                row[1] += 1
            cursor.updateRow(row)

arcpy.MakeFeatureLayer_management(netproj_bz, out_layer='subriver',
                                  where_clause='dangle IS NULL')
arcpy.CopyFeatures_management('subriver', netmain_bz)

arcpy.MakeFeatureLayer_management(netproj_bz, out_layer='subriver',
                                  where_clause='dangle = 1')
arcpy.CopyFeatures_management('subriver', netdangle_bz)
#Also cheked 1:250K dataset and provides nothing more
#https://geoftp.ibge.gov.br/cartas_e_mapas/bases_cartograficas_continuas/bc250/versao2019/geopackage/

#-------------------- Australia -----------------------------------
#Download and unzip database
url_au = "ftp://ftp.bom.gov.au/anon/home/geofabric/Geofabric_Metadata_GDB_V3_2.zip"
outzip_au = os.path.join(datdir_au, 'fabric_Metadata_GDB_V3_2.zip')
if not arcpy.Exists(outzip_au):
    print(url_au)
    dlfile(url=url_au, outpath=outzip_au, ignore_downloadable=True)
    unzip(outzip_au)

#Geofabric_National_V3_2_PRODUCT_README.txt
arcpy.Describe(netref_au).SpatialReference.name
[f.name for f in arcpy.ListFields(netref_au)]

#Export by drainage area size class (in m2)
for da in [[10, 100], [100, 1000], [10**3,10**4], [10**4,10**5], [10**5,10**6]]: #km2
    subriver_out = os.path.join(resdir_au,
                                '{0}_DA{1}_{2}'.format(os.path.split(netref_au)[1],
                                                        re.sub('[.]', '', str(da[0])),
                                                        re.sub('[.]', '', str(da[1]))))
    sqlexp = 'UpstrDArea >= {0} AND UpstrDArea <= {1}'.format(da[0]*10**6, da[1]*10**6) #Convert to m3
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(netref_au, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

arcpy.CopyFeatures_management(netref_au, os.path.join(resdir_au, os.path.split(netref_au)[1]))

#Select HyoSHEDS basins intersecting with Australian dataset
arcpy.Project_management( os.path.join(resdir_au, os.path.split(netref_au)[1]),
                          netproj_au, out_coor_system=hydrobasin12)
netbasinters_au = os.path.join(resdir_au, 'net_hydrobasins12')
arcpy.Intersect_analysis([netproj_au, hydrobasin12], out_feature_class=netbasinters_au, join_attributes='ALL')
arcpy.AddGeometryAttributes_management(netbasinters_au, Geometry_Properties='LENGTH_GEODESIC',
                                       Length_Unit='Kilometers')
arcpy.CopyRows_management(netbasinters_au, os.path.join(resdir, 'netbas12_inters_australia.csv'))


#-------------------- France - BD Topage® ----------------------------------------------------------------------------------------------
# bdtopage_url = http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDTopage/2019/BD_Topage_FXX_2019.zip
# bdtopo_url = ftp://BDTOPO_V3_ext:Aish3ho8!!!@ftp3.ign.fr/BDTOPO_3-0_2020-09-15/BDTOPO_3-0_HYDROGRAPHIE_SHP_LAMB93_FXX_2020-09-15.7z
# dlfile(url=url, outpath=os.path.join(datdir_us, "NHDPlus_hr"), ignore_downloadable=True)

#http://www.sandre.eaufrance.fr/?urn=urn:sandre:donnees:773::::::referentiel:3.1:html
#datdir_bdtopage = os.path.join(datdir_fr, 'bdtopage')
#innet_bdtopage = os.path.join(datdir_bdtopage, 'TronconHydrographique_FXX.shp')
#No network information. Would have to computer river orders


#-----------------------------------------------------------------------------------------------------------------------
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
# -------------------- FRANCE - ONDE EAU -------------------------------------------------------------------------------
for yr in range(2012, 2021):
    in_url = "https://onde.eaufrance.fr/sites/default/files/fichiers-telechargeables/onde_france_{}.zip".format(yr)
    dlfile(in_url,
           outpath=datdir_onde,
           outfile=os.path.split(in_url)[1],
           ignore_downloadable=True)

# Download Carthage Hydrographic network
dlfile(
    "http://services.sandre.eaufrance.fr/telechargement/geo/ETH/BDCarthage/FXX/2017/France_metropole_entiere/SHP/TronconHydrograElt_FXX-shp.zip",
    outpath=datdir_onde,
    outfile="TronconHydrograElt_FXX-shp.zip",
    ignore_downloadable=True)

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

#######################################################################################################################
# http://www.freshwaterplatform.eu/
# http://project.freshwaterbiodiversity.eu/
# http://www.lifetrivers.eu/
# http://irbas.inrae.fr/people/related-publications
# http://www.ub.edu/fem/index.php/en/inici-riunet-en
# https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e