"""
Download national river networks, model outputs, and on-the-ground observation of river intermittence
"""

from setup_localIRformatting import *
import math

#US dirs
datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

inbas_us = os.path.join(datdir_us, 'WBD_National_GDB.gdb', 'WBD', 'WBDHU8')
outnet_us = os.path.join(resdir_us, 'network')
bassel_us = os.path.join(resdir_us, 'hydrobasins12')
refnetsub_us = os.path.join(resdir_us, 'refnet_o10')

#France dirs
datdir_fr = os.path.join(compdatdir, 'France')
pathcheckcreate(datdir_fr)
resdir_fr = os.path.join(compresdir, 'france.gdb')
pathcheckcreate(resdir_fr)
outnet_fr = os.path.join(resdir_fr, 'network')
bassel_fr = os.path.join(resdir_fr, 'hydrobasins12')
refnetsub_fr = os.path.join(resdir_fr, 'refnet_o10')

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
arcpy.Merge_management(inputs=flsubl, output=refnetsub_us)

#Divide by river order
for ord in [[1,3], [4,5], [6,7], [8,9], [10, 11]]:
    subriver_out = os.path.join(resdir_us,
                                '{0}_ORD{1}_{2}'.format(os.path.split(refnetsub_us)[1], ord[0], ord[1]))
    if not arcpy.Exists(subriver_out):
        sqlexp = 'StreamOrde >= {0} AND StreamOrde <= {1}'.format(ord[0], ord[1])
        print(sqlexp)
        arcpy.MakeFeatureLayer_management(refnetsub_us, out_layer='subriver', where_clause=sqlexp)
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

for ord in [[1,3], [4,5], [6,7], [8,9], [10, 11]]:
    subriver_out = os.path.join(resdir_us,
                                '{0}_ORD{1}_{2}'.format(os.path.split(outflsubmr)[1], ord[0], ord[1]))
    if not arcpy.Exists(subriver_out):
        sqlexp = 'StreamOrde >= {0} AND StreamOrde <= {1}'.format(ord[0], ord[1])
        print(sqlexp)
        arcpy.MakeFeatureLayer_management(outflsubmr, out_layer='subriver', where_clause=sqlexp)
        arcpy.CopyFeatures_management('subriver', subriver_out)
        arcpy.Delete_management('subriver')

#-------------------- France ----------------------------------------------------------------------------------------------
#Data are personal communications from Snelder et al. 2013

#Join original French river network data with Snelder et al.' predictions
arcpy.MakeFeatureLayer_management(os.path.join(datdir_fr, 'rhtvs2_all_phi_qclass.shp'),
                                  out_layer='frlyr')
arcpy.AddJoin_management(in_layer_or_view='frlyr', in_field='ID_DRAIN',
                         join_table=os.path.join(datdir_fr, 'INT_RF.txt'), join_field='AllPred$ID_DRAIN')
arcpy.CopyFeatures_management(in_features='frlyr', out_feature_class=outnet_fr)

arcpy.DefineProjection_management(in_dataset=outnet_fr, coor_system=arcpy.SpatialReference(2192)) #ED_1950_France_EuroLambert

#Subsetlect all segments with drainage area > 10 km2
arcpy.MakeFeatureLayer_management(outnet_fr,
                                  out_layer='frlyr', where_clause='rhtvs2_all_phi_qclass_SURF_BV >= 10')
arcpy.CopyFeatures_management(in_features='frlyr', out_feature_class=refnetsub_fr)

#Divide french network by river order
for ord in [[1,3], [4,5], [6,7], [8,9]]:
    subriver_out = os.path.join(resdir_fr,
                                '{0}_ORD{1}_{2}'.format(os.path.split(refnetsub_fr)[1], ord[0], ord[1]))
    sqlexp = 'rhtvs2_all_phi_qclass_STRAHLER >= {0} AND rhtvs2_all_phi_qclass_STRAHLER <= {1}'.format(ord[0], ord[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(refnetsub_fr, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

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