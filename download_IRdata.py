"""
Download national river networks, model outputs, and on-the-ground observation of river intermittence
"""

from utility_functions import *

arcpy.env.qualifiedFieldNames = 'False'

compdatdir = os.path.join(datdir, 'Comparison_databases')
compresdir = os.path.join(resdir, 'Comparison_databases')

pathcheckcreate(compdatdir)
pathcheckcreate(compresdir)

hydrobasin12 = os.path.join(datdir, 'HydroATLAS', 'BasinATLAS_v10.gdb', 'BasinATLAS_v10_lev12')

#-------------------- USA ----------------------------------------------------------------------------------------------
#Extra links
#"https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/National/HighResolution/GDB/NHD_H_National_GDB.zip"
#"https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Hydrography/NHDPlusHR/Beta/GDB/"

datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

outnet_us = os.path.join(resdir_us, 'network')
bassel_us = os.path.join(resdir_us, 'hydrobasins12')

#Download NHDplus data
nhdplus_linklist = ["https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHDPlusHR/Beta/" \
                    "GDB/NHDPLUS_H_{}_HU4_GDB.zip".format(str(i).zfill(4)) for i in xrange(832, 2102)]
for url in nhdplus_linklist:
    print(url)
    dlfile(url=url, outpath=datdir_us, ignore_downloadable=True)

#Merge all flow lines
FLlist = [os.path.join(i, 'NHDFlowLine') for i in 
          getfilelist(dir=datdir_us, repattern='Hydrography', gdbf=True, nongdbf=False)]

arcpy.Merge_management(FLlist, output=outnet_us)

#Add attributes
#Drainage area
fldl_us = OrderedDict()
fldl_us['NHDPlusID']="LONG"
fldl_us['StreamOrde']="SHORT"
fldl_us['TotDASqKm']="FLOAT"
fldl_us['QAMA']='FLOAT'

VAlist = getfilelist(dir=datdir_us, repattern='NHDPlusFlowlineVAA', gdbf=True, nongdbf=False)
attridict_us = defaultdict(list)
for tab in VAlist[0:10]:
    print(tab)
    for row in arcpy.da.SearchCursor(tab, ['NHDPlusID', 'StreamOrde', 'TotDASqKm']):
        attridict_us[row[0]] = [row[1], row[2]]

#Discharge
EROMQAlist = getfilelist(dir=datdir_us, repattern='NHDPlusEROMMA', gdbf=True, nongdbf=False)
for tab in EROMQAlist[0:10]:
    print(tab)
    for row in arcpy.da.SearchCursor(tab, ['NHDPlusID', 'QAMA']):
        if row[0] in attridict_us:
            attridict_us[row[0]].append(row[1])

#Join attributes to merged dataset
for fl in fldl_us:
    arcpy.AddField_management(outnet_us, field_name=fl, field_type=fldl_us[fl])

with arcpy.da.UpdateCursor(outnet_us, fldl_us.keys()) as cursor:
    x=0
    for row in cursor:
        if x % 100000 == 0:
            print(x)

        if row[0] in attridict_us:
            row[1] = attridict_us[row[0]][0] #StreamOrder for NHDPlusID
            row[2] = attridict_us[row[0]][1] #TotDASqKm
            row[3] = attridict_us[row[0]][2] #QAMA

            cursor.updateRow(row)

        x += 1

#Create a subselection of HydroSHEDS river sections that overlap with the French dataset
arcpy.MakeFeatureLayer_management(hydrobasin12, 'hydrous')
arcpy.SelectLayerByLocation_management('hydrous', overlap_type='CONTAINS', select_features=outnet_us,
                                       selection_type='NEW_SELECTION')
arcpy.CopyFeatures_management('hydrous', bassel_us)

#-------------------- France ----------------------------------------------------------------------------------------------
#Data are personal communications from Snelder et al. 2013
datdir_fr = os.path.join(compdatdir, 'France')
pathcheckcreate(datdir_fr)
resdir_fr = os.path.join(compresdir, 'france.gdb')
pathcheckcreate(resdir_fr)
outnet_fr = os.path.join(resdir_fr, 'network')
bassel_fr = os.path.join(resdir_fr, 'hydrobasins12')

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

http://www.freshwaterplatform.eu/
http://project.freshwaterbiodiversity.eu/
http://www.lifetrivers.eu/
http://irbas.inrae.fr/people/related-publications
http://www.ub.edu/fem/index.php/en/inici-riunet-en
https://onde.eaufrance.fr/content/t%C3%A9l%C3%A9charger-les-donn%C3%A9es-des-campagnes-par-ann%C3%A9e