from utility_functions import *

arcpy.env.overwriteOutput= True
arcpy.env.qualifiedFieldNames = 'False'

# Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')
figsharedir = os.path.join(resdir, 'figshare')

#Input vars
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
riveratlas = os.path.join(outgdb, 'RiverATLAS_v10pred')

#Output var
gpreds1 = os.path.join(resdir, 'GRDCstations_predbasic800.shp') #main.GRDCstations_predbasic800
gpreds30 = os.path.join(resdir, 'GRDCstations_predbasic800_mdur30.shp') #main.GRDCstations_predbasic800_mdur30

######################## ANALYSIS ################################################
#Rename the fields for predictions (this step was added for sharing the network)
arcpy.AlterField_management(in_table = riveratlas, field='predbasic800', new_field_name='predprob1')
arcpy.AlterField_management(in_table = riveratlas, field='predbasic800_mdur30', new_field_name='predprob30')

#Export to figshare dir
riveratlas_figsharegdb = os.path.join(figsharedir, 'GIRES_v10.gdb', 'GIRES_v10_rivers')
arcpy.CopyFeatures_management(riveratlas, riveratlas_figsharegdb)
#Re-shuffle and format columns in ArcGIS (too cumbersome, ran out of time to write it here).
#Rename fields for export to .shp

#Divide into continents
pfaf1 = {1:'af', 8:'ar', 4:'as', 5:'au', 2:'eu', 9:'gr', 7:'na', 6:'sa', 3:'si'}
arcpy.AddField_management(riveratlas_figsharegdb, 'continent', 'text')
with arcpy.da.UpdateCursor(riveratlas_figsharegdb, ['continent', 'PFAF_ID05']) as cursor:
    for row in cursor:
        row[0] = pfaf1[int(str(row[1])[0])]
        cursor.updateRow(row)

#Export to shapefile
arcpy.SplitByAttributes_analysis(riveratlas_figsharegdb,
                                 Target_Workspace= os.path.join(figsharedir, 'GIRES_v10_shp'),
                                 Split_Fields='continent')
arcpy.DeleteField_management(riveratlas_figsharegdb, 'continent')

#Get metadata for river network
metafield(riveratlas_figsharegdb, os.path.join(figsharedir, 'riveratlas_IRESpred_metadata.csv'))

#Re-export gauges in shapefile in ArcGIS

#Join the pred1 and pred30
arcpy.AddField_management(gpreds1, 'predprob30', 'FLOAT')
arcpy.AddField_management(gpreds1, 'predcat30', 'SHORT')
arcpy.AddField_management(gpreds1, 'predres30', 'FLOAT')
pred30dict = {row[0]:[row[1], row[2], row[3]] for row in
              arcpy.da.SearchCursor(gpreds30, ['GAUGE_NO', 'IRpredprob', 'IRpredcat_', 'preduncert'])}
with arcpy.da.UpdateCursor(gpreds1, ['GAUGE_NO', 'predprob30', 'predcat30', 'predres30']) as cursor:
    for row in cursor:
        row[1] = pred30dict[row[0]][0]
        row[2] = pred30dict[row[0]][1]
        row[3] = pred30dict[row[0]][2]
        cursor.updateRow(row)

gtemplatepath = os.path.join(figsharedir, 'gauging_stations_template.csv')
if not os.path.exists(gtemplatepath):
    metafield(gpreds1, gtemplatepath)
gtemplate = pd.read_csv(gtemplatepath)
gtemplatedict = OrderedDict()
for row in gtemplate.iterrows():
    gtemplatedict[row[1]['fname']] = row[1]['fname_new']

#Create FieldMappings to subset, rename, and order dataset
gflist = [f.name for f in arcpy.ListFields(gpreds1)]
fms_gauges = arcpy.FieldMappings()
for forig in gtemplatedict:
    if forig in gflist and forig != "Shape":
        fnew = gtemplatedict[forig]
        print("{0} : {1}").format(forig, fnew)
        fm = arcpy.FieldMap()
        fm.addInputField(gpreds1, forig)
        of = fm.outputField
        of.name = fnew
        of.aliasName = fnew
        fm.outputField = of
        fms_gauges.addFieldMap(fm)

arcpy.FeatureClassToFeatureClass_conversion(gpreds1, out_path=os.path.join(figsharedir, 'GIRES_v10.gdb'),
                                            out_name='GIRES_v10_stations',
                                            field_mapping=fms_gauges)
arcpy.FeatureClassToFeatureClass_conversion(gpreds1, out_path=os.path.join(figsharedir, 'GIRES_v10_shp'),
                                            out_name='GIRES_v10_stations.shp',
                                            field_mapping=fms_gauges)

#Get metadata for gauges
gmeta =  os.path.join(figsharedir, 'gauging_stations_metadata.csv')
if not os.path.exists(gmeta):
    metafield(os.path.join(figsharedir, 'GIRES_v10_shp', 'GIRES_v10_stations'),
              gmeta)