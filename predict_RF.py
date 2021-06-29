import arcpy
import os
import pandas as pd
from utility_functions import *

arcpy.env.overwriteOutput= True
arcpy.env.qualifiedFieldNames = 'False'

# Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')


#Input var
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
riveratlas_orig = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
riveratlas_v11 = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v11.gdb', 'RiverATLAS_v11')
riveratlas_predtablist = getfilelist(resdir, 'RiverATLAS_predbasic800_[0-9]{8}[.]csv$')
riveratlas_predtab = pd.read_csv(riveratlas_predtablist[-1])
riveratlas_predtabmdur30list = getfilelist(resdir, 'RiverATLAS_predbasic800_mdur30_[0-9]{8}[.]csv$')
riveratlas_predtabmdur30 = pd.read_csv(riveratlas_predtabmdur30list[-1])

predvars_tab = pd.read_csv(os.path.join(resdir, 'predictor_variables.csv'))

insitu_resdir = os.path.join(resdir, 'Insitu_databases')
ondeobs_ipr = getfilelist(insitu_resdir, 'ondeobs_IPR_predbasic800cat[0-9]+[.]shp$')[-1]
pnwobs_ipr = getfilelist(insitu_resdir, 'pnwobs_IPR_predbasic800cat[0-9]+[.]shp$')[-1]

#Output vars
riveratlas_v10format = os.path.join(outgdb, 'RiverATLAS_v10_o0selcols')
riveratlas_v11format = os.path.join(outgdb, 'RiverATLAS_v11_o0selcols')
riveratlas = os.path.join(outgdb, 'RiverATLAS_v10pred')
ondeobs_ipr_netjoin = os.path.join(insitu_resdir, 'ondeobs_IPRnetjoin.shp')
pnwobs_ipr_netjoin = os.path.join(insitu_resdir, 'pnwobs_IPRnetjoin.shp')

#Set variables to include in output network
varlist = ["HYRIV_ID", "NEXT_DOWN", "LENGTH_KM", "CATCH_SKM", "UPLAND_SKM",
        "ORD_STRA", "INLAKEPERC", "HYBAS_ID", "PFAF_ID05", "MID_X", "MID_Y"]
varlist.extend(list(predvars_tab.sort_values(by=['Category', 'Source', 'varcode'])['varcode']))

#Creatge lists of fields in original River_ATLAS (v1.0) and new variables in River_ATLAS (v1.1, and format these)
f10l = [f.name for f in arcpy.ListFields(riveratlas_orig)]

f11l = OrderedDict()
for f in arcpy.ListFields(riveratlas_v11):
    if f.name[-3:] == '_11':
        f11l[f.name[:-3]] = f.name
    else:
        f11l[f.name] = f.name
f11l['HYRIV_ID'] = 'REACH_ID'

#Create FieldMappings for v10 to make sure it is correctly ordered.
fms_v10 = arcpy.FieldMappings()
for fsel in varlist:
    if (fsel in f10l and fsel not in f11l) or fsel == 'HYRIV_ID':
        #print(fsel)
        fm = arcpy.FieldMap()
        fm.addInputField(riveratlas_orig, fsel)
        fms_v10.addFieldMap(fm)
del fsel
del fm

#Create FieldMappings for v11 to make sure it is correctly ordered and named.
fms_v11 = arcpy.FieldMappings()
for fsel in varlist:
    if fsel in f11l:
        print("{0} : {1}".format(f11l[fsel], fsel))
        fm = arcpy.FieldMap()
        fm.addInputField(riveratlas_v11, f11l[fsel])
        of = fm.outputField
        of.name = fsel
        of.aliasName = fsel
        fm.outputField = of
        fms_v11.addFieldMap(fm)

# Create copy of RiverATLAS for display
if not arcpy.Exists(riveratlas_v10format):
    arcpy.FeatureClassToFeatureClass_conversion(
        in_features=riveratlas_orig,
        out_path=os.path.split(riveratlas_v10format)[0],
        out_name=os.path.split(riveratlas_v10format)[1],
        where_clause='dis_m3_pyr > 0',
        field_mapping = fms_v10
    )

if not arcpy.Exists(riveratlas_v11format):
    arcpy.FeatureClassToFeatureClass_conversion(
        in_features=riveratlas_v11,
        out_path=os.path.split(riveratlas_v11format)[0],
        out_name=os.path.split(riveratlas_v11format)[1],
        field_mapping = fms_v11
    )

#Join v10 and v11
arcpy.MakeFeatureLayer_management(riveratlas_v10format,  'riv10format')
arcpy.AddJoin_management(in_layer_or_view='riv10format', in_field='HYRIV_ID',
                         join_table=riveratlas_v11format, join_field='HYRIV_ID',
                         join_type='KEEP_COMMON')
arcpy.CopyFeatures_management('riv10format', riveratlas)

#Delete superfluous fields and formatted separate networks
for fdel in ['OBJECTID_1', 'HYRIV_ID_1', 'UPLAND_SKM_1', 'HYBAS_ID03']:
    if fdel in [f.name for f in arcpy.ListFields(riveratlas)]:
        arcpy.DeleteField_management(riveratlas, fdel)
arcpy.Delete_management(riveratlas_v10format)
arcpy.Delete_management(riveratlas_v11format)

# Add predictions (join is to slow)
if not 'predbasic800' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800', 'FLOAT')
    arcpy.AddField_management(riveratlas, 'predcat1', 'SHORT')
if not 'predbasic800_mdur30' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800_mdur30', 'FLOAT')
    arcpy.AddField_management(riveratlas, 'predcat30', 'SHORT')

# .astype({'predbasic800cat': int})  # Convert intermittent predicted class to integer
preddict = (riveratlas_predtab[~riveratlas_predtab['predbasic800'].isnull()][['HYRIV_ID', 'predbasic800']]  # only keep records that have prediction
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

preddict_mdur30 = (
    riveratlas_predtabmdur30[~riveratlas_predtabmdur30['predbasic800_mdur30'].isnull()][['HYRIV_ID', 'predbasic800_mdur30']]  # only keep records that have prediction
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'predbasic800', 'predcat1',
                                        'predbasic800_mdur30', 'predcat30']) as cursor: #'predbasic800_mdur30'
    x = 0
    for row in cursor:
        if x % 100000 == 0:
            print(x)
        try:
            row[1] = preddict[row[0]]
            row[2] = int(preddict[row[0]] >= 0.5)
            row[3] = preddict_mdur30[row[0]]
            row[4] = int(preddict_mdur30[row[0]] >= 0.5)
            cursor.updateRow(row)
        except:
            pass
        x += 1

# [f.name for f in arcpy.ListFields(riveratlas)]
# ord = {row[0] for row in arcpy.da.SearchCursor(riveratlas, 'ORD_STRA')}

#Export predictions by streamflow size class
for dis in [[0,0.1], [0.1,1], [1,10], [10,100], [100,1000], [1000, 10000], [10000, 1000000]]:
    subriver_out = os.path.join(os.path.split(riveratlas)[0],
                                '{0}_DIS{1}_{2}'.format(os.path.split(riveratlas)[1],
                                                        re.sub('[.]', '', str(dis[0])),
                                                        re.sub('[.]', '', str(dis[1]))))
    sqlexp = 'dis_m3_pyr >= {0} AND dis_m3_pyr <= {1}'.format(dis[0], dis[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(riveratlas, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')

#Export predictions by drainage area size class
for da in [[0, 10], [10, 100], [100, 1000], [10**3,10**4], [10**4,10**5], [10**5,10**6], [10**6,10**7]]:
    subriver_out = os.path.join(os.path.split(riveratlas)[0],
                                '{0}_DA{1}_{2}'.format(os.path.split(riveratlas)[1],
                                                        re.sub('[.]', '', str(da[0])),
                                                        re.sub('[.]', '', str(da[1]))))
    sqlexp = 'UPLAND_SKM >= {0} AND UPLAND_SKM <= {1}'.format(da[0], da[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(riveratlas, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')


#Join ONDE observations + preds to network
print('Joinin {} to network'.format(ondeobs_ipr))
arcpy.MakeFeatureLayer_management(riveratlas, 'riveratlaslyr')
arcpy.AddJoin_management('riveratlaslyr', in_field='HYRIV_ID',
                         join_table=ondeobs_ipr, join_field = 'HYRIVID',
                         join_type='KEEP_COMMON')
arcpy.CopyFeatures_management('riveratlaslyr', out_feature_class=ondeobs_ipr_netjoin)
arcpy.Delete_management('riveratlaslyr')

#Join PNW observations + preds to network
print('Joinin {} to network'.format(pnwobs_ipr))
arcpy.MakeFeatureLayer_management(riveratlas, 'riveratlaslyr')
arcpy.AddJoin_management('riveratlaslyr', in_field='HYRIV_ID',
                         join_table=pnwobs_ipr, join_field = 'HYRIV_I',
                         join_type='KEEP_COMMON')
arcpy.CopyFeatures_management('riveratlaslyr', out_feature_class=pnwobs_ipr_netjoin)
arcpy.Delete_management('riveratlaslyr')