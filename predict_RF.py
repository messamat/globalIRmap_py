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
riveratlas_predtablist = getfilelist(resdir, 'RiverATLAS_predbasic800_[0-9]{8}[.]csv$')
riveratlas_predtab = pd.read_csv(riveratlas_predtablist[-1])
riveratlas_predtabmdur30list = getfilelist(resdir, 'RiverATLAS_predbasic800_mdur30_[0-9]{8}[.]csv$')
riveratlas_predtabmdur30 = pd.read_csv(riveratlas_predtabmdur30list[-1])

insitu_resdir = os.path.join(resdir, 'Insitu_databases')
ondeobs_ipr = getfilelist(insitu_resdir, 'ondeobs_IPR_predbasic800cat[0-9]+[.]shp$')[-1]
pnwobs_ipr = getfilelist(insitu_resdir, 'pnwobs_IPR_predbasic800cat[0-9]+[.]shp$')[-1]

#Output vars
riveratlas = os.path.join(outgdb, 'RiverATLAS_v10pred')
ondeobs_ipr_netjoin = os.path.join(insitu_resdir, 'ondeobs_IPRnetjoin.shp')
pnwobs_ipr_netjoin = os.path.join(insitu_resdir, 'pnwobs_IPRnetjoin.shp')


# Create copy of RiverATLAS for display
if not arcpy.Exists(riveratlas):
    arcpy.MakeFeatureLayer_management(riveratlas_orig, 'riveratlas_origlyr')
    arcpy.SelectLayerByAttribute_management('riveratlas_origlyr', selection_type='NEW_SELECTION',
                                            where_clause='dis_m3_pyr > 0')
    arcpy.CopyFeatures_management('riveratlas_origlyr', riveratlas)

# Add predictions (join is to slow)
if not 'predbasic800' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800', 'FLOAT')
if not 'predbasic800_mdur30' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800_mdur30', 'FLOAT')

# .astype({'predbasic800cat': int})  # Convert intermittent predicted class to integer
preddict = (riveratlas_predtab[~riveratlas_predtab['predbasic800'].isna()][['HYRIV_ID', 'predbasic800']]  # only keep records that have prediction
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

preddict_mdur30 = (
    riveratlas_predtabmdur30[~riveratlas_predtabmdur30['predbasic800_mdur30'].isna()][['HYRIV_ID', 'predbasic800_mdur30']]  # only keep records that have prediction
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'predbasic800', 'predbasic800_mdur30']) as cursor:
    x = 0
    for row in cursor:
        if x % 100000 == 0:
            print(x)
        try:
            row[1] = preddict[row[0]]
            row[2] = preddict_mdur30[row[0]]
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