import arcpy
import os
import pandas as pd
from utility_functions import *

arcpy.env.overwriteOutput= True

# Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

#Input var
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
riveratlas_orig = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
riveratlas_predtab_u10o10 = os.path.join(resdir, 'RiverATLAS_predbasic800_20201027_u10o10.csv')
riveratlas_predtab_u10nodiso10 = os.path.join(resdir, 'RiverATLAS_predbasic800_20201027_u10nodiso10.csv')
riveratlas_predtab_u10o1 = os.path.join(resdir, 'RiverATLAS_predbasic800_20201028_u10o1.csv')

riveratlas = os.path.join(outgdb, 'RiverATLAS_v10pred')

riveratlas_predtab = pd.read_csv(riveratlas_predtab_u10o1)

# Create copy of RiverATLAS for display
if not arcpy.Exists(riveratlas):
    arcpy.MakeFeatureLayer_management(riveratlas_orig, 'riveratlas_origlyr')
    arcpy.SelectLayerByAttribute_management('riveratlas_origlyr', selection_type='NEW_SELECTION',
                                            where_clause='dis_m3_pyr > 0')
    arcpy.CopyFeatures_management('riveratlas_origlyr', riveratlas)

# Add predictions (join is to slow)
if not 'predbasic800' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800', 'FLOAT')

# .astype({'predbasic800cat': int})  # Convert intermittent predicted class to integer
preddict = (riveratlas_predtab[~riveratlas_predtab['predbasic800'].isna()][['HYRIV_ID', 'predbasic800']]  # only keep records that have prediction
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'predbasic800']) as cursor:
    x = 0
    for row in cursor:
        if x % 100000 == 0:
            print(x)
        try:
            row[1] = preddict[row[0]]
            cursor.updateRow(row)
        except:
            pass
        x += 1

# [f.name for f in arcpy.ListFields(riveratlas)]
# ord = {row[0] for row in arcpy.da.SearchCursor(riveratlas, 'ORD_STRA')}

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
