import arcpy
import os
import pandas as pd

arcpy.env.overwriteOutput= True

# Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

#
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
riveratlas_orig = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
riveratlas_predtab = pd.read_csv(os.path.join(resdir, 'RiverATLAS_predbasic800.csv'))
riveratlas = os.path.join(outgdb, 'RiverATLAS_v10pred')

# Create copy of RiverATLAS for display
if not arcpy.Exists(riveratlas):
    arcpy.CopyFeatures_management(riveratlas_orig, riveratlas)

# Add predictions (join is to slow)
if not 'predbasic800' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800', 'FLOAT')

riveratlas_predtab.columns
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

for ord in [[1,2], [3,4], [5,6], [7,8], [9,10]]:
    subriver_out = os.path.join(os.path.split(riveratlas)[0],
                                '{0}_ORD{1}_{2}'.format(os.path.split(riveratlas)[1], ord[0], ord[1]))
    sqlexp = 'ORD_STRA >= {0} AND ORD_STRA <= {1}'.format(ord[0], ord[1])
    print(sqlexp)
    arcpy.MakeFeatureLayer_management(riveratlas, out_layer='subriver', where_clause=sqlexp)
    arcpy.CopyFeatures_management('subriver', subriver_out)
    arcpy.Delete_management('subriver')
