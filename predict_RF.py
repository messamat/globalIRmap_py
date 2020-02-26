import arcpy
import os
import pandas as pd

# Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')
resdir = os.path.join(rootdir, 'results')

#
outgdb = os.path.join(resdir, 'spatialoutputs.gdb')
riveratlas_orig = os.path.join(datdir, 'HydroATLAS', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
riveratlas_predtab = pd.read_csv(os.path.join(resdir, 'RiverATLAS_predbasic800.csv'))
riveratlas = os.path.join(outgdb, 'RiverATLAS_v10')

# Create copy of RiverATLAS for display
if not arcpy.Exists(riveratlas):
    arcpy.CopyFeatures_management(riveratlas_orig, riveratlas)

# Add predictions (join is to slow)
if not 'predbasic800cat' in [f.name for f in arcpy.ListFields(riveratlas)]:
    arcpy.AddField_management(riveratlas, 'predbasic800cat', 'SHORT')

preddict = (riveratlas_predtab[~riveratlas_predtab['predbasic800cat'].isna()]  # only keep records that have prediction
            .astype({'predbasic800cat': int})  # Convert intermittent predicted class to integer
            .set_index('HYRIV_ID')  # Set index
            .squeeze()  # Convert to panda series
            .to_dict())  # Convert to dictionary for fast access

with arcpy.da.UpdateCursor(riveratlas, ['HYRIV_ID', 'predbasic800cat']) as cursor:
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
