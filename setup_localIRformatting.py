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
hydroacc = os.path.join(datdir, 'Bernhard\HydroATLAS\HydroATLAS_Geometry\Accu_area_grids\upstream_area_skm_15s.gdb',
                        'up_area_skm_15s')

#USwide dirs
datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

inbas_us = os.path.join(datdir_us, 'WBD_National_GDB.gdb', 'WBD', 'WBDHU8')
