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

hydrodir = os.path.join(datdir, 'Bernhard\HydroATLAS')
hydrogeomdir = os.path.join(hydrodir, 'HydroATLAS_Geometry')
hydrobasin12 = os.path.join(hydrodir, 'HydroATLAS_v10_final_data\\', 'BasinATLAS_v10_shp\\BasinATLAS_v10_lev12.shp')
riveratlas = os.path.join(hydrodir, 'HydroATLAS_v10_final_data\\', 'RiverATLAS_v10.gdb', 'RiverATLAS_v10')
hydroacc = os.path.join(hydrogeomdir, 'Accu_area_grids\upstream_area_skm_15s.gdb', 'up_area_skm_15s')
hydrolink = os.path.join(hydrogeomdir, 'Link_zone_grids', 'link_stream.gdb', 'link_str_arc')
hydrodisras = os.path.join(hydrodir, 'HydroATLAS_Data\Hydrology\discharge_wg22_1971_2000.gdb',
                                  'dis_nat_wg22_ls_year')

#USwide dirs
datdir_us = os.path.join(compdatdir, 'US')
pathcheckcreate(datdir_us)
resdir_us = os.path.join(compresdir, 'us.gdb')
pathcheckcreate(resdir_us)

inbas_us = os.path.join(datdir_us, 'WBD_National_GDB.gdb', 'WBD', 'WBDHU8')
