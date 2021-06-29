"""Developer: Mathis Messager
Purpose:
 - defines folder structure for formatting data to compare modeled estimates of global flow intermittence to national
   hydrographic datasets (Comparison_databases) and to in-situ/field-based observations of flow intermittence.
 - defines functions used in formatting data for the comparisons
"""

from utility_functions import *

rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
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

#Create panda dataframe from ArcGIS table
def get_arctab_df(in_tab, IDfield='OID@', valuefields='*', resetid = False):
    if IDfield == 'OID@':
        IDfield = arcpy.Describe(in_tab).OIDFieldName
    if valuefields == '*':
        valuefields = [f.name for f in arcpy.ListFields(in_tab)]
        valuefields.remove(IDfield)

    if isinstance(valuefields, list):
        fl = valuefields
        fl.insert(0, IDfield)
    else:
        fl = [IDfield, valuefields]

    ddict = {row[0]: row[1:] for row in arcpy.da.SearchCursor(in_tab, fl)}

    df = pd.DataFrame.from_dict(ddict, 'index')
    if resetid:
        df.reset_index(inplace=True)
    else:
        fl.remove(IDfield)

    df.columns = fl

    return(df)

#Compute average bearing for each line
def projdistbearpt(x1, y1, x2, y2):
    # Distance between the two points
    dist = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    # Angle (counterclockwise)
    if y2-y1 > 0:
        angle = ((math.atan((x2 - x1) / (y2 - y1)) * 180 / math.pi) - 90)
    else:
        angle = -90

    # Revert direction (clockwise) to get bearing from north
    if angle > 0:
        bearing = 360 - angle
    else:
        bearing = -angle

    return ([dist, bearing])

#Append XY position of all vertices for each line
def linevert_wbearing(in_vert, in_IDfield):
    vertazdict = defaultdict(list)
    with arcpy.da.SearchCursor(in_vert, [in_IDfield, 'SHAPE@XY']) as cursor:
        for row in cursor:
            vertazdict[row[0]].append(row[1])

    #Replace vertazdict values with average bearing (out of 360)
    for k, v in vertazdict.iteritems():
        wanglesum = 0
        distsum = 0
        for v1, v2 in zip(v[:-1], v[1:]):
            distbear = projdistbearpt(v1[0], v1[1], v2[0], v2[1])
            wanglesum += distbear[0] * distbear[1]
            distsum += distbear[0]
        vertazdict[k] = wanglesum/float(distsum)

    return(vertazdict)

#Snap to specific features
def indivsnap(in_features, in_field, join_field, indivsnap_env):
    arcpy.MakeFeatureLayer_management(in_features, 'obslyr')
    arcpy.MakeFeatureLayer_management(indivsnap_env[0][0], 'netlyr')

    oidfn = arcpy.Describe(in_features).OIDFieldName

    try:
        with arcpy.da.SearchCursor(in_features, [oidfn, in_field]) as cursor:
            for row in cursor:
                print(row[0])
                arcpy.SelectLayerByAttribute_management('obslyr', 'NEW_SELECTION', '{0} = {1}'.format(oidfn, row[0]))
                arcpy.SelectLayerByAttribute_management('netlyr', 'NEW_SELECTION', '{0} = {1}'.format(join_field, row[1]))
                indivsnap_env[0][0] = 'netlyr'

                arcpy.Snap_edit('obslyr', indivsnap_env)
    except Exception:
        print("Exception in user code:")
        traceback.print_exc(file=sys.stdout)
        arcpy.Delete_management('obslyr')
        arcpy.Delete_management('netlyr')

def expand_extent(in_extent, dist):
    return (
        arcpy.Extent(in_extent.XMin - dist,
                     in_extent.YMin - dist,
                     in_extent.XMax + dist,
                     in_extent.YMax + dist)
    )

def extract_obsdisacc(in_obs, in_net, obsin_field, netjoin_field, in_diss, in_acc, resdir, delintermediate=False):
    try:
        if not all(v in [f for f in arcpy.ListFields(in_obs)] for v in ['POINTdis', 'POINTDA']):
            arcpy.AddField_management(in_obs, 'POINTdis', 'DOUBLE')
            arcpy.AddField_management(in_obs, 'POINTDA', 'DOUBLE')

        oidfn = arcpy.Describe(in_obs).OIDFieldName

        arcpy.MakeFeatureLayer_management(in_net, 'netlyr')

        arcpy.env.snapRaster = in_acc
        refres = arcpy.Describe(in_acc).meanCellWidth
        refsr = arcpy.Describe(in_acc).SpatialReference

        templine = os.path.join(resdir, 'templine')
        temppt = os.path.join(resdir, 'temppt')
        tempras = os.path.join(resdir, 'tempnetras')
        subacc = os.path.join(resdir, 'subacc')
        subdis = os.path.join(resdir, 'subdis')

        # def get_accdis():
        with arcpy.da.SearchCursor(in_obs, [oidfn, obsin_field, 'SHAPE@XY', 'POINTdis']) as cursor:
            for row in cursor:
                outtab = os.path.join(resdir, 'stat{}'.format(row[0]))
                if (not arcpy.Exists(outtab)): #or (row[3] is None and arcpy.Exists(outtab)):
                    if row[1] is not None:
                        print(row[0])
                        arcpy.CopyFeatures_management(arcpy.PointGeometry(arcpy.Point(row[2][0], row[2][1])),
                                                      temppt)
                        arcpy.DefineProjection_management(temppt, refsr)

                        arcpy.SelectLayerByAttribute_management('netlyr', 'NEW_SELECTION',
                                                                '{0} = {1}'.format(netjoin_field, row[1]))
                        arcpy.CopyFeatures_management('netlyr', templine)

                        arcpy.env.extent = expand_extent(arcpy.Describe(templine).extent, dist=refres * 2)
                        arcpy.PolylineToRaster_conversion(templine,
                                                          value_field='HYRIV_ID',
                                                          out_rasterdataset=tempras,
                                                          cellsize=refres)

                        EucAllocation(Int(0.5 + 100 * Con(Raster(tempras) == row[1], in_acc)),
                                      maximum_distance=refres * 2).save(subacc)
                        EucAllocation(Int(0.5 + 100 * Con(Raster(tempras) == row[1], in_diss)),
                                      maximum_distance=refres * 2).save(subdis)

                        Sample(in_rasters=[subacc, subdis],
                               in_location_data=temppt,
                               out_table=outtab,
                               resampling_type='NEAREST'
                               )

                        arcpy.ClearEnvironment('extent')

        exvaldict = {}
        for tab in getfilelist(resdir, 'stat[0-9]+'):
            print(tab)
            for row in arcpy.da.SearchCursor(tab, ['subdis_Band_1', 'subacc_Band_1']):
                exvaldict[int(os.path.split(tab)[1][4:])] = [row[0], row[1]]

        with arcpy.da.UpdateCursor(in_obs, [oidfn, 'POINTdis', 'POINTDA']) as cursor:
            for row in cursor:
                if row[0] in exvaldict:
                    print(row[0])
                    row[1] = exvaldict[row[0]][0]
                    row[2] = exvaldict[row[0]][1]
                    cursor.updateRow(row)

        arcpy.ClearEnvironment('snapRaster')

        if delintermediate == True:
            print('Deleting all tables...')
            for tab in getfilelist(resdir, 'stat[0-9]+'):
                arcpy.Delete_management(tab)

            arcpy.Delete_management('netlyr')
            arcpy.Delete_management(templine)
            arcpy.Delete_management(temppt)
            arcpy.Delete_management(tempras)
            arcpy.Delete_management(subacc)
            arcpy.Delete_management(subdis)


    except:
        traceback.print_exc()
        arcpy.ClearEnvironment('snapRaster')
        arcpy.ClearEnvironment('extent')

        if arcpy.Exists('netlyr'):
            arcpy.Delete_management('netlyr')
        if arcpy.Exists(templine):
            arcpy.Delete_management(templine)
        if arcpy.Exists(temppt):
            arcpy.Delete_management(temppt)
        if arcpy.Exists(tempras):
            arcpy.Delete_management(tempras)
        if arcpy.Exists(subacc):
            arcpy.Delete_management(subacc)
        if arcpy.Exists(subdis):
            arcpy.Delete_management(subdis)
