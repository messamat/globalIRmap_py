import os
import re
import requests
import pandas as pd
import zipfile
import gzip
import io
import csv
import arcpy
import itertools
import traceback
import urlparse
import sys
import ftplib
import numpy as np

def divbb(bbox, res, divratio):
    box_lc_x, box_lc_y, box_rc_x, box_rc_y = bbox
    coln = (box_rc_x - box_lc_x) / float(res)
    rown = (box_rc_y - box_lc_y) / float(res)

    xbblist = np.arange(box_lc_x, box_rc_x + (float(res) * coln / divratio),
                        float(res) * coln / divratio)
    ybblist = np.arange(box_lc_y, box_rc_y + (float(res) * rown / divratio),
                        float(res) * rown / divratio)

    if abs(xbblist[-1]) > abs(box_rc_x):
        xbblist[-1] = box_rc_x

    if abs(ybblist[-1]) > abs(box_rc_y):
        ybblist[-1] = box_rc_y

    xbblist = np.unique(xbblist)
    ybblist = np.unique(ybblist)

    fullbblist = []
    for pairx in zip(xbblist[:-1], xbblist[1:]):
        for pairy in zip(ybblist[:-1], ybblist[1:]):
            fullbblist.append((pairx[0], pairy[0], pairx[1], pairy[1]))

    return (fullbblist)

def getwkspfiles(dir, repattern):
    arcpy.env.workspace = dir
    filenames_list = (arcpy.ListDatasets() or []) + (
            arcpy.ListTables() or [])  # Either LisDatsets or ListTables may return None so need to create empty list alternative
    if not repattern == None:
        filenames_list = [os.path.join(dir, filen)
                          for filen in filenames_list if re.search(repattern, filen)]
    return (filenames_list)
    arcpy.ClearEnvironment('workspace')


def getfilelist(dir, repattern=None, gdbf=True, nongdbf=True):
    """Function to iteratively go through all subdirectories inside 'dir' path
    and retrieve path for each file that matches "repattern"
    If the provided path is an ArcGIS workspace"""
    try:
        if arcpy.Describe(dir).dataType == 'Workspace':
            if gdbf == True:
                print('{} is ArcGIS workspace...'.format(dir))
                getwkspfiles(dir, repattern)
            else:
                raise ValueError(
                    "A gdb workspace was given for dir but gdbf=False... either change dir or set gdbf to True")
        else:
            filenames_list = []

            if gdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for in_dir in dirnames:
                        fpath = os.path.join(dirpath, in_dir)
                        if arcpy.Describe(fpath).dataType == 'Workspace':
                            print('{} is ArcGIS workspace...'.format(fpath))
                            filenames_list.extend(getwkspfiles(dir=fpath, repattern=repattern))

            if nongdbf == True:
                for (dirpath, dirnames, filenames) in os.walk(dir):
                    for file in filenames:
                        if repattern is None:
                            filenames_list.append(os.path.join(dirpath, file))
                        else:
                            if re.search(repattern, file):
                                filenames_list.append(os.path.join(dirpath, file))
        return (filenames_list)

    # Return geoprocessing specific errors
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages(2))
    # Return any other type of error
    except:
        # By default any other errors will be caught here
        e = sys.exc_info()[1]
        print(e.args[0])


def pathcheckcreate(path, verbose=True):
    """"Function that takes a path as input and:
      1. Checks which directories and .gdb exist in the path
      2. Creates the ones that don't exist"""

    dirtocreate = []
    # Loop upstream through path to check which directories exist, adding those that don't exist to dirtocreate list
    while not os.path.exists(os.path.join(path)):
        dirtocreate.append(os.path.split(path)[1])
        path = os.path.split(path)[0]

    dirtocreate.reverse()

    # After reversing list, iterate through directories to create starting with the most upstream one
    for dir in dirtocreate:
        # If gdb doesn't exist yet, use arcpy method to create it and then stop the loop to prevent from trying to create anything inside it
        if os.path.splitext(dir)[1] == '.gdb':
            if verbose:
                print('Create {}...'.format(dir))
            arcpy.CreateFileGDB_management(out_folder_path=path,
                                           out_name=dir)
            break

        # Otherwise, if it is a directory name (no extension), make a new directory
        elif os.path.splitext(dir)[1] == '':
            if verbose:
                print('Create {}...'.format(dir))
            path = os.path.join(path, dir)
            os.mkdir(path)

# Function to download and unzip miscellaneous types of files
# Partly inspired from https://www.codementor.io/aviaryan/downloading-files-from-urls-in-python-77q3bs0un
def getfilelist(dir, repattern):
    return [os.path.join(dirpath, file)
            for (dirpath, dirnames, filenames) in os.walk(dir)
            for file in filenames if re.search(repattern, file)]


def mergedel(dir, repattern, outfile, delete=False, verbose=False):
    flist = getfilelist(dir, repattern)
    pd.concat([pd.read_csv(file, index_col=[0], parse_dates=[0])
               for file in flist],
              axis=0) \
        .sort_index() \
        .to_csv(outfile)
    print('Merged and written to {}'.format(outfile))

    if delete == True:
        for tab in flist:
            os.remove(tab)
            if verbose == True:
                print('Delete {}'.format(tab))


def is_downloadable(url):
    """
    Does the url contain a downloadable resource
    """
    h = requests.head(url, allow_redirects=True)
    header = h.headers
    content_type = header.get('content-type')
    if 'html' in content_type.lower():
        return False
    return True


def get_filename_from_cd(url):
    """
    Get filename from content-disposition
    """
    r = requests.get(url, allow_redirects=True)
    cd = r.headers.get('content-disposition')
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)
    if len(fname) == 0:
        return None
    return fname[0]


def unzip(infile):
    # Unzip folder
    if zipfile.is_zipfile(infile):
        print('Unzipping {}...'.format(os.path.split(infile)[1]))
        with zipfile.ZipFile(infile) as zipf:
            zipfilelist = [info.filename for info in zipf.infolist()]
            listcheck = [f for f in zipfilelist if os.path.exists(os.path.join(infile, f))]
            if len(listcheck) > 0:
                print('Overwriting {}...'.format(', '.join(listcheck)))
            zipf.extractall(os.path.split(infile)[0])
        del zipf
    else:
        raise ValueError('Not a zip file')


def dlfile(url, outpath, outfile=None, fieldnames=None):
    """Function to download file from URL path and unzip it.
    URL (required): URL of file to download
    outpath (required): the full path including
    outfile (optional): the output name without file extension, otherwise gets it from URL
    fieldnames (optional): fieldnames in output table if downloading plain text"""

    try:
        if is_downloadable(url):  # check that url is not just html
            # Get output file name
            if outfile is None:
                outfile = get_filename_from_cd(url)
                if outfile is not None:
                    out = os.path.join(outpath, outfile)
                else:
                    out = os.path.join(outpath, os.path.split(url)[1])
            else:
                out = os.path.join(outpath, outfile + os.path.splitext(url)[1])
            del outfile

            # http request
            f = requests.get(url, allow_redirects=True)
            print "downloading " + url

            # Open local file for writing
            if not os.path.exists(out):
                if 'csv' in f.headers.get('content-type').lower():  # If csv file
                    df = pd.read_csv(io.StringIO(f.text))
                    df.to_csv(out, index=False)

                elif f.headers.get('content-type').lower() == 'text/plain':  # If plain text
                    dialect = csv.Sniffer().sniff(f.text)
                    txtF = csv.DictReader(f.text.split('\n'),
                                          delimiter=dialect.delimiter,
                                          fieldnames=fieldnames)
                    with open(out, "wb") as local_file:
                        writer = csv.DictWriter(local_file, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in txtF:
                            writer.writerow(row)

                elif f.headers.get('content-type').lower() == 'application/x-gzip':
                    outunzip = os.path.splitext(out)[0]

                    # Very inelegant. But trying to download and decompress in memory always messes up files
                    response = requests.get(url, stream=True)
                    if response.status_code == 200:
                        with open(out, 'wb') as f:
                            f.write(response.raw.read())
                    with gzip.GzipFile(out, 'rb') as input:
                        s = input.read()
                        with open(outunzip, 'wb') as output:
                            output.write(s)

                else:  # Otherwise, just try reading
                    try:  # Try writing to local file
                        with open(out, "wb") as local_file:
                            local_file.write(f.read())
                        # Unzip downloaded file
                        try:
                            unzip(out + '.zip')
                        except:
                            z = zipfile.ZipFile(io.BytesIO(f.content))
                            if isinstance(z, zipfile.ZipFile):
                                z.extractall(os.path.split(out)[0])
                    except Exception:  # If fails and is zip, directly download zip in memory
                        print('Try downloading zip in memory...')
                        os.remove(out)
                        if os.path.splitext(url)[1] == '.zip':
                            z = zipfile.ZipFile(io.BytesIO(f.content))
            else:
                print('{} already exists...'.format(out))
        else:
            print('File not downloadable...')

    # handle errors
    except requests.exceptions.HTTPError, e:
        print "HTTP Error:", e.code, url
    except Exception:
        traceback.print_exc()
        if os.path.exists(out):
            os.remove(out)


def APIdownload(baseURL, workspace, basename, itersize, IDlist, geometry):
    IDrange = range(IDlist[0], IDlist[1], itersize)
    arcpy.env.workspace = workspace

    for i, j in itertools.izip(IDrange, IDrange[1:]):
        itertry = itersize
        downlist = []
        # It seems like REST API server also has a limitation on the size of the downloads so sometimes won't allow
        # Therefore, if fails to download, try smaller increments until reaches increments of 2 if still fails at increments
        # of 2, then throw a proper error and break
        while True:
            try:
                IDrangetry = list(
                    sorted(set(range(i, j + 1, itertry) + [j])))  # To make sure that the list goes until the maximum
                # Loop with smaller increment within range
                for k, l in itertools.izip(IDrangetry, IDrangetry[1:]):
                    print('From {0} to {1}'.format(k, l))
                    where = "OBJECTID>={0} AND OBJECTID<{1}".format(k, l)
                    # &geometryType=esriGeometryPoint
                    query = "?where={0}&returnGeometry={1}&f=json&outFields=*".format(where, str(geometry).lower())
                    fsURL = baseURL + query
                    if geometry == True:
                        fs = arcpy.FeatureSet()
                    elif geometry == False:
                        fs = arcpy.RecordSet()
                    else:
                        raise ValueError('Invalid geometry argument: only boolean values are accepted')
                    fs.load(fsURL)
                    if long(arcpy.GetCount_management(fs)[0]) > 0:
                        outname = '{0}_{1}_{2}'.format(basename, k, l)
                        downlist.append(outname)
                        if geometry == True:
                            arcpy.CopyFeatures_management(fs, outname)
                        else:
                            arcpy.CopyRows_management(fs, '{}.csv'.format(outname))
                        print(outname)
                    else:
                        print('No data from OBJECTID {0} to OBJECTID {1}'.format(k, l))
                break
            except:
                if itertry > 5:
                    print(
                    'Count not download, delete previous {0} datasets, try downloading in smaller increments'.format(
                        len(downlist)))
                    if len(downlist) > 0:
                        for fc in downlist:
                            arcpy.Delete_management(fc)
                        downlist = []
                    itertry = itertry / 2
                else:
                    e = sys.exc_info()[1]
                    print('Exit with error: ' + e.args[0])
                    # sys.exit(1)
                    break


def downloadroads(countyfipslist, year=None, outdir=None):
    if year is None:
        year = 2018
    if year < 2008:
        raise ValueError("Roads data are not currently available via FTP for years prior to 2008.")
    if outdir == None:
        print('Downloading to {}...'.format(os.getcwd()))
        outdir = os.getcwd()
    if not os.path.exists(outdir):
        print('Creating {}...'.format(outdir))
        os.mkdir(outdir)

    # Open ftp connection
    try:
        failedlist = []
        rooturl = "ftp://ftp2.census.gov/geo/tiger/TIGER{0}/ROADS".format(year)
        urlp = urlparse.urlparse(rooturl)
        ftp = ftplib.FTP(urlp.netloc)
        ftp.login()
        ftp.cwd(urlp.path)

        # Iterate over county fips
        x = 0
        N = len(countyfipslist)

        for county_code in countyfipslist:
            outfile = "tl_{0}_{1}_roads.zip".format(year, county_code)
            out = os.path.join(outdir, outfile)
            if not os.path.exists(out):
                # ftp download
                print "downloading " + outfile
                try:
                    with open(out, 'wb') as fobj:  # using 'w' as mode argument will create invalid zip files
                        ftp.retrbinary('RETR {}'.format(outfile), fobj.write)
                except Exception:
                    traceback.print_exc()
                    failedlist.append(county_code)
                    pass

                ######DID NOT WORK
                # urllib.urlretrieve(url, out)
                ######KEPT HANGING
                # try:
                #     r = urllib2.urlopen(url)
                #     print "downloading " + url
                #     with open(out, 'wb') as f:
                #         shutil.copyfileobj(r, f)
                # finally:
                #     if r:
                #         r.close()

                ######KEPT HANGING
                # with contextlib.closing(urllib2.urlopen(url)) as ftprequest:
                #     print "downloading " + url
                #     with open(out, 'wb') as local_file:
                #         shutil.copyfileobj(ftprequest, local_file)

            else:
                print('{} already exists... skipping'.format(outfile))
            x += 1
            print('{}% of county-level data downloaded'.format(100 * x / N))
    finally:
        if ftp:
            ftp.quit()
        if len(failedlist) > 0:
            print('{} failed to download...'.format(','.join(failedlist)))


def downloadNARR(folder, variable, years, outdir=None):
    if outdir == None:
        print('Downloading to {}...'.format(os.getcwd()))
        outdir = os.getcwd()
    if not os.path.exists(outdir):
        print('Creating {}...'.format(outdir))
        os.mkdir(outdir)

    # Open ftp connection
    try:
        failedlist = []
        rooturl = "ftp://ftp.cdc.noaa.gov/Datasets/NARR/{0}".format(folder)
        urlp = urlparse.urlparse(rooturl)
        ftp = ftplib.FTP(urlp.netloc)
        ftp.login()
        ftp.cwd(urlp.path)

        # Iterate over county fips
        x = 0
        N = len(years)

        for year in years:
            outfile = "{0}.{1}.nc".format(variable, year)
            out = os.path.join(outdir, outfile)
            if not os.path.exists(out):
                # ftp download
                print "downloading " + outfile
                try:
                    with open(out, 'wb') as fobj:  # using 'w' as mode argument will create invalid zip files
                        ftp.retrbinary('RETR {}'.format(outfile), fobj.write)
                except Exception:
                    traceback.print_exc()
                    failedlist.append(outfile)
                    # Remove empty dataset
                    if os.path.exists(out) and os.path.getsize(out) == 0L:
                        os.remove(out)
                    pass
            else:
                print('{} already exists... skipping'.format(outfile))
            x += 1
            print('{}% of data downloaded'.format(100 * x / N))
    finally:
        if ftp:
            ftp.quit()
        if len(failedlist) > 0:
            print('{} failed to download...'.format(','.join(failedlist)))