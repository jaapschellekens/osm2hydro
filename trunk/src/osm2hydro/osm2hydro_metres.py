#!/usr/bin/env python
"""
Execute the osm2hydro script.  Experimental version with different units. At some stage
this can run at the global level

Usage:
    osm2hydro_matres.py  -c inifilename [-E extent][-W working dir][-C caseFolder]
                  [-N caseName][-O osm_extract_file][-F osmfile]

    -c inifilename
    Optional (overrides the ini settings):
        -E  extent according to: [xmin, ymin, xmax, ymax]
        -W working_directory
        -C caseFolder
        -N caseName
        -O osm_extract_file
        -F osmfile
        
Description:
    osm2hydro generates hydrological and hydraulic model schematisations 
    from OpenStreetMap data

Status:
    Still under construction by Jaap Schellekens, Hessel Winsemius, Jan Talsma, 
    Ferdinand Diermanse, Ruben Dahm, Reinder Brolsma and Daniel Tollenaar
    TODO: cleanup script, move stuff to functions

Dependencies:
    For osm2hydro, it is mandatory that GDAL version 1.10 or higher is installed.
    
    

$Id: osm2hydro_metres.py 9857 2013-12-09 13:48:11Z schelle $
$Date: 2013-12-09 14:48:11 +0100 (Mon, 09 Dec 2013) $
$Author: schelle $
$Revision: 9857 $
$HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/python/applications/osm2hydro/src/osm2hydro/osm2hydro_metres.py $
$Keywords: $

"""

"""
 Copyright notice

    Copyright 2013, 2014 Hesssel Winsemius, Jaap Schellekens, Gennadii Donchyts

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    

 $Id: OSM2hydro.py 798 2013-10-02 07:11:23Z schelle $
 $Date: 2013-10-02 09:11:23 +0200 (Wed, 02 Oct 2013) $
 $Author: schelle $
 $Revision: 798 $
 $HeadURL: https://repos.deltares.nl/repos/Hydrology/trunk/OpenStreetMaps/src/OSM2hydro/OSM2hydro.py $
 $Keywords: $

"""

import osgeo.osr as osr
import osgeo.gdal as gdal
import numpy as np
import os
import glob
import logging
import logging.handlers
import pyproj
import sys
import platform
import ConfigParser
#import pylab
import getopt
import urllib
import zipfile
import gdal_merge
import poly_density
import osm2shp
import gdal_density
import dem_filter
import pdb
import scipy.signal as signal
import map2shape

debuglog = True
srs = None

def usage(*args):
    """
    Print usage information

    @param *args: command line arguments given

    """
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)



def configget(log,config,section,var,default):
    """   
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value 
    with the default in the config-file
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to get
        - default - default string
        
    Returns:
        - string - either the value from the config file or the default value
    """
    
    Def = False
    try:
        ret = config.get(section,var)
    except:
        Def = True
        ret = default
        log.info( "returning default (" + default + ") for " + section + ":" + var)
        configset(config,section,var,default, overwrite=False)
    
    default = Def
    return ret       

def configset(config,section,var,value, overwrite=False):
    """   
    Sets a string in the in memory representation of the config object
    Deos NOT overwrite existing values if overwrite is set to False (default)
    
    Input:
        - config - python ConfigParser object
        - section - section in the file
        - var - variable (key) to set
        - value - the value to set
        - overwrite (optional, default is False)
   
    Returns:
        - nothing
        
    """
    
    if not config.has_section(section):
        config.add_section(section)
        config.set(section,var,value)
    else:     
        if not config.has_option(section,var):
            config.set(section,var,value)
        else:
            if overwrite:
                config.set(section,var,value)

def iniFileSetUp(configfile):
    """
    Reads .ini file and sets default values if not present
    """
    # TODO: clean up wflwo specific stuff
    #setTheEnv(runId='runId,caseName='caseName)
    # Try and read config file and set default options
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config

def setlogger(logfilename, logReference):
    """
    Set-up the logging system. Exit if this fails
    input:
        logfilename:    string, referring to the logfile
        logReference:   string, referring to reference used in log lines
    output:
        ch:             handle, refer to logging object
        logger:         logger object
    """
    try:
        #create logger
        logger = logging.getLogger(logReference)
        logger.setLevel(logging.DEBUG)
        ch = logging.handlers.RotatingFileHandler(logfilename,maxBytes=10*1024*1024, backupCount=5)
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger, ch
    except IOError:
        print "ERROR: Failed to initialize logger with logfile: " + logfilename
        sys.exit(2)

def closeLogger(logger, ch):
    """
    Closes the logger
    """
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch

def recursive_glob(rootdir='.', suffix=''):
    """
    Prepares a list of files in location rootdir, with a suffix
    input:
        rootdir:    string, path-name
        suffix:     suffix of required files
    output:
        fileList:   list-strings, file paths and names
    """

    fileList = [os.path.join(rootdir, filename)
            for rootdir, dirnames, filenames in os.walk(rootdir)
            for filename in filenames if filename.endswith(suffix)]
    return fileList

def readMap(fileName, fileFormat):
    """ Read geographical file into memory
    """
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(fileFormat)
    mapFormat.Register()
    ds = gdal.Open(fileName)
    if ds is None:
        print 'Could not open ' + fileName + '. Something went wrong!! Shutting down'
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX    = geotrans[1]
    resY    = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2,originX+resX/2+resX*(cols-1),cols)
    y = np.linspace(originY+resY/2,originY+resY/2+resY*(rows-1),rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1) # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0,0,cols,rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal

def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

        # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif',data.shape[1],data.shape[0],1,gdal.GDT_Float32)
    # Give georeferences
    xul = x[0]-(x[1]-x[0])/2
    yul = y[0]+(y[0]-y[1])/2
    TempDataset.SetGeoTransform( [ xul, x[1]-x[0], 0, yul, 0, y[1]-y[0] ] )
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data,0,0)
    TempBand.FlushCache()
    # This seesm to happen somtimes
    if FillVal == None:
        FillVal = 1E31
        
    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    outDataset.SetProjection(srs)

    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'

def cutMap(xmin, ymin, xmax, ymax, origMap,transMap):
    """
    translateMap translates the tif formatted files from the srtm to PCRaster .map files
    origMap = your_map.tif
    transMap= your_map.map
    """
    command= 'gdal_translate -a_nodata -32768 -ot Float32 -projwin %f %f %f %f %s %s' % (xmin, ymax, xmax, ymin, origMap, transMap)
    os.system(command)

def projectMap(src_proj, trg_proj, origMap,transMap):
    """
    projectMap projects a source map to a targetmap with defined srs
    src_proj = projection of source data
    trg_proj = projection of target data
    origMap = your_map.tif
    transMap= your_map.map
    """
    command= 'gdalwarp -s_srs "%s" -t_srs "%s" -srcnodata -32768 -dstnodata -9999 %s %s' % (src_proj, trg_proj, origMap, transMap)
    print command
    os.system(command)

def filter_shape(attribute, value, origShape, transShape):
    """
    Filters a shapefile based on a selection of and attribute/value
    attribute:      string - indicating the attribute
    value:          string - the value or comma-separated list of values of the attribute, filtered out
    origShape:      strong - location of original shapefile
    transShape:     string - location of target shapefile
    """
    value_list = value.split(',')
    value_string = ''
    for num, value_entry in enumerate(value_list):
        value_string = str(value_string + '\'' + value_entry + '\'')
        if num!= len(value_list)-1:
            value_string = str(value_string + ',')
        
    # value_string = str(value_string + '\'')
    command = 'ogr2ogr -where "%s in (%s)" %s %s' % (attribute, value_string, transShape, origShape)
    os.system(command)

def burn_lines(shapeFile, mapFile, value, x, y):
    """
    Make a map with zeros and burned in values at locations where a line shape is present
    shapeFile:      string - Shapefile with line elements that need to be burned
    mapFile:        string - Target GeoTIFF map file that will be written in this step
    value:          string - The value that should be burned into the map
    x:              1-D array - x-axis of target map
    y:              1-D array - y-axis of target map
    """
    # prepare the destination map
    maskMap = np.zeros((len(y),len(x)))
    writeMap(mapFile, 'GTiff', x, y, maskMap, -9999.)
    # burn in values with gdal_rasterize
    layer = os.path.split(shapeFile)[1].split('.')[0]
    command = 'gdal_rasterize -burn %g -l %s %s %s' % (value, layer, shapeFile, mapFile)
    print command    
    # run the command
    os.system(command)

def removeFiles(wildCard):
    filelist = glob.glob(wildCard)
    for filename in filelist:
        os.unlink(filename)

def retrieve_SRTM(url_prefix, file_prefix, url_suffix, demLoc, case, xmin, ymin, xmax, ymax, logger):
    # url_prefix=http://droppr.org/srtm/v4.1/6_5x5_TIFs/
    # url_prefix=ftp://srtm.csi.cgiar.org/SRTM_v41/SRTM_Data_GeoTIFF/
    tileMinX=np.int((np.round(xmin) + 180 ) / 5 + 1)
    tileMaxX=np.int((np.round(xmax) + 180 ) / 5 + 1)
    tileMinY=np.int((60 - np.round(ymax) ) / 5 + 1)
    tileMaxY=np.int((60 - np.round(ymin) ) / 5 + 1)
    tileAvX=(xmin + xmax) / 2
    tileAvY=(ymin + ymax) / 2
    logger.info(str('Retrieving DEM tiles minX: %3.2f maxX: %3.2f, minY: %3.2f, maxY: %3.2f') % (tileMinX, tileMaxX, tileMinY, tileMaxY))
    # compute UTM zone
    utm_zone=np.int(np.round((np.round(tileAvX) + 180 ) / 6 + 1))
    if tileAvY > 0:
        proj4string='+proj=utm +zone=%g +ellps=WGS84 +units=m +no_defs ' % utm_zone
    else:
        proj4string='+proj=utm +zone=%g +south +ellps=WGS84 +units=m +no_defs ' % utm_zone
    logger.info(str('Projecting to UTM zone: ' + str(utm_zone) + ', proj4: ' + proj4string))
    tileLat=tileMinY-1
    tileLon=tileMinX-1
    for tileLon in range(tileMinX, tileMaxX+1):
        for tileLat in range(tileMinY, tileMaxY+1):
            fileName = str(file_prefix + '%02.f_%02.f' + url_suffix) % (tileLon, tileLat)
            url = url_prefix + fileName
            fileTarget = os.path.join(demLoc, fileName)
            logger.info('Retrieving ' + url)
            try:
                if os.path.exists(fileTarget):
                    logger.info("Skipping existing file: " + fileTarget)
                else:
                    urllib.urlretrieve(url,fileTarget)
                logger.info('Unzipping %s' %(fileTarget))
                zf = zipfile.ZipFile(fileTarget , 'r') 
                nameList = zf.namelist()
                for n in nameList:
                    outFile = open(os.path.join(demLoc, n), 'wb')
                    outFile.write(zf.read(n))
                    outFile.close()
                zf.close()
                os.unlink(fileTarget)
            except:
                print 'File ' + url + ' cannot be found or is corrupt, skipping!'
            # call gdal_merge to stitch everythin together
    temporary_dem = os.path.join(demLoc, 'temp.tif')
    cut_dem       = os.path.join(demLoc, case + '_latlon.tif')
    proj_dem      = os.path.join(demLoc, case + '_proj.tif')
    source_dems   = os.path.join(demLoc, 'srtm*.tif')
    gdal_merge.main(argv=['dummy','-o', temporary_dem, source_dems])
    cutMap(xmin, ymin, xmax, ymax, temporary_dem, cut_dem) # this is the final lat lon map
    projectMap('EPSG:4326',proj4string, cut_dem, proj_dem) # this is the final projected map
    # remove all redundant files
    removeFiles(os.path.join(demLoc, 'srtm*.tif'))
    removeFiles(os.path.join(demLoc, 'srtm*.hdr'))
    removeFiles(os.path.join(demLoc, 'srtm*.tfw'))
    os.unlink(os.path.join(demLoc, 'readme.txt'))
    os.unlink(os.path.join(demLoc, 'temp.tif'))
    return proj4string

def reduce_resolution(x, y, data, target_resolution):
    """
    Reduces resolution of a map (contained by x-axis, y-axis and an array) to
    a user-defined target resolution. This is done by a convolution and then picking
    one value per the amount of required resolution values.
    """
    resolution_original = x[1]-x[0]
    resolution_factor   = np.int(np.maximum(np.round(target_resolution/resolution_original), 1))
    # prepare convolution array
    conv_arr = np.ones((resolution_factor, resolution_factor))
    conv_arr = conv_arr/conv_arr.sum()
    data_averaged = signal.convolve2d(data, conv_arr, mode='same')
    y_idx = range(np.int(np.round(resolution_factor/2)), len(y),resolution_factor)
    x_idx = range(np.int(np.round(resolution_factor/2)), len(x),resolution_factor)
    xi, yi = np.meshgrid(x_idx,y_idx)
    data_reduced  = data_averaged[yi, xi]
    y_reduced     = y[y_idx]
    x_reduced     = x[x_idx]
    return x_reduced, y_reduced, data_reduced


def removeFinishedProcesses(processes):
    """ given a list of (commandString, process), 
        remove those that have completed and return the result 
    """
    newProcs = []
    for pollCmd, pollProc in processes:
        retCode = pollProc.poll()
        if retCode==None:
            # still running
            newProcs.append((pollCmd, pollProc))
        elif retCode!=0:
            # failed
            raise Exception("Command %s failed" % pollCmd)
        else:
            print "Command %s completed successfully" % pollCmd
    return newProcs

def runCommands(commands, maxCpu):
    """
    Runs a list of processes deviding
    over maxCpu
    """
    processes = []
    for command in commands:
        print "Starting process %s" % command
        proc =  subprocess.Popen(shlex.split(command))
        procTuple = (command, proc)
        processes.append(procTuple)
        while len(processes) >= maxCpu:
            time.sleep(.2)
            processes = removeFinishedProcesses(processes)

    # wait for all processes
    while len(processes)>0:
        time.sleep(0.5)
        processes = removeFinishedProcesses(processes)
    print "All processes (" + str(len(commands)) + ") completed."


"""
Main program starts below
"""
def main(opts):
    configFile = None
    caseFolder_cm = None
    caseName = None
    osmExtractName = None
    osmFileName = None
    workFolder   = os.path.dirname(__file__) # location of osm2hydro.py

    extent = None
    
    for o, a in opts:
        print o, a
        if o == '-E': exec("extent = " + a) # ,globals(), globals()
        if o == '-W': workFolder = a 
        if o == '-c': configFile =  a
        if o == '-C': caseFolder_cm =  a
        if o == '-N': caseName =  a
        if o == '-O': osmExtractName =  a
        if o == '-F': osmFileName =  a
    
    if configFile == None:
         usage()
         sys.exit(2) 
        
    
    #    noArgExitStr = str('****************      osm2hydro v. 1.0       *****************\n' + \
    #    '$Id: osm2hydro_metres.py 9857 2013-12-09 13:48:11Z schelle $\n' + \
    #    '$Date: 2013-12-09 14:48:11 +0100 (Mon, 09 Dec 2013) $\n' + \
    #    '$Author: schelle $\n' + \
    #    '$Revision: 9857 $\n' + \
    #    'Welcome to osm2hydro. With this software, you can generate hydrological\n ' + \
    #    'models and hydraulic models on-the-fly based on SRTM data and OpenStreetMaps.\n ' + \
    #    'As model concepts, WFLOW and SOBEK are used\n ' + \
    #    'osm2hydro runs as follows: osm2hydro <config-file.ini> -<option1> -<option2> ....\n' + \
    #    '<config-file.ini>:       a osm2hydro configuration file, generated by you,\n' + \
    #    '                         the end-user.\n' + \
    #    '<option1>:               a run-time argument, used to select specific\n' + \
    #    '                         run-time options.\n' + \
    #    '***********************************************************\n' + \
    #    'For osm2hydro, it is mandatory that GDAL version 1.10 or higher is installed.\n' + \
    #    '- For Linux, please build the required GDAL source or install through your package manager.\n' + \
    #    '- For windows, the correct installation file can be downloaded through this link: \n' + \
    #    'http://www.gisinternals.com/sdk/Download.aspx?file=release-1600-gdal-1-9-2-mapserver-6-2-0\gdal-19-1600-core.msi\n' + \
    #    '****************         Have fun!        *****************')
    #    print noArgExitStr
    #    sys.exit(2)
    
    #Defaults
    run_dem       = False
    run_osm       = True
    SkipShapes    = False
    SkipOsmExtract=False
    run_hydrology = False
    run_hydraulics = True
    method_poly = False
    mergeroads=True
    above_is_paved=90
    width_road_main=0.5
    width_road_sec=0.1
    width_road_small=0.1
    width_road_track=0.01
    width_waterway_river=0.05
    width_waterway_riverbank=0.025
    width_waterway_stream=0.008
    width_waterway_canal=0.008
    width_waterway_drain=0.004
    width_waterway_ditch=0.001
    # Fixed inputs DEM
    url_prefix='ftp://xftp.jrc.it/pub/srtmV4/tiff/'
    file_prefix='srtm_'
    url_suffix='.zip'
    # fixed inputs for Hydraulics
    filter_name='bior6.8'
    proj4string_src = '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'
    
    #### DEBUG LINES ####
    
    """ Read configuration file below"""
    
    try:
        config = iniFileSetUp(configFile)
        print 'Reading options from ' + configFile    
        # first read the location of the output folder. Start a log file in the output folder.
        try:
            caseFolder = os.path.abspath(config.get('case','caseFolder'))
        except:
            print 'Entry "caseFolder" not given by config file. Exiting!'
            sys.exit(2)
        # Override with command-line
        if caseFolder_cm:
            caseFolder = caseFolder_cm
        try:
            case = config.get('case','caseName')
        except:
            print 'Entry "caseName" not given by config file. Exiting!'
            sys.exit(2)
            # prepare output folders
        # Override with command-line
        if caseName:
            case = caseName
        
    
        outputFolder = os.path.join(caseFolder, case)
        if not os.path.isdir(outputFolder):
            try:
                os.makedirs(outputFolder)
            except:
                print 'outputFolder "' + outputFolder + '" is invalid. Please check availability of path. Exiting!'
                sys.exit(2)
        # set the osm2hydro logger
        logger, ch      = setlogger(os.path.join(outputFolder,'OMS2Hydro.log'),'osm2hydro')
        logger.info('Starting osm2hydro v. 1.0')
        logger.info('$Id: osm2hydro_metres.py 9857 2013-12-09 13:48:11Z schelle $')
        logger.info('$Date: 2013-12-09 14:48:11 +0100 (Mon, 09 Dec 2013) $')
        logger.info('$Author: schelle $')
        logger.info('$Revision: 9857 $')
        logger.info('Command-line options: ' + str(opts))
    
        maxCPU = int(configget(logger,config,'osm','maxCPU','2'))
        try:
            osmFile = os.path.abspath(configget(logger,config,'osm','osmFile','default.osm'))
            osmExtract = osmFile + '_extract.pbf'
        except:
            logger.error('Osmfile not found (' + osmFile + '). Exiting!');closeLogger(logger, ch);sys.exit(2)
            
        # override with command-line
        if osmExtractName:
            osmExtract= osmExtractName
        if osmFileName:
            osmFile = osmFileName
        try:
            osmConfig = os.path.abspath(configget(logger,config,'osm','osmConfig','osm2shp.ini'))
        except:
            logger.error('osmConfig (' + osmConfig +' ) not found Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            xmin = np.float(config.get('geography','xmin'))
        except:
            logger.error('Entry "xmin" not given by config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            xmax = np.float(config.get('geography','xmax'))
        except:
            logger.error('Entry "xmax" not given by config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            ymin = np.float(config.get('geography','ymin'))
        except:
            logger.error('Entry "ymin" not given by config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            ymax = np.float(config.get('geography','ymax'))
        except:
            logger.error('Entry "ymax" not given by config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            resolution = np.float(config.get('geography','resolution'))
        except:
            logger.error('Entry "ymax" not given by config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        # sections
        try:
            exec("run_hydrology = " + configget(logger,config,'hydrology','run',str(run_hydrology)))
        except:
            logger.error('Problem with Entry "run" in section "hydrology" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            exec("run_osm = " + configget(logger,config,'osm','run',str(run_osm)))
        except:
            logger.error('Problem with Entry "run" in section "osm" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            exec("run_dem = " + configget(logger,config,'dem','run',str(run_dem)))
        except:
            logger.error('problem with Entry "run" in section "dem" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)        
        try:
            exec("run_hydraulics = " + configget(logger,config,'hydraulics','run',str(run_hydraulics)))
        except:
            logger.error('Problem with Entry "run" in section "hydraulics" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            streams = configget(logger,config,'hydraulics','streams','river')
        except:
            logger.error('Problem with Entry "streams"  in section "hydraulics" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            filter_weight = np.float(configget(logger,config,'hydraulics','filter_weight','1.0'))
        except:
            logger.error('Problem with Entry "filter_weight"  in section "hydraulics" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        try:
            burn_value = np.abs(np.float(configget(logger,config,'hydraulics','burn_value','-999.0')))
        except:
            logger.error('Problem with Entry "burn_value" in section "hydraulics" in config file. Exiting!');closeLogger(logger, ch);sys.exit(2)
        # Now read all shapefiles into a list with tuples
        elev_feat         = []
        try:
            elev_featStr      = config.items('elevated_features')             # Names of the variables, as they occur in the NetCDF files of the GCM model outputs
            for feat in elev_featStr:
                elev = np.float(feat[1])
                # write output unit, input unit, slope and offset to list of unitMaps
                elev_feat.append((feat[0] + '.shp',elev))
                # add units to full list of units
        except:
            logger.warning('No elevated features found. So elevation model will not be manipulated to include elevated feats such as roads and dikes')
        
    except:
        logger.error( 'Config file "' + configFile + '" does not exist or is improperly formatted.')
        sys.exit(2)
        
    # Override ini file settings with command-line options
    if extent != None:
        #[xmin, ymin, xmax, ymax]
        xmin = extent[0]
        ymin = extent[1]
        xmax = extent[2]
        ymax = extent[3]
    

    if configget(logger,config,'osm','poly',str(method_poly)) == 'True':
        method_poly=True

    if configget(logger,config,'osm','SkipShapes',str(SkipShapes)) == 'True':
        SkipShapes=True

    if configget(logger,config,'osm','SkipOsmExtract',str(SkipOsmExtract)) == 'True':
        SkipOsmExtract=True

    print 'skip OSM' + str(SkipOsmExtract)
    mergeroads=bool(configget(logger,config,'osm','mergeroads',str(mergeroads)))
    
    resamp_grid_method=float(configget(logger,config,'osm','resamp',"10.0"))
    width_road_main=float(configget(logger,config,'osm','width_road_main',str(width_road_main)))
    width_road_sec=float(configget(logger,config,'osm','width_road_sec',str(width_road_sec)))
    width_road_small=float(configget(logger,config,'osm','width_road_small',str(width_road_small)))
    width_road_track=float(configget(logger,config,'osm','width_road_track',str(width_road_track)))
    
    width_waterway_river=float(configget(logger,config,'osm','width_waterway_river',str(width_waterway_river)))
    width_waterway_riverbank=float(configget(logger,config,'osm','width_waterway_riverbank',str(width_waterway_riverbank)))
    width_waterway_stream=float(configget(logger,config,'osm','width_waterway_stream',str(width_waterway_stream)))
    width_waterway_canal=float(configget(logger,config,'osm','width_waterway_canal',str(width_waterway_canal)))
    width_waterway_drain=float(configget(logger,config,'osm','width_waterway_drain',str(width_waterway_drain)))
    width_waterway_ditch=float(configget(logger,config,'osm','width_waterway_ditch',str(width_waterway_ditch)))
    above_is_paved=float(configget(logger,config,'osm','above_is_paved',str(above_is_paved)))
    
    """ Read configuration file done!!!"""
    # No convertt to size for the hires grid
#    width_road_main=width_road_main*resamp_grid_method
#    width_road_sec=width_road_sec*resamp_grid_method
#    width_road_track=width_road_track*resamp_grid_method
#    width_road_small=width_road_small*resamp_grid_method
#    width_waterway_river=width_waterway_river*resamp_grid_method
#    width_waterway_riverbank=width_waterway_riverbank*resamp_grid_method
#    width_waterway_stream=width_waterway_stream*resamp_grid_method
#    width_waterway_canal=width_waterway_canal*resamp_grid_method
#    width_waterway_drain=width_waterway_drain*resamp_grid_method
#    width_waterway_ditch=width_waterway_ditch*resamp_grid_method
#    
    """ Start osm2hydro"""
    
    # =========================================================================
    # establish root folder
    rootLoc        = os.path.split(sys.argv[0])[0]
    
    # establish run-time executables
    
    print platform.system()

    os.environ['OSM_CONFIG_FILE']=os.path.join(os.path.abspath(workFolder), 'osmconf.ini')
    
    if platform.system() == "Linux" or platform.system() == "Linux2":
        osm2shp_exe    = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','osm2shp'))
        osmconvert_exe = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','osmconvert'))
        osmfilter_exe  = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','osmfilter'))
        shptree_exe  = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','shptree'))
    else:
        osm2shp_exe    = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','osm2shp.exe'))
        osmconvert_exe = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','osmconvert.exe'))
        osmfilter_exe  = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','osmfilter.exe'))
        shptree_exe  = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','shptree.exe'))
        
    
    # prepare output folders
    demLoc         = os.path.abspath(os.path.join(outputFolder,'dem'))
    osmshapesLoc   = os.path.abspath(os.path.join(outputFolder,'osmshapes'))
    modelshapesLoc = os.path.abspath(os.path.join(outputFolder,'modelshapes'))
    gridsLoc       = os.path.abspath(os.path.join(outputFolder,'grids'))
    modelgridsLoc  = os.path.abspath(os.path.join(outputFolder,'modelgrids'))
    try:
        os.makedirs(demLoc)
    except:
        pass # happens if folder already exists
    try:
        os.makedirs(osmshapesLoc)
    except:
        pass # happens if folder already exists
    try:
        os.makedirs(modelshapesLoc)
    except:
        pass # happens if folder already exists
    try:
        os.makedirs(gridsLoc)
    except:
        pass # happens if folder already exists
    
    try:
        os.makedirs(modelgridsLoc)
    except:
        pass # happens if folder already exists
    
    # now run the SRTM downloader
    if run_dem:
        logger.info("Retrieving SRTM dem...")
        retrieve_SRTM(url_prefix, file_prefix, url_suffix, demLoc, case, xmin, ymin, xmax, ymax, logger)
    
    
    
    if run_osm:
        if not SkipOsmExtract:
          logger.info("Processing OSM file...")
          # first derive extracted part of osm file, -0.1 and + 0.1 are used to select within a larger area to minimize boundary effects --complex-ways --complete-ways do not do the job
          command = str('%s %s  -v  --out-pbf --hash-memory=1000 --drop-version --drop-relations --complex-ways --complete-ways --drop-author -b="%f,%f,%f,%f" -o=%s') % (osmconvert_exe, osmFile, xmin-0.1, ymin-0.1, xmax+0.1, ymax+0.1, osmExtract)
          print command

          #command = str('%s %s --drop-broken-refs  -v  -b="%f,%f,%f,%f" -o=%s') % (osmconvert_exe, osmFile, xmin, ymin, xmax, ymax, osmExtract)
          os.system(command)
          osmFile = osmExtract
        else:
          logger.info("Skipping osm extraction...")
          
          
        if not SkipShapes:
          logger.info("Creating shapes...")
          #command = osm2shp_exe + ' -c ' + osmConfig + ' -d ' + osmshapesLoc + ' ' + osmFile
          if os.path.isfile(osmFile):
              osm2shp.main(['-o',osmshapesLoc,'-O',osmFile,'-c',osmConfig,'-M',maxCPU])
          else:
              logger.error('OSM extract file "' + osmFile + '" is not available, closing...');closeLogger(logger, ch);sys.exit(2)
    
          #['-S',lu_water_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_water_tif, '-F', 'GTiff','-r',str(resamp_grid_method)]
          #os.system(command)
        else:
          logger.info("Skipping shapefile creation...")
          
        lu_water_shp   = os.path.join(osmshapesLoc, 'lu_water.shp')
        lu_paved_shp   = os.path.join(osmshapesLoc, 'lu_paved.shp')
        lu_unpaved_shp = os.path.join(osmshapesLoc, 'lu_unpaved.shp')
        lu_roads_main_shp  = os.path.join(osmshapesLoc, 'lu_roads_main.shp')
        lu_roads_track_shp  = os.path.join(osmshapesLoc, 'lu_roads_track.shp')
        lu_roads_small_shp  = os.path.join(osmshapesLoc, 'lu_roads_small.shp')
        lu_roads_sec_shp  = os.path.join(osmshapesLoc, 'lu_roads_sec.shp')
        lu_waterway_river_shp  = os.path.join(osmshapesLoc, 'lu_waterway_river.shp')
        lu_waterway_riverbank_shp  = os.path.join(osmshapesLoc, 'lu_waterway_riverbank.shp')
        lu_waterway_stream_shp  = os.path.join(osmshapesLoc, 'lu_waterway_stream.shp')
        lu_waterway_canal_shp  = os.path.join(osmshapesLoc, 'lu_waterway_canal.shp')
        lu_waterway_drain_shp  = os.path.join(osmshapesLoc, 'lu_waterway_drain.shp')
        lu_waterway_ditch_shp  = os.path.join(osmshapesLoc, 'lu_waterway_ditch.shp')
        lu_buildings_shp  = os.path.join(osmshapesLoc, 'lu_buildings.shp')
        
        # now create a list of all waterways, so that they can be merged together in one file
        lu_waterway_list = [lu_waterway_river_shp, lu_waterway_stream_shp, lu_waterway_canal_shp, lu_waterway_drain_shp, lu_waterway_ditch_shp]
        lu_waterway_all  = os.path.join(osmshapesLoc, 'lu_waterway_all.shp')
         # now create indexes to speed up processing
        if not SkipShapes:
            # now initiate a new file, containing the first layer of the list of waterway files
            command = 'ogr2ogr ' + lu_waterway_all + ' ' + lu_waterway_list[0]
            os.system(command)
            # now add all remaining files from the list
            for lu_waterway_layer in lu_waterway_list[1:]:
                command = 'ogr2ogr -update -append ' + lu_waterway_all + ' ' + lu_waterway_layer + ' -nln ' + os.path.split(lu_waterway_all)[1].split('.')[0]
                os.system(command)
       
      
        lu_water_tif   = os.path.join(gridsLoc, 'lu_water.tif')
        lu_paved_tif   = os.path.join(gridsLoc, 'lu_paved.tif')
        lu_unpaved_tif = os.path.join(gridsLoc, 'lu_unpaved.tif')
        roads_den_tif  = os.path.join(gridsLoc, 'roads_den.tif')
        lu_pavedpol_tif  = os.path.join(gridsLoc, 'lu_pavedpol.tif')
        lu_roads_main_tif = os.path.join(gridsLoc, 'lu_roads_main.tif')
        lu_roads_sec_tif = os.path.join(gridsLoc, 'lu_roads_sec.tif')
        lu_roads_small_tif = os.path.join(gridsLoc, 'lu_roads_small.tif')
        lu_roads_track_tif = os.path.join(gridsLoc, 'lu_roads_track.tif')
        lu_roads_tif = os.path.join(gridsLoc, 'lu_roads.tif')
        lu_waterway_river_tif = os.path.join(gridsLoc, 'lu_waterway_river.tif')
        lu_waterway_riverbank_tif = os.path.join(gridsLoc, 'lu_waterway_riverbank.tif')
        lu_waterway_stream_tif = os.path.join(gridsLoc, 'lu_waterway_stream.tif')
        lu_waterway_canal_tif = os.path.join(gridsLoc, 'lu_waterway_canal.tif')
        lu_waterway_ditch_tif = os.path.join(gridsLoc, 'lu_waterway_ditch.tif')
        lu_waterway_drain_tif = os.path.join(gridsLoc, 'lu_waterway_drain.tif')
        lu_paved_fromroads_tif = os.path.join(gridsLoc, 'lu_paved_fromroads.tif')
        lu_buildings_tif = os.path.join(gridsLoc, 'lu_buildings.tif')
    
        # convert shapes to fractional tiffs per land cover type
        logger.info("Converting shape files to grids...")
        
        if method_poly:
            poly_density.main(['-S',lu_water_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax), '-C',str(resolution),'-o', lu_water_tif, '-F', 'GTiff','-1'])
            poly_density.main(['-S',lu_paved_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax), '-C',str(resolution),'-o', lu_paved_tif, '-F', 'GTiff','-1'])
            poly_density.main(['-S',lu_unpaved_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax), '-C',str(resolution),'-o', lu_unpaved_tif, '-F', 'GTiff','-1'])
            gdal_density.main(['-S',lu_waterway_river_shp,'-b',str(width_waterway_river),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_river_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_waterway_riverbank_shp,'-b',str(width_waterway_riverbank),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_riverbank_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_waterway_stream_shp,'-b',str(width_waterway_stream),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_stream_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_waterway_canal_shp,'-b',str(width_waterway_canal),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_canal_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_waterway_ditch_shp,'-b',str(width_waterway_ditch),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_ditch_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_waterway_drain_shp,'-b',str(width_waterway_drain),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_drain_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_main_shp,'-b',str(width_road_main),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_main_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_sec_shp,'-b',str(width_road_sec),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_sec_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_small_shp,'-b',str(width_road_sec),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_small_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_main_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_main_tif + "_den.tif", '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_sec_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_sec_tif + "_den.tif", '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-S',lu_roads_small_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_small_tif + "_den.tif" , '-F', 'GTiff','-r',str(resamp_grid_method)])        
        else:
            gdal_density.main(['-G','True','-S',lu_water_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_water_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_paved_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_paved_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_unpaved_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_unpaved_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_river_shp,'-b',str(width_waterway_river),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_river_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_riverbank_shp,'-b',str(width_waterway_riverbank),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_riverbank_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_stream_shp,'-b',str(width_waterway_stream),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_stream_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_canal_shp,'-b',str(width_waterway_canal),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_canal_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_ditch_shp,'-b',str(width_waterway_ditch),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_ditch_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_waterway_drain_shp,'-b',str(width_waterway_drain),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_waterway_drain_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_roads_main_shp,'-b',str(width_road_main),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_main_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_roads_sec_shp,'-b',str(width_road_sec),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_sec_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_roads_small_shp,'-b',str(width_road_small),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_small_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-M','-S',lu_roads_track_shp,'-b',str(width_road_track),'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_track_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            
            gdal_density.main(['-G','True','-S',lu_roads_main_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_main_tif + "_den.tif", '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_roads_sec_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_sec_tif + "_den.tif", '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_roads_small_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_small_tif + "_den.tif" , '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_roads_track_shp,'-b',"1",'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_roads_track_tif + "_den.tif" , '-F', 'GTiff','-r',str(resamp_grid_method)])
            gdal_density.main(['-G','True','-S',lu_buildings_shp,'-E',str('[%f,%f,%f,%f]') % (xmin, ymin, xmax, ymax),'-C',str(resolution),'-o', lu_buildings_tif, '-F', 'GTiff','-r',str(resamp_grid_method)])
            
           
            
        logger.info("Creating corrrected LU maps...")
        x, y, lu_water, FillVal   = readMap(lu_water_tif,'GTiff')
        x, y, lu_waterway_river, FillVal   = readMap(lu_waterway_river_tif,'GTiff')
        x, y, lu_waterway_riverbank, FillVal   = readMap(lu_waterway_riverbank_tif,'GTiff')
        x, y, lu_waterway_stream, FillVal   = readMap(lu_waterway_stream_tif,'GTiff')
        x, y, lu_waterway_canal, FillVal   = readMap(lu_waterway_canal_tif,'GTiff')
        x, y, lu_waterway_drain, FillVal   = readMap(lu_waterway_drain_tif,'GTiff')
        x, y, lu_waterway_ditch, FillVal   = readMap(lu_waterway_ditch_tif,'GTiff')
        x, y, lu_paved, FillVal   = readMap(lu_paved_tif,'GTiff')
        x, y, lu_unpaved, FillVal = readMap(lu_unpaved_tif,'GTiff')
        x, y, lu_roads_main, FillVal = readMap(lu_roads_main_tif,'GTiff')
        x, y, lu_roads_sec, FillVal = readMap(lu_roads_sec_tif,'GTiff')
        x, y, lu_roads_small, FillVal = readMap(lu_roads_small_tif,'GTiff')
    
        x, y, den_roads_main, FillVal = readMap(lu_roads_main_tif + "_den.tif",'GTiff')
        x, y, den_roads_sec, FillVal = readMap(lu_roads_sec_tif + "_den.tif",'GTiff')
        x, y, den_roads_small, FillVal = readMap(lu_roads_small_tif + "_den.tif",'GTiff')

        if os.path.exists(lu_roads_track_tif + "_den.tif"):
            x, y, den_roads_track, FillVal = readMap(lu_roads_track_tif + "_den.tif",'GTiff')    
            roads_den = den_roads_main + den_roads_sec + den_roads_small + den_roads_track
        else:
            roads_den = den_roads_main + den_roads_sec + den_roads_small

        if os.path.exists(lu_roads_track_tif):
            x, y, lu_roads_track, FillVal = readMap(lu_roads_track_tif,'GTiff')
            lu_roads_total = lu_roads_small + lu_roads_sec + lu_roads_main + lu_roads_track
        else:
            lu_roads_total = lu_roads_small + lu_roads_sec + lu_roads_main

        ## determine percentile  (above %) of road density
        urbover = np.percentile(roads_den,above_is_paved)
        urbovermax = np.max(roads_den)
        alturb_idx = roads_den >= urbover
        alturb     = np.zeros((lu_water.shape))
        alturb[alturb_idx] = np.minimum(1.0 * (roads_den[alturb_idx] + 0.75 *(urbovermax-urbover)) /urbovermax,1.0)
        alturb = 0.8 * alturb
    
        global srs
        ds = gdal.Open(lu_water_tif, gdal.GA_ReadOnly)
        srs = ds.GetProjection()
        ds = None
        
        writeMap(lu_pavedpol_tif, 'GTiff', x, y, lu_paved, FillVal)    
        ###### Series of correction to get all to add-up to zero ##############
        # TODE: move this to seperate function
        # add roads to paved
        if (mergeroads): # Why the 0.8!!!!!!
            lu_paved = np.maximum(alturb,lu_paved)
        lu_zz = lu_paved + lu_roads_total
         # Limit pave landuse to 1 first
        lu_pv_idx = lu_zz > 1.0    
        # Subtract that from the roads
        #lu_paved = lu_zz
        lu_paved[lu_pv_idx] = 1.0
        lu_roads_total[lu_pv_idx] = 0.0
        
        
        lu_total     = lu_water + lu_paved + lu_unpaved 
        
        # Add waterways to water and max to 1.0
        lu_water = np.minimum(1.0,(lu_water + lu_waterway_river + lu_waterway_riverbank + lu_waterway_stream + lu_waterway_canal + lu_waterway_drain + lu_waterway_ditch ))
        # calculate correction for where there is water (water is most important)
        watercor_idx = np.all([lu_total > 1.0, lu_water > 0.0], axis=0)
        watercor     = np.zeros((lu_water.shape))
        watercor[watercor_idx] = np.minimum(lu_total[watercor_idx]-1.,lu_water[watercor_idx])
        
        # first try to get it from paved
        watercor_paved = np.minimum(watercor,lu_paved)
        lu_paved = lu_paved - watercor_paved
        watercor_unpaved = np.minimum((watercor - watercor_paved), lu_unpaved)        
        lu_unpaved = lu_unpaved - watercor_unpaved
        watercor = watercor - watercor_paved - watercor_unpaved
        #writeMap("un1.tif", 'GTiff', x, y, lu_unpaved, FillVal)
        lu_total        = lu_water + lu_paved + lu_unpaved
        
        # correct again, now assuming unpaved is more important
        pavedcor_idx = np.all([lu_total > 1.0],axis=0)
        pavedcor     = np.zeros((lu_total.shape))
        pavedcor[pavedcor_idx] = np.minimum(lu_total[pavedcor_idx]-1,lu_paved[pavedcor_idx])
        
        lu_paved     = lu_paved - pavedcor
        lu_total     = lu_water + lu_paved + lu_unpaved
        # now add unassigned surface
        unassigned = 1.0 - lu_total
        lu_unpaved = lu_unpaved + unassigned
        lu_total     = lu_water + lu_paved + lu_unpaved
        #writeMap("un2.tif", 'GTiff', x, y, lu_unpaved, FillVal)
        
        # Finally burn in roads in unpaved areas and waterways
        burnidx = np.all([(lu_total + lu_roads_total) > 1.0,lu_unpaved >0], axis=0)
        rdcor = np.zeros((lu_water.shape))
        rdcor[burnidx] = np.minimum(lu_unpaved[burnidx], lu_roads_total[burnidx])
        lu_unpaved = lu_unpaved - rdcor    
    
        lu_paved = lu_paved + rdcor    
        lu_total     = lu_water + lu_paved + lu_unpaved
        #writeMap("tot1.tif", 'GTiff', x, y, lu_total, FillVal)    
        # write to files
        writeMap(lu_water_tif, 'GTiff', x, y, lu_water, FillVal)
        writeMap(lu_paved_tif, 'GTiff', x, y, lu_paved, FillVal)
        writeMap(lu_unpaved_tif, 'GTiff', x, y, lu_unpaved, FillVal)
        writeMap(lu_roads_tif, 'GTiff', x, y, lu_roads_total, FillVal)
        writeMap(roads_den_tif, 'GTiff', x, y, roads_den, FillVal)
        writeMap(lu_paved_fromroads_tif, 'GTiff', x, y, alturb, FillVal)
    
    
    if run_hydraulics:
        logger.info("Running hydraulics...")
        x_av=(xmin + xmax) / 2
        y_av=(ymin + ymax) / 2
        # compute UTM zone
        utm_zone=np.int(np.round((np.round(x_av) + 180 ) / 6 + 1))
        if y_av > 0:
            proj4string_trg = '+proj=utm +zone=%g +ellps=WGS84 +units=m +no_defs ' % utm_zone
        else:
            proj4string_trg = '+proj=utm +zone=%g +south +ellps=WGS84 +units=m +no_defs ' % utm_zone
        logger.info(str('Projecting to UTM zone: ' + str(utm_zone) + ', proj4: ' + proj4string_trg))
    
        # define relevant shapefile names
        lu_waterways_all  = os.path.join(osmshapesLoc, 'lu_waterway_all.shp')
        streams_shp    = os.path.join(modelshapesLoc, 'streams.shp')
        # define relevant grid file names
        demFile        = os.path.join(demLoc, case + '_latlon.tif')
        demFilePCR     = os.path.join(demLoc, case + '_latlon.map')
        demBurnFile    = os.path.join(demLoc, case + '_burnIn.map')
        burnValueFile  = os.path.join(demLoc, case + '_burnValue.tif')
        demFilterFile  = os.path.join(demLoc, case + '_smooth.map')
        demReducedFile = os.path.join(demLoc, case + '_reduced.tif')
        dem_hydraulicsFile = os.path.join(modelgridsLoc, case + '_hydraulics.asc')
    
        # Filter out the significant streams from the waterways
        logger.info("Filtering waterways...")
        filter_shape('waterway', streams, lu_waterways_all, streams_shp)
        if not(os.path.isfile(demBurnFile)):
            # Load DEM in memory
            x, y, dem, FillVal = readMap(demFile,'GTiff')
            # write to PCRaster format
            writeMap(demFilePCR, 'PCRaster', x, y, dem, FillVal)
            dem[dem==FillVal]  = np.nan
            
            # burn rivers into the high-resolution elevation DEM
            burn_lines(streams_shp, burnValueFile, burn_value, x, y)
            # load burn map in memory and add to original DEM
            x, y, burnValueMap, FillValBurn = readMap(burnValueFile, 'GTiff')
            demBurn = dem - burnValueMap
            # save to file and remove the burn value map
            demBurn[np.isnan(demBurn)] = -9999.
            writeMap(demBurnFile, 'PCRaster', x, y, demBurn, -9999.)
            os.unlink(burnValueFile)
        # 
        # smooth the DEM!
        logger.info("Smoothing DEM...")
        if not(os.path.isfile(demFilterFile)):
            x, y, dem, FillVal = readMap(demFile,'GTiff')
            dem[dem==FillVal] = np.nan
            dem_smooth = dem_filter.dem_filter(dem, filter_name, filter_weight)
            writeMap(demFilterFile, 'PCRaster', x, y, dem_smooth, FillVal)
    
        # reduce DEM to required resolution
        if not(os.path.isfile(demReducedFile)):
            x, y, dem_smooth, FillVal = readMap(demFilterFile, 'PCRaster')
            x_reduced, y_reduced, dem_reduced = reduce_resolution(x, y, dem_smooth, resolution)
            dem[np.isnan(dem)] = FillVal
            dem_reduced[np.isnan(dem_reduced)] = FillVal
            writeMap(demReducedFile, 'GTiff', x_reduced, y_reduced, dem_reduced, FillVal)
        else:
            # in case the reduced DEM is not in memory, read it below!!
            x_reduced, y_reduced, dem_reduced, FillVal = readMap(demReducedFile, 'GTiff')
            dem_reduced[dem_reduced==FillVal] = np.nan
    
        burnValueTotal = np.zeros(dem_reduced.shape)
    
        # burn elevated features on top of smoothed and reduced elevation map
        if not(os.path.isfile(dem_hydraulicsFile)):
            for elev_feat_shp, elev_feat_value in elev_feat:
                burn_lines(os.path.join(osmshapesLoc,elev_feat_shp), burnValueFile, elev_feat_value, x_reduced, y_reduced)
                # load burn map in memory and add to original DEM
                x_reduced, y_reduced, burnValueMap, FillValBurn = readMap(burnValueFile, 'GTiff')
                # assume that the maximum value found must be used for burning
                burnValueTotal = np.maximum(burnValueMap, burnValueTotal)
                os.unlink(burnValueFile)
            dem_hydraulics = dem_reduced + burnValueTotal
            # save to file and remove the burn value map
            dem_hydraulics[np.isnan(dem_hydraulics)] = -9999.
            burnValueTotal[np.isnan(burnValueTotal)] = -9999.
            writeMap(dem_hydraulicsFile, 'AAIGrid', x_reduced, y_reduced, dem_hydraulics, -9999.)
        
        # based on the PCRaster maps generated, derive 2 shapefiles (profiles and computation points)
        rivermap      = os.path.join(demLoc, 'riversid.map')
        drainmap      = os.path.join(demLoc, 'drain.map')
        lddmap        = os.path.join(demLoc, 'ldd.map')
        depthmap      = os.path.join(demLoc, 'river_depth.map')
        widthmap      = os.path.join(demLoc, 'river_width.map')
        catchpointmap = os.path.join(demLoc, 'catch_point.map')
        catchsurfmap  = os.path.join(demLoc, 'catchments.map')
        shapedata     = os.path.join(modelshapesLoc, 'SOBEK_shapes')
        shapecatch    = os.path.join(shapedata, 'catchments.shp')
        
        # Check if ldd is present, if not: 
        # Derive a local drainage direction map from the elevation model. Separate process to prevent memory overload
        if not(os.path.isfile(lddmap)):
            pcrCommand = 'pcrcalc ' + lddmap + ' = lddcreate(' + demBurnFile + ', 1e31, 1e31, 1e31, 1e31)' 
            logger.info('Running PCRaster with command: ' + pcrCommand)
            os.system(pcrCommand)
        # check if riversid.map already exists, if not run the PCRaster script
        if not(os.path.isfile(rivermap)):
            # now start the PCRaster scripts that prepares the gridded maps of streams and so on.
            pcrModFile = os.path.join(workFolder, 'derive_catch.mod')
            pcrCommand = 'pcrcalc -f ' + pcrModFile + ' ' + demBurnFile + ' ' + demFilePCR + ' ' + lddmap + ' ' + drainmap + ' ' + demLoc
            logger.info('Running PCRaster with command: ' + pcrCommand)
            os.system(pcrCommand)
        logger.info('Deriving shape files for SOBEK model')
        map2shape.PCR_river2Shape(rivermap, drainmap, lddmap, depthmap, widthmap, demFilterFile, proj4string_src, proj4string_trg, shapedata, logger)
        map2shape.PCR_catch2shape(catchpointmap, catchsurfmap, proj4string_src, proj4string_trg, shapedata, logger)
        # convert raster to shape
        gdalCommand = 'gdal_polygonize.py ' + catchsurfmap + ' -f "ESRI Shapefile" ' + shapecatch + ' catchments CATCH'
        logger.info('Running GDAL with command: ' + gdalCommand)
        os.system(gdalCommand)
    else:
        logger.info("Skipping hydraulics...")
    # all is done,
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    del logger, ch
    print 'end'
        
    
if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], 'c:W:E:N:C:O:F:')

    main(opts)
        

