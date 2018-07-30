"""
osm2tiff.py - convert an osm|pbf file to a set of tiff using gdal_rasterize. Filtering with tags is controlled by a .ini file.

Usage:

::
  osm2tiff.py [-h][-c configfile] -O my_osm_file.osm|pbf -o outputdir

  -O osm_file - the osm|pbf file
  -o outputdir - the directory to store the output in
  -h - show this information
  -c configfile (default is osm2tiff.ini)
  -M maxproc - maximum number of gdal_rasterize processes to start (default = 1)



dependencies
~~~~~~~~~~~~

  - gdal_rasterize (tested with gdal >= 1.10.0)
  - the OSM_CONFIG_FILE environment var should be set and point to a valid gdal osmconf.ini file. The
    osm2hydro source code contains a version that should be used with this script.
  
ini file
~~~~~~~~

The ini file has a section for each shape file that will be created:

  ::
    [raster]
    resolution=0.01

    [lu_roads_sec]
    type=lines
    highway=motorway_link,secondary,"road"
    aeroway='runway','taxiway'
        
    [lu_roads_main]
    type=lines 
    highway='motorway','trunk','primary'
        
        
    [lu_roads_small]
    type=lines
    highway='footway','path','pedestrian','residential','service','tertiary','track','unclassified','cycleway'
        
    [lu_water]
    type=multipolygons
    natural='water','reservoir','basin','salt_pond'
    landuse='water','basin','reservoir','salt_pond'
       
      
$Author: $
$Id: $
$Rev: $
"""

"""
 Copyright notice

    Copyright 2018 Jaap Schellekens, VanderSat

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
"""

import os,sys,getopt
import platform
import ConfigParser
import subprocess
import shlex
import time
import gdal
import numpy as np


def readMap(fileName, fileFormat):
    """
    Read geographical file into memory

    Input:
        filename - string holding file name
        fileFormat - a file format string according to the GDAL list of formats

    Returns:
        x,y,data,FillVal

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
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols)
    y = np.linspace(originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)  # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    FillVal = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    return x, y, data, FillVal


def writeMap(fileName, fileFormat, x, y, data, FillVal):
    """ Write geographical data into file"""

    verbose = False
    srs = None
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(fileFormat)

    # Processing
    if verbose:
        print 'Writing to temporary file ' + fileName + '.tif'
    # Create Output filename from (FEWS) product name and date and open for writing
    TempDataset = driver1.Create(fileName + '.tif', data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    # Give georeferences
    xul = x[0] - (x[1] - x[0]) / 2
    yul = y[0] + (y[0] - y[1]) / 2
    TempDataset.SetGeoTransform([xul, x[1] - x[0], 0, yul, 0, y[1] - y[0]])
    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    # This seesm to happen somtimes
    if FillVal == None:
        FillVal = 1E31

    TempBand.SetNoDataValue(FillVal)
    # Create data to write to correct format (supported by 'CreateCopy')
    if verbose:
        print 'Writing to ' + fileName + '.map'
    outDataset = driver2.CreateCopy(fileName, TempDataset, 0)
    #outDataset.SetProjection(srs)

    TempDataset = None
    outDataset = None
    if verbose:
        print 'Removing temporary file ' + fileName + '.tif'
    os.remove(fileName + '.tif');

    if verbose:
        print 'Writing to ' + fileName + ' is done!'


def configget(config,section,var,default):
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
        print( "returning default (" + default + ") for " + section + ":" + var)
        configset(config,section,var,default, overwrite=False)
    
    default = Def
    return ret       




def configset(config,section,var,value, overwrite=False):
    """   
    Sets a string in the in memory representation of the config object
    Does NOT overwrite existing values if overwrite is set to False (default)
    
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

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config


def mkgdal_rasterizestr(config,sec,resolution=0.008):
    """
    Create the gdal_rasterize string from a ini file to select data using 
    sql statements
    - one statement for each section
    - multiple tags per feature using a comma separated list

    :returns list of files created:
    """

    thestr = ""
    nropt = 0
    stype = config.get(sec,'type')
    for opt in config.options(sec):
        if opt == 'type':
            stype = config.get(sec,opt) 
        elif opt == 'width':
            metresdegree = 110570.01
            burnval = metresdegree*resolution/float(config.get(sec,opt))
        else:
            tags = config.get(sec,opt).split(',')
            for i in tags:
                nropt = nropt +1
                tag = i.strip('"').strip("'")
                if nropt == 1:
                    thestr = thestr + opt +"=" + "'" + tag + "'"
                else:
                    thestr = thestr + " or " + opt +"=" + "'" + tag + "'"

    if stype == 'lines':
        thestr = "-burn " + str(burnval) + "  -sql \"select * from lines where " + thestr
    else: 
        if stype=='multipolygons':
            thestr = " -burn 1 -sql \"select * from multipolygons where " + thestr
        else:
                if stype=='other_relations':
                    thestr = " -burn 1 -sql \"select * from other_relations where " + thestr
                elif stype == 'multilinestrings':
                    thestr = " -burn 1 -sql \"select * from multilinestrings where " + thestr
                else:
                     print "unexpected type in ini file section " + sec + "(" + stype +")"
                     return None
    
    thestr = thestr + "\""
            
    return thestr
    
    
def extractlayers_new(config,osmfile,outputdir,maxprocesses=2,gdal_rasterize="gdal_rasterize"):
    """
    http://stackoverflow.com/questions/4992400/running-several-system-commands-in-parallel         
    """
        
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
            command = command.replace('\\','/') # otherwise shlex.split removes all path separators
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
        print "All gdal_rasterize processes (" + str(len(commands)) + ") completed."
        
            
    commands = []
    filesmade = []

    """ Create the commands here... """    
    res = float(configget(config,'raster','resolution','0.1'))
    for sec in config.sections():
        print sec
        if 'raster' not in sec:
            thestr = mkgdal_rasterizestr(config,sec,resolution=res)
            # check if the tiff file already exists
            filesmade.append(os.path.join(outputdir, sec + ".tif"))
            if os.path.exists(outputdir + "/" + sec + ".tif"):
                print "Removing existing tif: " + sec
                os.remove(outputdir + "/" + sec + ".tif")
            thestr = gdal_rasterize + " -at -add -co TILED=TRUE -co COMPRESS=LZW -tr " + str(res) +" " + str(res) + " " + osmfile + " " +  outputdir + "/" + sec + ".tif " + thestr

            commands.append(thestr)


    runCommands(commands,maxprocesses)

    return filesmade



def mergelayers(allfiles):
    """
    Merges the different water component layers into one map with water coverage maximizing at 1

    :return: one map
    """

    for thisfile in allfiles:
        x, y, thistif, fillval = readMap(thisfile,'GTiff')
        # set fillval to nan
        thistif[thistif == fillval] = np.nan
        if thisfile == allfiles[0]:
            outwater = thistif.copy()
        else:
            outwater = outwater + thistif

    outwater[outwater > 1.0] = 1.0

    #reset nan to fillval
    outwater[np.logical_not(np.isfinite(outwater))] = fillval

    return x, y, outwater, fillval






def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

def main(argv=None):
    """
    Main function, processes the command-line options, reads the config file
    and starts the extract function.
    """
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return    
    rootLoc = os.path.split(sys.argv[0])[0]

    # First try to find gdal_rasterize in a dist directory (part of the osm2hydro package). If
    # not found assume it is located in the search path
    if platform.system() == "Linux" or platform.system() == "Linux2":
        gdal_rasterize = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','gdal_rasterize'))
    else:
        gdal_rasterize = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','gdal_rasterize.exe'))

    if not os.path.exists(gdal_rasterize):
        gdal_rasterize = "/usr/bin/gdal_rasterize"
 
    osmfile= "not_set.osm"
    outputdir= "./"
    configfile = 'osm2tiff.ini'
    maxcpu = 1
   
    try:
        opts, args = getopt.getopt(argv, 'O:o:hc:M:')
    except getopt.error, msg:
        usage(msg)
        
    for o, a in opts:
        if o == '-O': osmfile = a
        if o == '-M': maxcpu  = int(a)
        if o == '-o': 
            outputdir = a
        if o == '-h':
            usage()
            return()

    # First create the directories if they do not already exists
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)    

    # Set interleaved reading
    os.environ['OGR_INTERLEAVED_READING']='YES'

    if not os.path.exists(os.environ['OSM_CONFIG_FILE']):
        print "Gdal osm config cannot be found: " + os.environ['OSM_CONFIG_FILE']
        print "... or the OSM_CONFIG_FILE environment variable is not set"
        
    config = iniFileSetUp(configfile)
    
    filesmade = extractlayers_new(config,osmfile,outputdir,maxprocesses=maxcpu,gdal_rasterize=gdal_rasterize)
        
    x, y, finalout, fillval = mergelayers(filesmade)

    writeMap(os.path.join(outputdir,'merged.tif'),'GTiff',x,y,finalout,fillval)



        
    
if __name__ == '__main__':
    sys.exit(main())
