"""
osm2shp - convert an osm|pbf file to a set of shapefiles using ogr2ogr. Filtering with tags is controlled by a .ini file.

Usage:

::

  osm2shp.py [-h][-c configfile] -O my_osm_file.osm|pbf -o outputdir

  -O osm_file - the osm|pbf file
  -o outputdir - the directory to store the output in
  -h - show this information
  -c configfile (default is osm2shp.ini)
  -M maxproc - maximum number of ogr2ogr processes to start (default = 2)


dependencies
~~~~~~~~~~~~

  - ogr2ogr (tested with gdal >= 1.10.0)
  - the OSM_CONFIG_FILE environment var should be set and point to a valid gdal osmconf.ini file
  
ini file
~~~~~~~~

The ini file has a section for each shape file created:

  ::

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
"""

import os,sys,getopt
import platform
import ConfigParser
import subprocess
import shlex
import time


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


def mkogr2ogrstr(config,sec):
    """
    Create the ogr2ogr string from a ini file to select data using 
    sql statements
    - one statement for each section
    - multiple tags per feature using a comma separated list
    """
    thestr = ""
    nropt = 0
    stype = config.get(sec,'type')
    for opt in config.options(sec):
        if opt == 'type':
            stype = config.get(sec,opt) 
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
        thestr = "-lco SHPT=ARC -sql \"select * from lines where " + thestr
    else: 
        if stype=='multipolygons':
            thestr = "-lco SHPT=POLYGON -sql \"select * from multipolygons where " + thestr
        else:
            print "unexpected type in ini file section " + sec + "(" + stype +")"
            return None
    
    thestr = thestr + "\""
            
    return thestr
    
    
def extractlayers_new(config,osmfile,outputdir,maxprocesses=2,ogr2ogr="ogr2ogr"):
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
        print "All ogr2ogr processes (" + str(len(commands)) + ") completed."
        
            
    commands = []

    """ Create the commands here... """    
    for sec in config.sections():
        thestr = mkogr2ogrstr(config,sec)
        # check if the shpe file already exists
        if os.path.exists(outputdir + "/" + sec + ".shp"):
            print "Removing existing shape: " + sec
            os.remove(outputdir + "/" + sec + ".shp")
            os.remove(outputdir + "/" + sec + ".dbf")
            os.remove(outputdir + "/" + sec + ".prj")
            os.remove(outputdir + "/" + sec + ".shx")
        thestr = ogr2ogr + " -f \"ESRI Shapefile\"  --config OGR_INTERLEAVED_READING YES " + outputdir + "/" + sec + ".shp " + osmfile + " " + thestr

        commands.append(thestr)


    runCommands(commands,maxprocesses)


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
    rootLoc        = os.path.split(sys.argv[0])[0]
    
    if platform.system() == "Linux" or platform.system() == "Linux2":
        ogr2ogr = os.path.abspath(os.path.join(rootLoc,'..','..','dist','linux','ogr2ogr'))
    else:
        ogr2ogr = os.path.abspath(os.path.join(rootLoc,'..','..','dist','win32','ogr2ogr.exe'))

    if not os.path.exists(ogr2ogr):
        ogr2ogr = "ogr2ogr"   
 
    osmfile= "not_set.osm"
    outputdir= "./"
    configfile = 'osm2shp.ini'
    maxcpu = 2
   
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

    os.environ['OGR_INTERLEAVED_READING']='YES'

    if not os.path.exists(os.environ['OSM_CONFIG_FILE']):
        print "Gdal osm config cannot be found: " + os.environ['OSM_CONFIG_FILE']
        
    config = iniFileSetUp(configfile)
    
    extractlayers_new(config,osmfile,outputdir,maxprocesses=maxcpu,ogr2ogr=ogr2ogr)
        
    return

        
    
if __name__ == '__main__':
    sys.exit(main())