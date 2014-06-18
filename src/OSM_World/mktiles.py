"""
Usage:

   mkworld_paved.py -I configfile


Steps:

Preperation
~~~~~~~~~~~
- in the first step the script take one or multiple polygons
  to cut-out parts of the planet file (made using use ogr2poly.py to convert to osm poly)
- make ini files fro osm2hydro for each polygon. The filenames
  should have the same name but a .ini extension
Running
~~~~~~~
- obtain a planet.osm file (use pbf or o5m format is possible)

Run the script:
- in the second step te splitter program us used to split the pbf in tiles
  for each polygon,
- next the osm2hydro program is run for each time
- finally the tiles are merged using gdal_merge 


Runtime
~~~~~~~
# Breaking up of the planet file in areas: +- 4hrs
# Splitting-up the areas in small tiles using splitter: +-
# Running osm2hydro on all tiles
# merging tiles into makes fro the world

"""

import ConfigParser, sys, os, getopt
import time, subprocess, shlex

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
        print "All processes (" + str(len(commands)) + ") completed."




def iniFileSetUp(configfile):
    """
    Reads .ini file and sets default values if not present
    """

    #setTheEnv(runId='runId,caseName='caseName)
    # Try and read config file and set default options
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(configfile)
    return config


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
        #log.info( "returning default (" + default + ") for " + section + ":" + var)
        configset(config,section,var,default, overwrite=False)

    default = Def
    return ret



def run_splitter(osmfile,boundary,outputdir,startnumber=1,predareas=None,java="c:\Windows\System32\java.exe"):
    """
    Calls splitter to break-up the planet
    - boundary - osmpoly file (make with ogr2poly)
    - predareas - areas.list file from previous run
    """
    # Should be put in the python sccript
    if predareas:
        ret = os.system(java + " -Xmx10000m -jar ../../third-party/splitter/splitter.jar  --split-file=" + predareas + " --output-dir=" + outputdir + " --description=\"OSM World splitter results\" --keep-complete=true --mapid=" + str(startnumber) + " --output=pbf  --max-nodes=3200000 --polygon-file=" +boundary+ " --write-kml=areas.kml " + osmfile)
    else:
        ret = os.system(java + " -Xmx10000m -jar ../../third-party/splitter/splitter.jar --overlap=5000  --output-dir=" + outputdir + " --description=\"OSM World splitter results\" --keep-complete=true --mapid=" + str(startnumber) + " --output=pbf  --max-nodes=3200000  --polygon-file=" +boundary+ " --write-kml=areas.kml " + osmfile)

    if ret:
        print "Error, command returned non zero exit code:" + str(ret)
        exit(ret)


def read_splitter_areas(fname):
    """
    reads a splitter .list file with tile min max
    returns list of lists:
        tilename,ymin,xmin,ymax,xmax
        
    Example file:
    -------------
    # List of areas
    # Generated Thu Nov 07 15:12:37 CET 2013
    #
    00000001: -1933312,5134336 to -538624,5818368
    #       : -41.484375,110.170898 to -11.557617,124.848633
    
    00000002: -985088,5818368 to -387072,7188480
    #       : -21.137695,124.848633 to -8.305664,154.248047
    
    """
    def my_split(s, seps):
        res = [s]
        for sep in seps:
            s, res = res, []
            for seq in s:
                res += seq.split(sep)
        return res
    
    ret = []
    fs = open(fname,'r')
    # Convert the splitter map units back to decimal degree
    corconv = 360.0/2.0**24
    
    for line in fs:
        if '#' not in line:  
            zz = my_split(line.strip(),[':',',',' to '])
            if len(zz) == 5:
                new = []
                for it in zz:
                    new.append(it.strip())
                new =[]
                new.append(zz[0])
                for it in zz[1:len(zz)]:
                    new.append(int(it) * corconv)
                    
                ret.append(new)    
    fs.close()   
    return ret


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)




def main(argv=None):


    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    opts, args = getopt.getopt(sys.argv[1:], 'I:')
    configfile = "planet.ini"

    for o, a in opts:
        print o, a
        if o == '-I': configfile = a

    config = iniFileSetUp(configfile)

    java = configget(None,config,"osmworld","java","c:\Windows\System32\java.exe")
    outputdir = configget(None,config,"osmworld","outputdir","world")
    osmfile = configget(None,config,"osmworld","osmfile","planet.o5m")
    _boundaries = configget(None,config,"osmworld","boundaries",'["area_0.poly","area_1.poly","area_2.poly","area_3.poly","area_4.poly","area_5.poly","area_6.poly"]')
    exec "boundaries = " + _boundaries
    osmconvert = configget(None,config,"osmworld","osmconvert","..\..\dist\win32\osmconvert.exe")
    osm2hydro = configget(None,config,"osmworld","osm2hydro", "python ..\osm2hydro\osm2hydro_metres.py ")


    # Main loop
    #----------------------------------
    start=1 # to keep track of the number of tiles between splitter runs
    # Loop over all regions (define by the boundary polygons)
    for poly in boundaries:
        # Use osmconvert to cut out a pbf for the polygon if is not already esists
        thisodir = outputdir + "_" + poly
        if not os.path.exists(poly+".pbf"):
            print "Cutting " + osmfile + " using polygon : " + poly
            execstr = osmconvert + " " + osmfile +" -o="+poly+".pbf -B=" + poly + " -v  --out-pbf --hash-memory=1000 --drop-broken-refs  --drop-version --drop-relations --complex-ways --complete-ways --drop-author"
            ret = os.system(execstr)
            if ret:
                print "Error, command returned non zero exit code:" + str(ret)
                exit(ret)
        else:
            print "Skipping extracting polygon cutout from planet: " + poly+ ".pbf"

    for poly in boundaries:
        thisodir = outputdir + "_" + poly
        # Run the splitter program for this polygon. This will split it up in smaller tiles
        if not os.path.exists( os.path.join(thisodir,"areas.list")):
            run_splitter(poly+".pbf",poly,thisodir,startnumber=start,predareas=None)
        else:
            print "Skipping splitting of file: " + poly+ ".pbf"

        # read the tile boundaries from the splitter output (areas.list).
        # We need this to call osm2hydro for all the tiles that have been constructed
        res = read_splitter_areas(os.path.join(thisodir,"areas.list"))
        #start=int(res[-1][0]) + 1


    for poly in boundaries:
        thisodir = outputdir + "_" + poly
        # now run OSM2hydro for each tile
        osm2hydrocmd=[]
        res = read_splitter_areas(os.path.join(thisodir,"areas.list"))
        start=int(res[-1][0]) + 1
        print "Last tile number: " + str(start)
        osm2ini = os.path.splitext(os.path.basename(poly))[0] + ".ini"
        for tile in res:
            tilefile = tile[0]
            if not os.path.exists(os.path.join(thisodir,str(tilefile),"OMS2Hydro.log")):
                ymin = tile[1]
                xmin = tile[2]
                ymax = tile[3]
                xmax = tile[4]
                # This was for liux
                #osm2hydrostr = osm2hydro + " -N " + tilefile + " -C " + thisodir+ " -c " + osm2ini + "  -E \'[" + str(xmin) +"," + str(ymin) + "," + str(xmax) + "," + str(ymax) +"]\' -O  " + thisodir+   "/" + str(tilefile) + ".osm.pbf"
                # This for windows
                osm2hydrostr = osm2hydro + " -N " + tilefile + " -C " + thisodir+ " -c " + osm2ini + "  -E [" + str(xmin) +"," + str(ymin) + "," + str(xmax) + "," + str(ymax) +"] -O  " + thisodir+   "/" + str(tilefile) + ".osm.pbf"
                print "Starting: " + osm2hydrostr
                ret = os.system(osm2hydrostr)
                if ret:
                    print "Error, command returned non zero exit code:" + str(ret)
                    exit(ret)
            else:
                print "Skipping osm2hydro for tile: " + str(tile)



    exit()
    mergefiles = ["roads_den.tif","lu_paved.tif","lu_pavedpol.tif","lu_water.tif","lu_unpaved.tif","lu_roads.tif"]


    for mf in mergefiles:
        print "merging " + mf + "....."
        execstr = "gdal_merge.py -co 'COMPRESS=PACKBITS' -o merge_" + mf + " `find . -name " + mf + "`"
        ret = os.system(execstr)
        if ret:
            print "Error, command returned non zero exit code:" + str(ret)
            exit(ret)

    #CAN WE USE 0cutline in gdalworp??

    # Now make a catchment mask
    execstr ="pcrcalc mask.map = " + mergefiles[0] + "* 0.0"
    os.system(execstr)
    execstr="gdal_translate mask.map -co 'COMPRESS=PACKBITS' -ot Float32 mask.tif"
    os.system(execstr)
    execstr = "gdal_rasterize -burn 1 -l TM-WORLD_BORDERS-0.3 TM-WORLD_BORDERS-0.3.shp mask.tif "
    os.system(execstr)


    for mf in mergefiles:
        execstr = "gdal_translate.py  merge_" + mf + " -of PCRaster merge_" + mf +".map"
        os.system(execstr)
        execstr = "resample -r10 merge_" + mf + ".map small_merge_" + mf + ".map"
        os.system(execstr)

if __name__ == "__main__":
    main()