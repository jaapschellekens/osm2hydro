import os
import math
from numpy import *

"""

Steps:

- in the first step the script take one or multiple polygons
  to cut-out parts of the planet file (made using use ogr2poly.py to convert to osm poly)
- in the second step te splitter program us used to split the pbf in tiles
  for each polygon,
- next the osm2hydro program is run fro each time
- finally the tiles are merged using gdal_merge 

"""


outputdir="3-tiles_world"
#osmfile="planet-131009.osm.pbf"
# An o5m file is quicker to process
osmfile="~/world/planet.o5m"
# These are the seperate boundaries to run splitter on
boundaries=["tileselect_0.poly","tileselect_1.poly","tileselect_2.poly","tileselect_3.poly","tileselect_4.poly","tileselect_5.poly","tileselect_6.poly"]
#boundaries=["areas.poly"]
osmconvert = "~/src/osm2hydro/dist/linux/osmconvert"
osm2hydro = "python ~/src/osm2hydro/src/osm2hydro/osm2hydro_metres.py "
#osm2hydro = "python ~/src/osm2hydro/applications/osm2hydro/src/osm2hydro/osm2hydro.py "



def run_splitter(osmfile,boundary,outputdir,startnumber=1,predareas=None):
    """
    Calls splitter to break-up the planet
    - boundary - osmpoly file (make with ogr2poly)
    - predareas - areas.list file from previous run
    """
    # Should be put in the python sccript
    if predareas:
        os.system("java -Xmx5000m -jar ~/src/splitter/dist/splitter.jar  --split-file=" + predareas + " --output-dir=" + outputdir + " --description=\"OSM World splitter results\" --keep-complete=true --mapid=" + str(startnumber) + " --output=pbf  --max-nodes=3200000 --polygon-file=" +boundary+ " --write-kml=areas.kml " + osmfile)
    else:
        os.system("java -Xmx5000m -jar ~/src/splitter/dist/splitter.jar --overlap=5000  --output-dir=" + outputdir + " --description=\"OSM World splitter results\" --keep-complete=true --mapid=" + str(startnumber) + " --output=pbf  --max-nodes=3200000  --polygon-file=" +boundary+ " --write-kml=areas.kml " + osmfile)


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





#----------------------------------

start=1 # to keep track of the number of tiles between splitter runs


# Loop over all regions (define by the boundary polygons)
for poly in boundaries:

    # Use osmconvert to cut out a pbf for the polygon if is not already esists
    thisodir = outputdir + "_" + poly
    if not os.path.exists(poly+".pbf"):
        print "Cutting " + osmfile + " using polygon : " + poly
        os.system(osmconvert + " " + osmfile+" -o="+poly+".pbf -B=" + poly + " -v  --out-pbf --hash-memory=1000 --drop-broken-refs  --drop-version --drop-relations --complex-ways --complete-ways --drop-author")
    else:
        print "Skipping extracting polygon cutout from planet: " + poly+ ".pbf"

    # Run the splitter program for this polygon. This will split it up in smaller tiles
    if not os.path.exists( os.path.join(thisodir,"areas.list")):
        run_splitter(poly+".pbf",poly,thisodir,startnumber=start,predareas=None)
    else:
        print "Skipping splitting of file: " + poly+ ".pbf"

    # read the tile boundaries from the splitter output (areas.list).
    # We need this to call osm2hydro for all the tiles that have been constructed
    res = read_splitter_areas(os.path.join(thisodir,"areas.list"))
    start=int(res[-1][0]) + 1
    print "Last tile number: " + str(start)

    # now run OSM2hydro for each tile
    osm2hydrocmd=[]
    for tile in res:
        tilefile = tile[0]
        if not os.path.exists(os.path.join(thisodir,str(tilefile),"OMS2Hydro.log")):
            tilefile = tile[0]
            ymin = tile[1]
            xmin = tile[2]
            ymax = tile[3]
            xmax = tile[4]
            osm2hydrostr = osm2hydro + " -N " + tilefile + " -C " + thisodir+ " -c planet_metres.ini -E \'[" + str(xmin) +"," + str(ymin) + "," + str(xmax) + "," + str(ymax) +"]\' -O  " + thisodir+   "/" + str(tilefile) + ".osm.pbf"
            print "Starting: " + osm2hydrostr
            os.system(osm2hydrostr)
        else:
            print "Skipping osm2hydro for tile: " + str(tile)
            

        
            
# Now collect the tiles and merge

mergefiles = ["roads_den.tif","lu_paved.tif","lu_pavedpol.tif","lu_water.tif","lu_unpaved.tif","lu_roads.tif"]


for mf in mergefiles:
    print "merging " + mf + "....."
    execstr = "gdal_merge.py -co 'COMPRESS=PACKBITS' -o merge_" + mf + " `find . -name " + mf + "`"
    os.system(execstr)

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

