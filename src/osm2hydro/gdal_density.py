"""
Created on Fri Dec  7 15:00:29 2012

Usage:

:: 

 gdal_density -S shapefile -E extent -C cellsize -o outputfile 
              -F outputformat [-r resamplefactor (default = 10)]
              [-b burnvalue (default=1)][-G True|False][-D True|False]
              [-t tempdir][-M]
 
 -G Use gdalwarp instead of pcraster resample (default= False)
 -D delete high resolution files after processing (Default=True)
 -E [a,b,c,d] following gdal conventions
 -t tempdir (to store the high resolution temporary files, default is the current dir
  or the GDAL_DENSITY_TMP environment variable)
 -M if specified the burn value is assumed to be given in metres
 
Converts a shape to a grid. The resulting grid holds (for each cell) the
fraction (0-1) covered by the shapefile. The resamplefactor determines the 
accuracy of the final results. If r is one you get a maps with only zeros and
ones (two possibilities). If r is two you get 5 possible fractions 
(0,0.25,0.5,0.75,1), if r is ten you get a map with 101 possible fractions 
etc...


requirements:

 - gdal_warp
 - gdal_rasterize
 - pcraster 4.0 + python bindings
 
 
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
import os
import getopt
import sys
import ogr
import gdal
import pdb
from pcraster import *
import osgeo.osr as osr
import osgeo.gdal as gdal
import math
import tempfile

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)



def lattometres(lat):
    """"
    Determines the length of one degree lat/long at a given latitude (in meter).
    Code taken from http:www.nga.mil/MSISiteContent/StaticFiles/Calculators/degree.html
    Input: map with lattitude values for each cell
    Returns: length of a cell lat, length of a cell long
    """
    #radlat = spatial(lat * ((2.0 * math.pi)/360.0))
    #radlat = lat * (2.0 * math.pi)/360.0
    radlat = spatial(lat) # pcraster cos/sin work in degrees!
    
    
    m1 = 111132.92        # latitude calculation term 1
    m2 = -559.82        # latitude calculation term 2
    m3 = 1.175            # latitude calculation term 3
    m4 = -0.0023        # latitude calculation term 4
    p1 = 111412.84        # longitude calculation term 1
    p2 = -93.5            # longitude calculation term 2
    p3 = 0.118            # longitude calculation term 3
    # # Calculate the length of a degree of latitude and longitude in meters
    
    latlen = m1 + (m2 * cos(2.0 * radlat)) + (m3 * cos(4.0 * radlat)) + (m4 * cos(6.0 * radlat))
    longlen = (p1 * cos(radlat)) + (p2 * cos(3.0 * radlat)) + (p3 * cos(5.0 * radlat))
        
    return latlen, longlen  
    


def detDegreeLen(metres,map):
    """
    in length in metres
    out length in degree
    """
    
    aa = ycoordinate(boolean(map + 1.0))
    yl, xl = lattometres(aa)
       
    xl = metres/xl
    yl = metres/yl
    
    #xl = xl * celllength()
    #yl = yl * celllength()
    
    reallength = (xl + yl) * 0.5
    
    return reallength
     

def makeMultMap(outfile,pcrout,metres,gdal_translate="gdal_translate"):
    """
    Make a map in metres to multiply with and perform the multiplication
    """
    
     # Metres
    os.system(gdal_translate + " -of PCRaster " + outfile + " " + pcrout)
    setclone(pcrout)
    orgmap = readmap(pcrout)
    reallength = detDegreeLen(metres,orgmap)

    ttmap = os.path.join(os.path.dirname(pcrout),"_tt.map")
    report(orgmap * reallength/celllength(),ttmap)
    os.system(gdal_translate + " -of PCRaster " + ttmap + " " + outfile)
    
    
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return    
    

    gdal_rasterize = "gdal_rasterize"
    gdal_translate ="gdal_translate"
    gdal_warp = "gdalwarp"
    delete_hires=True
    rsamp = 10.0
    burn=1.0
    resamplewithgdal = False
    format = "GTiff"
    verbose = True
    burninmetres = False
   
    try:
        opts, args = getopt.getopt(argv, 'G:S:E:o:C:F:hr:b:D:t:M')
    except getopt.error, msg:
        usage(msg)
    
    try:
      tt = os.environ['GDAL_DENSITY_TMP']
      tmpdir = tt
    except:
      tmpdir = '.'
    
    for o, a in opts:
        if o == '-S': poly_ds = a
        if o == '-t': tmpdir = a
        if o == '-E': exec("extent = " + a) # ,globals(), globals()
        if o == '-G': exec("resamplewithgdal = " + a) 
        if o == '-D': exec("delete_hires = " + a)
        if o == '-M': burninmetres = True
        if o == '-C': 
            cellsize = float(a)
        if o == '-r': rsamp = float(a)
        if o == '-b': burn = float(a)
        if o == '-o':
            outfile = a           
        if o == '-F': format = a
        if o == '-h':
            usage()
            return()
    
    if tmpdir is not ".":
      tmp_name=tempfile.mkstemp(dir=tmpdir)
      outfilehires = tmp_name[1]  + "_hires.tif"
    else:
      tmp_name=tempfile.mkstemp(dir=tmpdir)
      outfilehires = tmp_name[1] +  "_hires.tif"

      
    print "gdal_density starting with options: " + str(argv)
    print "tmpfile: " + tmp_name[1]

    if not os.path.exists(poly_ds):
        return

    # Open shape to check if it is not empty
    ds = ogr.Open(poly_ds)
    lyr = ds.GetLayer(0)    
 
    #lyr = 1    
    hirescellsize = cellsize/rsamp
    pcroutfilehires= outfilehires + ".map"
    pcroutfile=outfile + ".map"
    width = (extent[2]-extent[0])/cellsize
    height = (extent[3]-extent[1])/cellsize
    width = int(math.ceil(width))
    height = int(math.ceil(height))

    if lyr.GetFeatureCount() > 0:

        if burninmetres:
            exestr = gdal_rasterize + " -tr " + str(hirescellsize) + " " + str(hirescellsize) + " -te " + str(extent[0]) + " " + str(extent[1]) + " " + str(extent[2]) + " " + str(extent[3]) + " -burn " + str(1.0 * rsamp) + " " + poly_ds + " " + outfilehires
            os.system(exestr)
        else:
            exestr = gdal_rasterize + " -tr " + str(hirescellsize) + " " + str(hirescellsize) + " -te " + str(extent[0]) + " " + str(extent[1]) + " " + str(extent[2]) + " " + str(extent[3]) + " -burn " + str(burn) + " " + poly_ds + " " + outfilehires
            
            if verbose:
                print "starting: " + exestr
            
            os.system(exestr)
        
        if resamplewithgdal:
            execstr = gdal_warp + " -ot Float32 -r average -ts " + str(width) + " " + str(height) + " -te " + str(extent[0]) + " " + str(extent[1]) + " " + str(extent[2]) + " " + str(extent[3]) + " " + outfilehires + ' ' + outfile
            os.execstr_org = gdal_warp + " -ot Float32 -r average -tr "+ str(cellsize) + " " + str(cellsize) + " " + outfilehires + ' ' + outfile

            os.system(execstr)
            if burninmetres:
                makeMultMap(outfile,"_temp.map",burn,gdal_translate="gdal_translate")
        else:
            os.system(gdal_translate + " -of PCRaster " + outfilehires + " " + pcroutfilehires)
            os.system("resample -e 20 -r "+ str(rsamp) + " " + pcroutfilehires + " " + pcroutfile)
            os.system(gdal_translate + " -of " + format + " " + pcroutfile + " " + outfile)
            if burninmetres:
                makeMultMap(outfile,"_temp.map",burn,gdal_translate="gdal_translate")

        
        if delete_hires:
            os.remove(outfilehires)
    else: # create an empty file9
        geotransform = (extent[0], hirescellsize, 0.0, extent[1], 0.0, -hirescellsize)
        raster = gdal.GetDriverByName('GTiff')
        ds2 = raster.Create(outfile, width, height, 1, gdal.GDT_Float32)
        ds2.SetGeoTransform(geotransform)
        ds2.GetRasterBand(1).SetNoDataValue(0.0)

        srs = lyr.GetSpatialRef()
        ds2.SetProjection(srs.ExportToWkt())
        ds2 = None

    lyr = None
    ds = None # close

    
if __name__ == '__main__':
    sys.exit(main()) 
