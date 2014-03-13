#!/usr/bin/env python
"""
 poly_density.py
 Calculate density of polygon data as a raster surface.
 Each raster cell contains a value indicating the cover of the underlying polygon (from 0 to 1).
 
 To get decent performance on large vector datasets, the input vector dataset must
 have a gdal-recognized spatial index (ie a .qix file for shapefiles as created by shtree)

 Usage:

 poly_density -S shapefile -E extent -C cellsize -o outputfile -F outputformat [-n] [-1]
 
 -n if set only the first feature encountered is used in the grid. Can be usefull
    if you have overlapping feutures
    
 -1 max value truncated to one (for overlapping polygons)
  

 Author: Matthew T. Perry, Adjusted by J. Schellekens, Deltares

 License: You are free to use or modify this code for any purpose.
          This license grants no warranty of any kind, express or implied. 
"""
import ogr
import sys
import numpy
import gdal
import getopt
import pdb

def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)

def getOpts():
    poly_ds = "lu_water.shp"
    poly_lyr = 0
    extent=[5.833,51.939,5.999,52.017]
    #extent = [-180., -90., 180., 90.]
    cellsize = 0.003 
    outfile = "lu_water.tif"
    format = "GTiff"
    return (poly_ds,poly_lyr,extent,cellsize,outfile,format)    
 
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
# if __name__ == "__main__":
    # Get the inputs
    (poly_ds,poly_lyr,extent,cellsize,outfile,format) = getOpts()  
    onlyOne = False
    max1 = False

    try:
        opts, args = getopt.getopt(argv, 'S:E:o:C:F:hn1')
    except getopt.error, msg:
        usage(msg)
    for o, a in opts:
        if o == '-S': poly_ds = a
        if o == '-E': exec("extent = " + a) # ,globals(), globals()
        if o == '-C': cellsize = float(a)
        if o == '-o': outfile = a
        if o == '-F': format = a
        if o == '-n': onlyOne=True
        if o == '-1': max1=True
        if o == '-h': 
		usage()
		return()

    # Get the input layer
    ds = ogr.Open(poly_ds)
    lyr = ds.GetLayer(poly_lyr)
    
    # TODO: Confirm dataset is polygon and extents overlap 

    ydist = extent[3] - extent[1]
    xdist = extent[2] - extent[0]
    
    #xcount = int(xdist/cellsize) + 1
    #ycount = int(ydist/cellsize) + 1
    xcount = int(round(xdist/cellsize))
    ycount = int(round(ydist/cellsize))


    # Create output raster  
    driver = gdal.GetDriverByName( format )
    dst_ds = driver.Create( outfile, xcount, ycount, 1, gdal.GDT_Float32 )

    # the GT(2) and GT(4) coefficients are zero,     
    # and the GT(1) is pixel width, and GT(5) is pixel height.     
    # The (GT(0),GT(3)) position is the top left corner of the top left pixel
    gt = (extent[0],cellsize,0,extent[3],0,(cellsize*-1.))
    dst_ds.SetGeoTransform(gt)
    
    dst_band = dst_ds.GetRasterBand(1)
    dst_band.SetNoDataValue( -9999 )

    pixelnum = 0
    
    for ypos in range(ycount):
        # Create output line array
        outArray = numpy.zeros( (1, xcount) )
        for xpos in range(xcount):
            # create a 4-item list of extents 
            minx = xpos * cellsize + extent[0] 
            maxy = extent[3] - ypos * cellsize 
            miny = maxy - cellsize
            maxx = minx + cellsize
            
            # Create Polygon geometry from BBOX
            wkt = 'POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))' \
               % (minx, miny, minx, maxy, maxx, maxy, maxx, miny, minx, miny)
            g = ogr.CreateGeometryFromWkt(wkt)

            # Set spatial filter
            lyr.SetSpatialFilter(g)


            #print wkt, lyr.GetFeatureCount()
            
            #continue

            # Loop through all features/geometries w/in filter
            feat = lyr.GetNextFeature()
            area = 0
            while feat is not None:

                # Intersect with polygon lyr
                sg = feat.GetGeometryRef().Intersection(g)
                
                if sg:
                    if onlyOne: # Hack to only count one feature if you have overlapping shapes
                        if area < 1E-10:
                            area = area + sg.GetArea()
                    else:
                        area = area + sg.GetArea()
                    
                feat = lyr.GetNextFeature()
            
            lyr.ResetReading()

            # Calculate area of intersection
            pct_cover = area / (cellsize*cellsize)
            
            if max1:
                pct_cover = min(1.0,pct_cover)
                
            #Assign percent areal cover as value in line array
            numpy.put( outArray, xpos, pct_cover )

            pixelnum += 1
 
        print '%.2f pct complete' % (float(pixelnum)/(xcount*ycount) * 100.)
        dst_band.WriteArray(outArray,0,ypos)

if __name__ == '__main__':
    sys.exit(main()) 
