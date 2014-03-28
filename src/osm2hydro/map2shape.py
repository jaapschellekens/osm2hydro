#!/usr/bin/env python
"""
This script converts gridded river sections into a shape file. 

Input:

    rivermap:       PCRaster map, containing the river sections as integers
    drainmap:       Accumulated drainage map over the pixels, contained in rivermap.__add__
    
Both inputs are produced by a script called 'srtm_burn_process.bat' (no shell-version available yet)
This script can be run from within srtm_burn_process.bat eventually
"""

"""
 Copyright notice

    Copyright 2013, 2014 Hesssel Winsemiuss

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

import csv
from datetime import datetime
import numpy as np
import sys
import scipy
import osgeo.ogr as ogr
import osgeo.gdal as gdal
import osr
import pyproj
import pdb
import os

def readMap(fileName, fileFormat):
    """ 
    Read geographical file into memory
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


def convertCoord(proj_src, proj_trg, x, y):
    """
    Convert a list of x,y pairs in a certain projection to another projection

    input:

        proj_src:   string, EPSG or proj4 string referring to projection of source coordinates
        proj_trg:   string, EPSG or proj4 string referring to projection of target coordinates
        x:          NumPy array, vector or 2D array of x-coordinates (source)
        y:          NumPy array, vector or 2D array of y-coordinates (source)

    output:

        X:          NumPy array, vector or 2D array of x-coordinates (target)
        Y:          NumPy array, vector or 2D array of y-coordinates (target)
    """
    srs1 = pyproj.Proj(proj_src) # OPT['proj4_params'])
    srs2 = pyproj.Proj(proj_trg) # wgs84
    X,Y  = pyproj.transform(srs1, srs2, x,y) # Do add 0. to avoid trunc issues.
    return X,Y

def degree_grid_to_surface(R, lati, loni):
    # retrieve y-axis
    grid=True
    if len(lati.shape) == 2:
        lat = lati[:,0]
        lon = loni[0,:]
    elif len(lati.shape) == 1:
        grid=False
        lat = lati
        lon = loni
    else:
        print 'Number of dimensions in y-axis is incorrect, must be 1 or 2'
        #break
    # free some memory
    del lati, loni
    lon_res = (lon.max()-lon.min())/(len(lon)-1)
    lat_res = (lat.max()-lat.min())/(len(lat)-1)
    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0
    xdist_per_deg = (degrees_to_radians*np.cos(lat*degrees_to_radians))*R
    ydist_per_deg = degrees_to_radians*R
    xdist_per_cell = xdist_per_deg*lon_res
    ydist_per_cell = ydist_per_deg*lat_res
    
    # compute over y-axis, the surface per grid cell
    surf_per_cell  = xdist_per_cell*ydist_per_cell
    # spread out over full grid
    if grid:
        surf_reshape   = surf_per_cell.reshape(len(surf_per_cell),1)
        surf           = surf_reshape.repeat(len(lon),axis=1)
    else:
        surf           = surf_per_cell
        
    return surf


def PCR_catch2shape(catchpointmap, catchsurfmap, proj_src, proj_trg, SHP_FILENAME, logger):
    # prepare osr object for writing .prj file
    outSpatialRef_latlon = osr.SpatialReference()    
    outSpatialRef_latlon.ImportFromProj4(proj_src)
    outSpatialRef_latlon.MorphToESRI()
    
    outSpatialRef = osr.SpatialReference()    
    outSpatialRef.ImportFromProj4(proj_trg)
    outSpatialRef.MorphToESRI()
    logger.info('Reading ' + catchpointmap)
    x, y, catchpoints, FillVal = readMap(catchpointmap, 'PCRaster');catchpoints[catchpoints==FillVal]= -1
    logger.info('Reading ' + catchsurfmap)
    x, y, catchsurf, FillVal = readMap(catchsurfmap, 'PCRaster');catchsurf[catchsurf==FillVal]= -1
    xi,yi = np.meshgrid(x,y)
    
    # compute surface of the cells, only compute the values over one column of the y-values, all the other are the same!!!
    surf = degree_grid_to_surface(6371e3, yi[:,0], xi[0,:])
    #pdb.set_trace()
    [iiy,iix] = np.where(catchpoints > 0)
    
#    longitudes = xi[iiy,iix]
#    latitudes  = yi[iiy,iix]
    catchPoint    = catchpoints[iiy,iix]
#    drainage   = drain[iiy,iix]
#    depths     = depth[iiy,iix]
#    widths     = width[iiy,iix]
    maxCatchPoint = catchPoint.max()
    
    # Create new shapefile
    ogr.UseExceptions()
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(SHP_FILENAME)
    layer_point = ds.CreateLayer("laterals", None, ogr.wkbPoint)
    layer_point_latlon = ds.CreateLayer("laterals_latlon", None, ogr.wkbPoint)
    point_ID = ogr.FieldDefn()
    point_ID.SetName('LATERAL')
    point_ID.SetType(ogr.OFTString)
    point_ID.SetWidth(15)
    layer_point.CreateField(point_ID)
    layer_point_latlon.CreateField(point_ID)
    
    point_ID = ogr.FieldDefn()
    point_ID.SetName('SURFACE')
    point_ID.SetType(ogr.OFTString)
    point_ID.SetWidth(15)
    layer_point.CreateField(point_ID)
    layer_point_latlon.CreateField(point_ID)
    
    # Create a new line geometry per river segment
    for id in np.arange(1, maxCatchPoint + 1):
        logger.info('Writing point element "' + str(id) + '" of "' + str(maxCatchPoint) + '"')
        
        y_idx, x_idx = np.where(catchpoints == id)
        lat_select = yi[y_idx, x_idx]
        lon_select = xi[y_idx, x_idx]
        # convert lat-lon coordinates to projected coordinates
        x_proj_select, y_proj_select = convertCoord(proj_src, proj_trg, lon_select, lat_select)
        # compute the areal surface of the found lateral/catchment
        y_surf, x_surf = np.where(catchsurf == id)
        #pdb.set_trace()
        catcharea = surf[y_surf].sum()
        # make point features and write to file
        point_latlon = ogr.Geometry(type=ogr.wkbPoint)
        point_latlon.SetPoint(0, np.float64(lon_select), np.float64(lat_select))
        feature_point_latlon = ogr.Feature(layer_point_latlon.GetLayerDefn())
        feature_point_latlon.SetGeometryDirectly(point_latlon)
        # write the attributes to the current feature
        #pdb.set_trace()
        feature_point_latlon.SetField('LATERAL', 'LAT' + str(id))
        feature_point_latlon.SetField('SURFACE', str(np.round(catcharea)))
        layer_point_latlon.CreateFeature(feature_point_latlon)
        feature_point_latlon.Destroy()

        point = ogr.Geometry(type=ogr.wkbPoint)
        #pdb.set_trace()
        point.SetPoint(0, np.float64(x_proj_select[0]), np.float64(y_proj_select[0]))
        feature_point = ogr.Feature(layer_point.GetLayerDefn())
        feature_point.SetGeometryDirectly(point)
        # write the attributes to the current feature
        feature_point.SetField('LATERAL', 'LAT' + str(id))
        feature_point.SetField('SURFACE', str(np.round(catcharea)))
        layer_point.CreateFeature(feature_point)
        feature_point.Destroy()

    ds.Destroy()
    # finally write .prj files
    catchpoint_prj = os.path.join(SHP_FILENAME, 'laterals.prj')
    catchpoint_latlon_prj = os.path.join(SHP_FILENAME, 'laterals_latlon.prj')
    catchpoint_prj_file = open(catchpoint_prj, 'w');catchpoint_prj_file.write(outSpatialRef.ExportToWkt());catchpoint_prj_file.close()
    catchpoint_latlon_prj_file = open(catchpoint_latlon_prj, 'w');catchpoint_latlon_prj_file.write(outSpatialRef.ExportToWkt());catchpoint_latlon_prj_file.close()

def PCR_river2Shape(rivermap, drainmap, lddmap, depthmap, widthmap, demmap, proj_src, proj_trg, SHP_FILENAME, logger):
    # prepare osr object for writing .prj file
    outSpatialRef_latlon = osr.SpatialReference()    
    outSpatialRef_latlon.ImportFromProj4(proj_src)
    outSpatialRef_latlon.MorphToESRI()
    
    outSpatialRef = osr.SpatialReference()    
    outSpatialRef.ImportFromProj4(proj_trg)
    outSpatialRef.MorphToESRI()
    profile_interval = 10 # make a profile every xx pixels, NOW FIXED, MAKE CONFIGURABLE
    logger.info('Reading ' + rivermap)
    x, y, riversid, FillVal = readMap(rivermap, 'PCRaster');riversid[riversid==FillVal]= -1
    logger.info('Reading ' + drainmap)
    x, y, drain, FillVal = readMap(drainmap, 'PCRaster');drain[drain==FillVal]= np.nan
    logger.info('Reading ' + lddmap)
    x, y, ldd, FillVal = readMap(lddmap, 'PCRaster');
    logger.info('Reading ' + depthmap)
    x, y, depth, FillVal = readMap(depthmap, 'PCRaster');depth[depth==FillVal]= np.nan
    logger.info('Reading ' + widthmap)
    x, y, width, FillVal = readMap(widthmap, 'PCRaster');width[width==FillVal]= np.nan
    logger.info('Reading ' + demmap)
    x, y, dem, FillVal = readMap(demmap, 'PCRaster');dem[dem==FillVal]= np.nan
    xi,yi = np.meshgrid(x,y)
    
    # mesh of surrounding pixels
    xi_window, yi_window = np.meshgrid(range(-1,2),range(-1,2))
    # mesh of ldd grid values
    ldd_values = np.array([[7, 8, 9],[4, 5, 6],[1, 2, 3]])
    
    
    [iiy,iix] = np.where(riversid>0)
    
#    longitudes = xi[iiy,iix]
#    latitudes  = yi[iiy,iix]
    riverId    = riversid[iiy,iix]
#    drainage   = drain[iiy,iix]
#    depths     = depth[iiy,iix]
#    widths     = width[iiy,iix]
    maxRiverId = riverId.max()
    
    # Create new shapefile
    ogr.UseExceptions()
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(SHP_FILENAME)
    layer_line = ds.CreateLayer("rivers", None, ogr.wkbLineString)
    layer_point = ds.CreateLayer("profiles", None, ogr.wkbPoint)
    layer_line_latlon = ds.CreateLayer("rivers_latlon", None, ogr.wkbLineString)
    layer_point_latlon = ds.CreateLayer("profiles_latlon", None, ogr.wkbPoint)
    point_ID = ogr.FieldDefn()
    point_ID.SetName('POINT_ID')
    point_ID.SetType(ogr.OFTInteger)
    point_ID.SetWidth(4)
    layer_point.CreateField(point_ID)
    layer_point_latlon.CreateField(point_ID)
    
    river_ID = ogr.FieldDefn()
    river_ID.SetName('RIVER_ID')
    river_ID.SetType(ogr.OFTInteger)
    river_ID.SetWidth(4)
    layer_point.CreateField(river_ID)
    layer_point_latlon.CreateField(river_ID)
    layer_line.CreateField(river_ID)
    layer_line_latlon.CreateField(river_ID)
    
    point_PROF = ogr.FieldDefn()
    point_PROF.SetName('PROF')
    point_PROF.SetType(ogr.OFTString)
    point_PROF.SetWidth(14)
    layer_point.CreateField(point_PROF)
    layer_point_latlon.CreateField(point_PROF)
    
    point_width = ogr.FieldDefn()
    point_width.SetName('WIDTH')
    point_width.SetType(ogr.OFTReal)
    point_width.SetWidth(6)
    layer_point.CreateField(point_width)
    layer_point_latlon.CreateField(point_width)
    
    point_depth = ogr.FieldDefn()
    point_depth.SetName('DEPTH')
    point_depth.SetType(ogr.OFTReal)
    point_depth.SetWidth(6)
    layer_point.CreateField(point_depth)
    layer_point_latlon.CreateField(point_depth)
    
    point_surf = ogr.FieldDefn()
    point_surf.SetName('SURFACE')
    point_surf.SetType(ogr.OFTReal)
    point_surf.SetWidth(6)
    layer_point.CreateField(point_surf)
    layer_point_latlon.CreateField(point_surf)
    
    # Create a new line geometry per river segment
    for id in np.arange(0,maxRiverId + 1):
        logger.info('Writing line element "' + str(id) + '"')
        # pdb.set_trace()
        y_idx, x_idx = np.where(riversid == id)
        drain_idx = drain[y_idx, x_idx]
        lat_select = yi[y_idx, x_idx]
        lon_select = xi[y_idx, x_idx]
        width_select=width[y_idx, x_idx]
        depth_select=depth[y_idx, x_idx]
        dem_select=dem[y_idx, x_idx]
        # now find ascending drainage order
        order = drain_idx.argsort()
        lat_select = lat_select[order]
        lon_select = lon_select[order]
        width_select=width_select[order]
        depth_select=depth_select[order]
        dem_select=dem_select[order]
        
        # convert lat-lon coordinates to projected coordinates
        x_proj_select, y_proj_select = convertCoord(proj_src, proj_trg, lon_select, lat_select)
        
        # create new line segments in shapefiles
        line_latlon = ogr.Geometry(type=ogr.wkbLineString)
        line = ogr.Geometry(type=ogr.wkbLineString)
        # add points sequentially to line segment
        count = profile_interval # always write a profile in the first
        for nr in range(0,len(lat_select)):
            line_latlon.AddPoint(np.float64(lon_select[nr]), np.float64(lat_select[nr]))
            line.AddPoint(np.float64(x_proj_select[nr]), np.float64(y_proj_select[nr]))
            if nr > 0:
                if count == profile_interval:
                    point_latlon = ogr.Geometry(type=ogr.wkbPoint)
                    point_latlon.SetPoint(0, np.float64(lon_select[nr]), np.float64(lat_select[nr]))
                    feature_point_latlon = ogr.Feature(layer_point_latlon.GetLayerDefn())
                    feature_point_latlon = ogr.Feature(layer_point_latlon.GetLayerDefn())
                    feature_point_latlon.SetGeometryDirectly(point_latlon)
                    # write the attributes to the current feature
                    feature_point_latlon.SetField('POINT_ID', nr)
                    feature_point_latlon.SetField('RIVER_ID', int(id))
                    feature_point_latlon.SetField('WIDTH', 100*np.float64(width_select[nr]))
                    feature_point_latlon.SetField('DEPTH', 100*np.float64(depth_select[nr]))
                    feature_point_latlon.SetField('SURFACE', 100*np.float64(dem_select[nr]))
                    feature_point_latlon.SetField('PROF', str('PROF_%04.f_%04.f') % (np.int8(id), nr))
                    layer_point_latlon.CreateFeature(feature_point_latlon)
                    feature_point_latlon.Destroy()

                    point = ogr.Geometry(type=ogr.wkbPoint)
                    point.SetPoint(0, np.float64(x_proj_select[nr]), np.float64(y_proj_select[nr]))
                    feature_point = ogr.Feature(layer_point.GetLayerDefn())
                    feature_point.SetGeometryDirectly(point)
                    # write the attributes to the current feature
                    feature_point.SetField('POINT_ID', nr)
                    feature_point.SetField('RIVER_ID', int(id))
                    feature_point.SetField('WIDTH', 100*np.float64(width_select[nr]))
                    feature_point.SetField('DEPTH', 100*np.float64(depth_select[nr]))
                    feature_point.SetField('SURFACE', 100*np.float64(dem_select[nr]))
                    feature_point.SetField('PROF', str('PROF_%04.f_%04.f') % (np.int8(id), nr))
                    layer_point.CreateFeature(feature_point)
                    feature_point.Destroy()
                    count = 0
                    if nr == 1:
                        print 'Writing profile in %gst point in river section "%g"' % (nr, id)
                    elif nr == 2:
                        print 'Writing profile in %gnd point in river section "%g"' % (nr, id)
                    elif nr == 3:
                        print 'Writing profile in %grd point in river section "%g"' % (nr, id)
                    else:
                        print 'Writing profile in %gth point in river section "%g"' % (nr, id)
                count += 1
                
        # now find the point downstream of the last pixel from the ldd, which is connected with the downstream river
        try:
            xi_select       = xi[y_idx[order][-1] + yi_window, x_idx[order][-1] + xi_window]
            yi_select       = yi[y_idx[order][-1] + yi_window, x_idx[order][-1] + xi_window]
            # convert to projected coordinates
            xi_proj_select, yi_proj_select = convertCoord(proj_src, proj_trg, x, y)
            
            ldd_at_pos      = ldd[y_idx[order][-1], x_idx[order][-1]]
            ldd_y, ldd_x    = np.where(ldd_values==ldd_at_pos)
            downstream_y    = yi_select[ldd_y, ldd_x]
            downstream_x    = xi_select[ldd_y, ldd_x]
            downstream_y_proj = yi_proj_select[ldd_y, ldd_x]
            downstream_x_proj = xi_proj_select[ldd_y, ldd_x]
            
            line_latlon.AddPoint(np.float64(downstream_x), np.float64(downstream_y))
            line.AddPoint(np.float64(downstream_x_proj), np.float64(downstream_y_proj))
        except:
            # most downstream point of segment is on the boundary of the map, so skip this step
            print 'River segment id: %g is on boundary of the map' % id
        # Add line as a new feature to the shapefiles
        feature_latlon = ogr.Feature(feature_def=layer_line_latlon.GetLayerDefn())
        feature_latlon.SetGeometryDirectly(line_latlon)
        feature_latlon.SetField('RIVER_ID', int(id))
        layer_line_latlon.CreateFeature(feature_latlon)
        # Cleanup
        feature_latlon.Destroy()

        feature = ogr.Feature(feature_def=layer_line.GetLayerDefn())
        feature.SetGeometryDirectly(line)
        feature.SetField('RIVER_ID', int(id))
        layer_line.CreateFeature(feature)
        # Cleanup
        feature.Destroy()




    ds.Destroy()
    # finally write .prj files
    rivers_prj = os.path.join(SHP_FILENAME, 'rivers.prj')
    profiles_prj = os.path.join(SHP_FILENAME, 'profiles.prj')
    rivers_prj_latlon = os.path.join(SHP_FILENAME, 'rivers_latlon.prj')
    profiles_prj_latlon = os.path.join(SHP_FILENAME, 'profiles_latlon.prj')
    rivers_prj_file = open(rivers_prj, 'w');rivers_prj_file.write(outSpatialRef.ExportToWkt());rivers_prj_file.close()
    profiles_prj_file = open(profiles_prj, 'w');profiles_prj_file.write(outSpatialRef.ExportToWkt());profiles_prj_file.close()
    rivers_prj_latlon_file = open(rivers_prj_latlon, 'w');rivers_prj_latlon_file.write(outSpatialRef_latlon.ExportToWkt());rivers_prj_latlon_file.close()
    profiles_prj_latlon_file = open(profiles_prj_latlon, 'w');profiles_prj_latlon_file.write(outSpatialRef_latlon.ExportToWkt());profiles_prj_latlon_file.close()

#catchpointmap   = r'd:\osm2hydro_cases\limpopo\limpopo\dem\catch_point.map'
#catchsurfacemap = r'd:\osm2hydro_cases\limpopo\limpopo\dem\catchments.map'
#proj_trg        = '+proj=utm +zone=36 +south +ellps=WGS84 +units=m +no_defs '
#proj_src        = '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'
#SHP_FILENAME    = r'd:\osm2hydro_cases\limpopo\limpopo\modelshapes\SOBEK_shapes'
#PCR_catch2shape(catchpointmap, catchsurfacemap, proj_src, proj_trg, SHP_FILENAME)