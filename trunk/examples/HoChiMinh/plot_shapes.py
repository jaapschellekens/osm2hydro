# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:39:26 2013

@author: Hessel Winsemius
"""
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import LineCollection
import matplotlib.font_manager as font
from matplotlib import cm
import numpy as np
import osgeo.ogr as ogr
import pdb


def load_shape_point(shapePath):
    'Given a shapePath, return a list of points in GIS coordinates'
    # Open shapeData
    shapeData = ogr.Open(shapePath)
    # Get the first layer
    layer = shapeData.GetLayer()
    # Initialize
    X = []
    Y = []
    # For each point,
    for index in xrange(layer.GetFeatureCount()):
        # Get
        feature = layer.GetFeature(index)
        geometry = feature.GetGeometryRef()
        # Get pointCoordinates
        pointCoordinates = geometry.GetX(), geometry.GetY()
        # Append
        X.append(pointCoordinates[0])
        Y.append(pointCoordinates[1])
        # Cleanup
        feature.Destroy()
    # Get spatial reference as proj4
    proj4 = layer.GetSpatialRef().ExportToProj4()
    # Cleanup
    shapeData.Destroy()
    # Return
    return X, Y, proj4


def setup_Basemap(xmin, xmax, ymin, ymax, projection='cyl', resolution='c', \
    drawcoastlines=True, drawrivers=True,drawcountries=True, \
    meridiantickdist=0, paralleltickdist=0):
    m = Basemap(projection='cyl', llcrnrlon=xmin, \
            urcrnrlon=xmax, llcrnrlat=ymin, urcrnrlat=ymax, resolution=resolution)
    #xi,yi = m(*np.meshgrid(x, y))
    if drawcoastlines:
        m.drawcoastlines(linewidth=0.25)
    if drawrivers:
        m.drawrivers(color='b', linewidth=0.25)
    if drawcountries:
        m.drawcountries(linewidth=0.25)
    
    m.drawmapboundary(linewidth=0.5)
    if meridiantickdist > 0:
        m.drawmeridians(np.arange(xmin,xmax+1,meridiantickdist),labels=[1,0,0,1],linewidth=0.5, fontsize=10)
    if paralleltickdist > 0:
        m.drawparallels(np.arange(ymin,ymax+1,paralleltickdist),labels=[1,0,0,1],linewidth=0.5, fontsize=10)
    return m
#im = m.contourf(xi, yi, data, norm, cmap=cmap, vmin=vmin, vmax=vmax, extend='max')

xmin = 98
xmax = 102
ymin = 12
ymax = 18
OSMshapefile = r'd:\osm2hydro_cases\Thailand\osmshapes\lu_waterway_all'
OSMfiltered  = r'd:\osm2hydro_cases\Thailand\modelshapes\streams'
SOBEKnetwork = r'd:\osm2hydro_cases\Thailand\modelshapes\test'
laterals     = r'd:\osm2hydro_cases\Thailand\modelshapes\SOBEK_shapes\laterals_latlon.shp'

figTitles  = ['a. OSM all waterways' , 'b. OSM main waterways', 'c. Hydraulic schematisation']

plt.close('all')
fig = plt.figure(figsize=(12,7))
fig.subplots_adjust(left=0.05, right=0.95, wspace=0.25)
n = 0
plt.subplot(1,3,n+1)
m = setup_Basemap(xmin, xmax, ymin, ymax, resolution='h', drawrivers=False, meridiantickdist=1., paralleltickdist=1.)
m.readshapefile(OSMshapefile, 'waterway', color='b', linewidth=.2)

plt.title(figTitles[n], fontsize=12)

n += 1
plt.subplot(1,3,n+1)
m = setup_Basemap(xmin, xmax, ymin, ymax, resolution='h', drawrivers=False, meridiantickdist=1., paralleltickdist=1.)
m.readshapefile(OSMfiltered, 'waterway', color='b', linewidth=.2)
plt.title(figTitles[n], fontsize=12)

n += 1
plt.subplot(1,3,n+1)
m = setup_Basemap(xmin, xmax, ymin, ymax, resolution='h', drawrivers=False, meridiantickdist=1., paralleltickdist=1.)
p1 = m.readshapefile(SOBEKnetwork, 'RIVER_ID', color='r', linewidth=1.)

lons, lats, proj4 = load_shape_point(laterals)
x, y = m(lons,lats)
m.scatter(x,y,30,marker='o',color='k')

#p2 = m.readshapefile(laterals, 'POINT_ID', color='g')
plt.title(figTitles[n], fontsize=12)
plt.savefig('Thailand.pdf')
plt.savefig('Thailand.png')

plt.show()
