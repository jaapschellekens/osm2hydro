# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 11:41:58 2013

@author: Hessel Winsemius
"""

import numpy as np
import os, sys
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.font_manager as font
from matplotlib import cm
### PARAMETERS FOR MATPLOTLIB :
import matplotlib as mpl


def basemap_imshow(x, y, data, xmin, xmax, ymin, ymax, projection='merc', \
                resolution='c', vmin=None, vmax=None, norm=None, cmap=cm.jet, \
                drawcoastlines=True, drawrivers=True, drawcountries=True, \
                linewidth=0.5, meridian_dist=30, parallel_dist=15, figTitle='', shapefile=None):
                    
                    
    m           = Basemap(projection=projection, llcrnrlon=xmin, \
            urcrnrlon=xmax, llcrnrlat=ymin, urcrnrlat=ymax, resolution=resolution)
    xi, yi      = np.meshgrid(x, y)
    data_trans  = m.transform_scalar(data,xi[0,:],yi[:,0],xi.shape[1],yi.shape[0], order=0)
    im          = m.imshow(data_trans,cmap,vmin=vmin,vmax=vmax, norm=norm, interpolation='nearest')
    if drawcoastlines:
        m.drawcoastlines(linewidth=linewidth)
    if drawrivers:
        m.drawrivers(color='b', linewidth=linewidth)
    if drawcountries:
        m.drawcountries(linewidth=linewidth)
    m.drawmapboundary(linewidth=linewidth)
    m.drawmeridians(np.arange(xmin,xmax+1,meridian_dist),labels=[1,0,0,1],linewidth=0.5, fontsize=10)
    m.drawparallels(np.arange(ymin,ymax+1,parallel_dist),labels=[1,0,0,1],linewidth=0.5, fontsize=10)
    if shapefile:
        m.readshapefile(shapefile, 'state',color='r')

    plt.title(figTitle, fontsize=10)
    return m, im

