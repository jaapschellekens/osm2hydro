# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:40:09 2013

@author: schelle
"""

from gdal_readmap import gdal_readmap
from basemap_imshow import basemap_imshow
import numpy as np
import matplotlib as mpl
from matplotlib.pyplot import *



xx,yy,osm,zz = gdal_readmap('lu_paved_lux.tif','GTiff')
x,y,corina,zz = gdal_readmap('pavedcorina_lux.map','PCRaster')
corina = np.round(corina)

cmap = cm.jet
carray = np.array([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

filt = osm < carray.min()
osm[filt] = np.nan
filt = corina < carray.min()
corina[filt] = np.nan
Norm = mpl.colors.BoundaryNorm(carray,cmap.N)

fig = figure()
subplot(2,1,1)

minx = 6
miny = 49.55
maxx = 6.25
maxy = 49.65

basemap_imshow(x,np.flipud(y),np.flipud(corina),minx,maxx,miny,maxy,shapefile='lu_roads_main',norm=Norm,meridian_dist=0.05,parallel_dist=0.05,resolution='h')
title('Corine',fontsize=11)
subplot(2,1,2)
m,im=basemap_imshow(xx,np.flipud(yy),np.flipud(osm),minx,maxx,miny,maxy,shapefile='lu_roads_main',norm=Norm,meridian_dist=0.05,parallel_dist=0.05,resolution='h')
m.drawmapscale(6.2, 49.57, 6.2, 49.57, 10, barstyle='fancy')
title('OSM',fontsize=11)


axx = fig.add_axes([0.17,0.05,0.68,0.02])
fig.colorbar(im,cax=axx,orientation='horizontal', extend='min')
savefig('lux-two.pdf',dpi=300)
savefig('lux-two.png',dpi=300)
