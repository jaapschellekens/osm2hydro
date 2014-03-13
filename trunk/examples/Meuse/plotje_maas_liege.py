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

x,y,pavedpol,zz = gdal_readmap('lu_pavedpol.map','PCRaster')
x,y,pavedrd,zz = gdal_readmap('lu_paved_fromroads.map','PCRaster')
x,y,paved,zz = gdal_readmap('lu_paved.map','PCRaster')
x,y,corina,zz = gdal_readmap('pavedcorinaresamp.map','PCRaster')


cmap = cm.jet
carray = np.array([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
filt = pavedpol < carray.min()
pavedpol[filt] = np.nan
filt = pavedrd < carray.min()
pavedrd[filt] = np.nan
filt = paved < carray.min()
paved[filt] = np.nan
filt = corina < carray.min()
corina[filt] = np.nan
Norm = mpl.colors.BoundaryNorm(carray,cmap.N)

fig = figure(figsize=(9,8))
subplot(2,2,1)

minx = 5.1
miny = 50.3
maxx = 5.9
maxy = 50.9

basemap_imshow(x,np.flipud(y),np.flipud(corina),minx,maxx,miny,maxy,norm=Norm,meridian_dist=0.2,parallel_dist=0.2,resolution='h')
title('Corine',fontsize=11)
subplot(2,2,2)
basemap_imshow(x,np.flipud(y),np.flipud(pavedpol),minx,maxx,miny,maxy,norm=Norm,meridian_dist=0.2,parallel_dist=0.2,resolution='h')
title('OSM Polygons',fontsize=11)
subplot(2,2,3)
basemap_imshow(x,np.flipud(y),np.flipud(pavedrd),minx,maxx,miny,maxy,norm=Norm,meridian_dist=0.2,parallel_dist=0.2,resolution='h')
title('OSM Roads',fontsize=11)
subplot(2,2,4)
m, im = basemap_imshow(x,np.flipud(y),np.flipud(paved),minx,maxx,miny,maxy,norm=Norm,meridian_dist=0.2,parallel_dist=0.2,resolution='h')
m.drawmapscale(5.3, 50.38, 5.3, 50.38, 20, barstyle='fancy')
title('OSM combined',fontsize=11)
axx = fig.add_axes([0.15,0.03,0.7,0.02])
fig.colorbar(im,cax=axx,orientation='horizontal',extend='min')
savefig('liege-four.pdf',dpi=300)
savefig('liege-four.png',dpi=300)