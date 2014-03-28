# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:22:36 2012

 Copyright notice

    Copyright 2013, 2014 Hesssel Winsemius

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

# test for wavelet filtering of SRTM elevation models

#create some data
import pywt
import numpy as np
import pdb
#file = 'Bangkok_cut.tif'
#filter_name='bior6.8'
# find misses and fill with a zero

# max. level of waveforms
#levels = int(np.floor(np.log2(zi2.shape[0])))

def dem_filter(dem, filter_name, filter_weight):
    xBlockSize = 1000
    yBlockSize = 1000
    # load full DEM in memory
    # x, y, dem, FillVal = readMap(file, 'GTiff')
    ii = np.isnan(dem)
    dem[ii] = 0.
    # if dem bigger than a certain area, chop the procedure in small pieces
    dem_new = np.zeros(dem.shape)
    rows = dem.shape[0]
    cols = dem.shape[1]
    
    for i in range(0,rows, yBlockSize):
        if i + yBlockSize < rows:
            numRows = yBlockSize
        else:
            numRows = rows - i
            if numRows % 2 == 1:
                # round to a even number
                numRows -= 1
        i2 = i + numRows
        for j in range(0, cols, xBlockSize):
            if j + xBlockSize < cols:
                numCols = xBlockSize
            else:
                numCols = cols - j
                if numCols % 2 == 1:
                    # round to a even number
                    numCols -= 1
            j2 = j + numCols
            print 'Filtering data-block y: %g -- %g; x: %g -- %g' % (i, i2, j, j2)
            coeffs = pywt.wavedec2(dem[i:i2, j:j2], filter_name, level=5) # level taken from paper
            filter_weight = np.float64(filter_weight)
            # now remove some thresholds
             # noiseSigma*sqrt(2*log2(image.size))
            NewWaveletCoeffs = map (lambda x: pywt.thresholding.soft(x,filter_weight),coeffs)
            dem_new[i:i2, j:j2] = pywt.waverec2(NewWaveletCoeffs, filter_name)
    
    # now bring back missing values
    dem_new[ii] = np.nan
    return dem_new
    # dem[ii] = np.nan
    # writeMap('Bangkok_cut_filter.tif', 'GTiff', x, y, dem_filter, FillVal)


# now reconstruct the original:
# zi_filter = pywt.waverec2(coeffs, 'bior1.1')

#pyplot.figure()
#pyplot.imshow(np.minimum(dem,10),interpolation='nearest')
#
#pyplot.figure()
#pyplot.imshow(np.minimum(dem_filter,10),interpolation='nearest')
#
#pyplot.figure()
#pyplot.plot(dem[650, 1000:1200],label='unfiltered')
#pyplot.plot(dem_filter[650, 1000:1200],label=str('filter %s level=%2.f thres=%2.2f') % (filter_name, level, threshold))
#pyplot.legend()
#from numpy import *
#import scipy
#import pywt
#
#image = scipy.misc.lena().astype(float32)
#
#noiseSigma = 16.0
#image += random.normal(0, noiseSigma, size=image.shape)
#
#wavelet = pywt.Wavelet('haar')
#levels  = int(floor( log2(image.shape[0]) ))
#
#
#WaveletCoeffs = pywt.wavedec2( image, wavelet, level=levels)