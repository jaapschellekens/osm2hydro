# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:03:41 2013

@author: schelle
"""

import numpy

try:
    from  wflow.wf_DynamicFramework import *
except ImportError:
    from  wf_DynamicFramework import *
try:
    from  wflow.wflow_funcs import *
except ImportError:
    from  wflow_funcs import *
    
    

corina = readmap("pavedcorinaresamp.map")
roads_den=readmap("roads_den.map")
lu_water=readmap("lu_water.map")
lu_unpaved=readmap("lu_unpaved.map")
lu_roads=readmap("lu_roads.map")
lu_paved=readmap("lu_paved.map")
lu_paved_fromroads=readmap("lu_paved_fromroads.map")
lu_pavedpol=readmap("lu_pavedpol.map")
subcatch=cover(readmap("wflow_subcatch.map"),0)

idar = pcr2numpy(areaaverage(scalar(subcatch),subcatch),0)

rd_=unique((pcr2numpy(areaaverage(lu_roads,subcatch),0) * 1000 + 10 * (idar + 1)).astype(int))
cor_=unique((pcr2numpy(areaaverage(corina,subcatch),0) * 1000 + 10 * (idar + 1)).astype(int))
osm_=unique((pcr2numpy(areaaverage(lu_paved,subcatch),0) * 1000 + 10 * (idar + 1)).astype(int))
pol_=unique((pcr2numpy(areaaverage(lu_pavedpol,subcatch),0) * 1000 + 10 * (idar + 1)).astype(int))

difar = unique((idar + 1) * 10)
print difar
print cor_ - difar 
print osm_ - difar
print rd_ - difar
print pol_ - difar

report(areaaverage(corina,subcatch),"corina_avg.map")
report(areaaverage(lu_paved,subcatch),"lu_paved_avg.map")


