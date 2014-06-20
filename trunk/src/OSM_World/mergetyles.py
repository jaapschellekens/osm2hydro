"""
Simple script to merge the tiles into single maps using gdal_merge
"""

import glob
import fnmatch
import os

regions = ['world_11.poly']

matches = []
patterns = ["roads_den.tif","lu_paved.tif","lu_pavedpol.tif","lu_paved_fromroads.tif","lu_roads.tif"]

def getflist(subdir,pattern):
    for dirName, subdirList, fileList in os.walk(subdir):
        for filename in fnmatch.filter(fileList, pattern):
            matches.append(os.path.join(dirName,filename))
    return matches


for pattern in patterns:
    for region in regions:
        matches = getflist(region,pattern)

        a = "\n"
        inflist = a.join(matches)

        fp = open('flist.txt','w')
        fp.write(inflist)
        fp.close()
        ofile = region + pattern
        mergestr = "gdal_merge.py -o " + ofile + " --optfile flist.txt"
        print mergestr
        os.system(mergestr)

