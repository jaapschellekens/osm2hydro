A set of hydrological tools for OpenStreetMap data.


What is included?
=================

dist/.................... binary version of all tools (+ generated from Python using Freezer)
  win32/ 
  linux/

examples/
  Arhnem/................ test case

src/
  osm2hydro/............. Python tools
    gdal_density.py...... ???
    gdal_merge.py........ part of GDAL + 1 bug fixed (better contribute it back and remove)
    map2shape.py......... ???
    osm2hydro.py......... ???
    osm2shp.py........... 
    poly_density.py...... convert shape to fraction coverage per pixel
  osm2shp/............... modified version of osm2shp, extended with filtering
  osmconvert/............ part of OSM, extracts subset of OSM file
  osmfilter/............. part of OSM
  setup.py

third-party/............. external tools, libraries (included, but not modified)


Installation
============

>> TODO <<

Do NOT install osm2hydro in a directory with a space in the name (e.g. "Program Files")

Required libraries

- GDAL >- 10.1.0 (Both the executables and the python bindings)
- PyWavelets - http://www.pybytes.com/pywavelets/
- PyProj - http://code.google.com/p/pyproj/


Usage
=====

Most of the functionality of the osm2hydro tools can be found in the examples/.



Developing on Windows (only needed for active development)
==========================================================

Download and Instally the following Python packages and development tools:

* PythonXY 2.7.x: https://code.google.com/p/pythonxy/wiki/Downloads
* Pyproj: http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyproj
* Python Tools for Visual Studio 2.0: https://pytools.codeplex.com/releases/view/103102
* GDAL 1.10: http://www.gisinternals.com/sdk/PackageList.aspx?file=release-1600-gdal-1-10-mapserver-6-4.zip
* GDAL 1.10 Python: http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
* basemap for matplotlib: http://sourceforge.net/projects/matplotlib/files/matplotlib-toolkits/

* Shapely 1.2.18: http://www.lfd.uci.edu/~gohlke/pythonlibs/
* Kartograph: http://kartograph.org/docs/kartograph.py/
* tinycss: pip install tinycss

If you plan to view examples/ using Python notebook (run in python):

> from IPython.external import mathjax; mathjax.install_mathjax()

Run the following in console to link osm2hydro sources to your python installation

> cd src
> python setup.py develop

At this point osm2hydro can be used as a command-line tool or as a library:

> import osm2hydro


Enjoy!


Developing on UNIX
==================

>> TODO: Jaap <<


Using IPython Notebook
======================

Make sure that osm2hydro is installed or linked (using "python setup.py develop")

> cd examples/Arnhem
> ipython notebook --pylab=inline

Then you should be able to browse and use notebooks in the current directory (files with .ipynb extension).

Code Contributors
=================

Hessel.Winsemius@deltares.nl
Jaap.Schellekens@deltares.nl
Gennadii.Donchyts@deltares.nl


Citations
=========

>> TODO: update final name once the article is submitted <<

Building rapid assessment hydrological models using open data and OpenStreetMap
J. Schellekens, R.J. Brolsma, R.J. Dahm, G.V. Donchyts, H.C. Winsemius


