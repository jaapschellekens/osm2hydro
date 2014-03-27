# used to develop and package (egg), list of most commands: http://pythonhosted.org/setuptools/setuptools.html, http://docs.python.org/2/distutils/setupscript.html, http://peak.telecommunity.com/DevCenter/PythonEggs
#
# python setup.py freeze
# python setup.py develop
# python setup.py bdist_wininst
# python setup.py sdist
# python setup.py bdist_egg


from setuptools import setup, find_packages
from distutils.cmd import Command
from distutils.core import setup
import os

version = '1.0'

class freeze(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        from bbfreeze import Freezer
        f = Freezer("dist/win32/osm2hydro") 
        f.addScript("../src/osm2hydro/dem_filter.py")
        f.addScript("../src/osm2hydro/gdal_density.py")
        f.addScript("../src/osm2hydro/gdal_merge.py")
        f.addScript("../src/osm2hydro/osm2hydro.py")
        f.addScript("../src/osm2hydro/poly_density.py")
        f()    # starts the freezing process
      

setup(name='osm2hydro',
      version=version,
      description="A set of tools for hydrological processing using OSM and other open data sources",
      long_description="""\
osm2hydro consist of a set of tools used to generate input for water-related models (models which simulate hydrographic features such as rivers, lakes, watersheeds, etc.)""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='hydro osm tool gdal',
      author='Hessel Winsemius',
      author_email='hessel.winsemius@deltares.nl',
      url='http://osm2hydro.googlecode.com',
      license='LGPL',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      cmdclass = {'freeze': freeze},
      install_requires=[
          # -*- Extra requirements: -*-
      ]
      )


