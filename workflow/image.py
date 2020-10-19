"""


"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.wcs import WCS
import glob
import sys

sys.path.append('/cherokee1/turner/phangs/cf/utils')

from image_utils import *

# set path of the data dir
data_dir = '/cherokee1/turner/phangs/'

# read in master list of galaxies 
galaxy_list = ascii.read('galaxy.list')


all_galaxies_flc_region_files(galaxy_list, data_dir)