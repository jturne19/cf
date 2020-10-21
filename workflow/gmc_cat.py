"""
main script for compiling the cluster catalog info and making class 1 and 2 catalogs and whatnot

2020-10-13
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
import sys

sys.path.append('/cherokee1/turner/phangs/cf/utils')

from gmc_cat_utils import *

# set path of the data dir
data_dir = '/cherokee1/turner/phangs/'

# read in list of galaxies 
galaxy_list = ascii.read('master_galaxy.list')

all_galaxies_gmc_region_files(galaxy_list, data_dir)