"""
main script for compiling the cluster catalog info and making class 1 and 2 catalogs and whatnot

2020-10-05
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
import sys

sys.path.append('/cherokee1/turner/phangs/cf/utils')

from cluster_cat_utils import *


# set path of the data dir
data_dir = '/cherokee1/turner/phangs/'

# read in master list of galaxies 
galaxy_list = ascii.read('master_galaxy.list')

# print out number of clusters and make the new fits files if needed
# all_galaxies_clust_cats(galaxy_list, data_dir, mkfits=False)

# all_galaxies_clust_region_files(galaxy_list, data_dir, radius_pix=10)
