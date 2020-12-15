"""
main script for compiling the GMC and cluster catalog info and creating all the ncessary ds9 region files

2020-10-27
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
import sys

sys.path.append('/cherokee1/turner/phangs/cf/utils')
from gmc_cat_utils import *
from cluster_cat_utils import *
from image_utils import *
from plot_utils import *

# set path of the data dir
data_dir = '/cherokee1/turner/phangs/cf/data/'

# read in list of galaxies 
master_galaxy_list = ascii.read('master_galaxy.list')
galaxy_list = ascii.read('galaxy.list')

# create ds9 region files for all the galaxies
# all_galaxies_gmc_region_files(master_galaxy_list, data_dir)

# # print out number of clusters and make the new fits files if needed
# all_galaxies_clust_cats(master_galaxy_list, data_dir, mkfits=False)

# # create ds9 region files for all the galaxies
# all_galaxies_clust_region_files(master_galaxy_list, data_dir, radius_pix=10)

# # create ds9 region files for all the HST footprints
# all_galaxies_flc_region_files(galaxy_list, data_dir)


# play around with radii for plotting the outline plots
#   	  ngc0628 ngc1433 ngc1559 ngc1566 ngc1792 ngc3351 ngc3627 ngc4535
radius = [0.04, 0.03, 0.024, 0.03, 0.025, 0.025, 0.04, 0.025]

# all_galaxies_outline_plots(galaxy_list, data_dir, sc_class='class12', radius=radius)

# now with the 3-color trilogy images as the background
# all_galaxies_outline_plots(galaxy_list, data_dir, sc_class='class12', radius=radius, bkgd='trilogy')

# outline plots color-coded by cluster age
# radius = [0.025]
all_galaxies_outline_plots(galaxy_list, data_dir, sc_class='class12', radius=radius, color_code='age')
