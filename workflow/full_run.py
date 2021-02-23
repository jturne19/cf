"""
script to hopefully go through a full run with proper directory naming and what not

2021-02-22
"""
import sys

import numpy as np 
from astropy.io import ascii
from astropy.table import Table
import matplotlib

sys.path.append('/cherokee1/turner/phangs/cf/utils')
from utils import *

# path to the directory with all the galaxies and data
data_dir = '/cherokee1/turner/phangs/cf/data/'

# parameters to change for each run

# file with the list of galaxies
galaxy_list_file = 'galaxy.list'

# name of the run
run_name = 'class12'

# do you want matplotlib interactive plots? (probably not)
interactive_plots = False

# which GMC cat are you using?
gmc_cat_suffix = '_12m+7m+tp_co21_nativeres_nativenoise_props'

# which base cluster catalog are you using?
sc_base_cat_suffix = '_phangshst_base_catalog'

# for analyses, class 1,2 or class 1,2,3?
sc_class = 'class12'

# perform catalog_setup steps? (making ds9 region files, outline plots)
do_catalog_setup = False

# perform the nearest neighbor analysis?
do_nn = False

# perform the star cluster - gmc association analysis?
do_assoc = False

# perform the two point corelation function analysis?
do_tpcf = False

# perform the cross-correlation analysis?
do_cf = False

# radius of the field of view used in the overview plots
#   	  ngc0628 ngc1365 ngc1433 ngc1559 ngc1566 ngc1792 ngc3351 ngc3627 ngc4535 ngc4548 ngc4571
radius = [0.04,   0.03,   0.03,   0.024,  0.03,   0.025,  0.025,  0.04,   0.025,  0.025,  0.025]
radius_no_3627 = [0.04, 0.03, 0.03, 0.024, 0.03, 0.025, 0.025, 0.025, 0.025, 0.025]

# set number of bins for the two-point correlation functions
nbins_tpcf = 10

# set number of bins for the cross correlation functions
nbins_crosscf = 10

# -- do not change parameters below this ---
for arg in (sys.argv)[1:]:
	exec(arg)

# read in the list of galaxies
galaxy_list = ascii.read(galaxy_list_file)

galaxy_list_no_3627 = Table(galaxy_list)
galaxy_list_no_3627.add_index('id')

try:
	row_3627 = galaxy_list_no_3627.loc_indices['ngc3627']
	galaxy_list_no_3627.remove_row(row_3627)
except:
	print('ngc3627 is not in the galaxy list so no need to worry about the trilogy image for this run')

if interactive_plots:
	matplotlib.use('Qt5agg')
else:
	matplotlib.use('agg')

print('')
print('running %s for %s clusters'%(run_name, sc_class))
print('')

if do_catalog_setup:

	# create ds9 region files for all the galaxies
	print('\nmaking gmc region files')
	print('-----------------------------------------------')
	all_galaxies_gmc_region_files(galaxy_list, data_dir, run_name, gmc_cat_suffix=gmc_cat_suffix)
	
	# print out number of star clusters and make the new fits files if needed
	print('\nmaking cluster catalogs for %s'%sc_class)
	print('-----------------------------------------------')
	all_galaxies_clust_cats(galaxy_list, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix, sc_class=sc_class, mkfits=True)
	
	# create ds9 region files for all the galaxies
	print('\nmaking cluster region files')
	print('-----------------------------------------------')
	all_galaxies_clust_region_files(galaxy_list, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix, sc_class=sc_class, radius_pix=10)
	
	# create gmc catalog for just the hst-alma overlap region
	print('\nmaking gmc catalog for the hst-alma overlap mask')
	print('-----------------------------------------------')
	generate_gmc_cat_masked(galaxy_list, data_dir, run_name, gmc_cat_suffix=gmc_cat_suffix)

	# create the basic outline plots 
	print('\nmaking basic outline plots for %s'%sc_class)
	print('-----------------------------------------------')
	all_galaxies_outline_plots(galaxy_list, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix, 
							   sc_class=sc_class, radius=radius)
	
	# now with the 3-color trilogy images as the background (ngc3627's trilogy image is too big for this)
	print('\nmaking outline plots with trilogy background for %s'%sc_class)
	print('-----------------------------------------------')
	all_galaxies_outline_plots(galaxy_list_no_3627, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix, 
							   sc_class=sc_class, radius=radius_no_3627, bkgd='trilogy')
	
	# outline plots color-coded by cluster age
	print('\nmaking outline plots colored by cluster age for %s'%sc_class)
	print('-----------------------------------------------')
	all_galaxies_outline_plots(galaxy_list, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix,
							  sc_class=sc_class, radius=radius, color_code='age')
	
	# outline plots color-coded by cluster mass
	print('\nmaking outline plots colored by cluster mass for %s'%sc_class)
	print('-----------------------------------------------')
	all_galaxies_outline_plots(galaxy_list, data_dir, run_name, sc_base_cat_suffix=sc_base_cat_suffix, gmc_cat_suffix=gmc_cat_suffix,
							  sc_class=sc_class, radius=radius, color_code='mass')
	
if do_nn:

	# nearest-neighbor analysis
	print('\nnearest neighbor analysis')
	print('-----------------------------------------------')
	all_galaxies_sc_gmc_sep(galaxy_list, data_dir, run_name, mkhist=True, sc_cat_suffix='%s.%s'%(sc_base_cat_suffix,sc_class), 
							gmc_cat_suffix=gmc_cat_suffix, output_cat_suffix='_cluster_catalog_in_mask_%s'%sc_class)

if do_assoc:

	# star cluster - gmc association analysis
	print('\nstar cluster - gmc association analysis')
	print('-----------------------------------------------')
	all_galaxies_sc_gmc_assoc(galaxy_list, data_dir, run_name, sc_mask_cat_suffix='_cluster_catalog_in_mask_%s'%sc_class, 
							  gmc_cat_suffix=gmc_cat_suffix, sc_class=sc_class)
if do_tpcf:

	# two point correlation function analysis
	print('\ntwo point correlation function analysis')
	print('-----------------------------------------------')
	all_galaxies_tpcf(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_%s_assoc_gmc'%sc_class, sc_class=sc_class, 
					  nbins=nbins_tpcf )

if do_cf:

	# two point correlation function analysis
	print('\ncross correlation function analysis')
	print('-----------------------------------------------')
	all_galaxies_crosscf(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_%s_assoc_gmc'%sc_class, sc_class=sc_class, 
						 nbins=nbins_crosscf )
