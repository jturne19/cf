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

# read in list of galaxies 
galaxy_list = ascii.read('galaxy.list')

gal_id 		= galaxy_list['id']
gal_alt_id 	= galaxy_list['alt_id']
gal_dist 	= galaxy_list['dist']


# loop through all galaxies
for i in range(len(galaxy_list)):

	gal_name = gal_id[i]
	gal_alt_name = gal_alt_id[i]

	print(gal_name)

	# check if the base catalog fits files exists
	try:
		base_cat = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits'%(gal_name,gal_name))[1].data

	except FileNotFoundError:
		print(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits not found, skipping'%(gal_name,gal_name))
		continue
	
	n_total    = len(base_cat)
	print('')
	print('%s total number of clusters = %i'%(gal_name, n_total))
	print('')


	# skip ngc4548, ngc4569, ngc4571, ngc4689, ngc5248
	if i < 9 :
	# make new fits files with just class 1,2,3 and 1,2
		new_clust_fits_from_classes(base_cat, [1,2,3], data_dir + 'hst/%s/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))
		new_clust_fits_from_classes(base_cat, [1,2], data_dir + 'hst/%s/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))

		# read those in
		cat123 = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))[1].data
		cat12  = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))[1].data
		
		n_class123 = len(cat123)
		n_class12  = len(cat12)
	
		print('%s number of class 1,2,3 clusters = %i'%(gal_name, n_class123))
		print('%s total number of clusters = %i'%(gal_name, n_class12))
		print('')
