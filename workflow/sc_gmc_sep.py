"""
main script for running the 

2020-12-16
"""
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import aplpy as ap
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u


import sys
sys.path.append('/cherokee1/turner/phangs/cf/utils')
from utils import *


# set path of the data dir
data_dir = '/cherokee1/turner/phangs/cf/data/'

# read in list of galaxies 
master_galaxy_list = ascii.read('master_galaxy.list')
galaxy_list = ascii.read('galaxy.list')

gal_id 		= galaxy_list['id']
gal_alt_id 	= galaxy_list['alt_id']
gal_dist 	= galaxy_list['dist']

# loop through all the galaxies in the list
for i in range(len(galaxy_list)):

	# galaxy props
	gal_name = gal_id[i]
	dist = gal_dist[i]
	print('')
	print(gal_name)

	# read in star cluster catalog [class 1 and 2 only for now]
	sc_cat = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class12.fits'%(gal_name, gal_name))[1].data

	# grab star cluster positions
	sc_x, sc_y,   = sc_cat['PHANGS_X'], sc_cat['PHANGS_Y']
	sc_ra, sc_dec = sc_cat['PHANGS_RA'], sc_cat['PHANGS_DEC']
	# grab star cluster ages
	sc_age = sc_cat['PHANGS_AGE_MINCHISQ']

	# read in GMC catalog 
	gmc_cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits'%(gal_name, gal_name))[1].data

	# grab center positions of the GMCs 
	gmc_ra, gmc_dec = gmc_cat['XCTR_DEG'], gmc_cat['YCTR_DEG']

	# read in the overlap mask
	mask_hdu = fits.open(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name))
	mask = mask_hdu[0].data
	mask_header = mask_hdu[0].header

	# convert star cluster x,y postions to integer pixels
	sc_x_int = np.array([int(np.round(x)) for x in sc_x ])
	sc_y_int = np.array([int(np.round(y)) for y in sc_y ])

	# check if the clusters are within the overlap mask
	sc_in_mask = np.array([True if mask[y,x] == 1 else False for y,x in zip(sc_y_int, sc_x_int)])
	wfalse = np.where(sc_in_mask == False)[0]

	# drop clusters outside of the mask
	sc_ra  = np.delete(sc_ra, wfalse)
	sc_dec = np.delete(sc_dec, wfalse)
	sc_x   = np.delete(sc_x, wfalse)
	sc_y   = np.delete(sc_y, wfalse)

	sc_age = np.delete(sc_age, wfalse)

	# convert to SkyCoords
	sc_coords = SkyCoord(ra=sc_ra*u.deg, dec=sc_dec*u.deg, frame='icrs', distance=dist*u.Mpc)
	gmc_coords = SkyCoord(ra=gmc_ra*u.deg, dec=gmc_dec*u.deg, frame='fk5', distance=dist*u.Mpc)

	# get the 3 nearest neighbor separations and distances 
	nn3_sep, nn3_dist = sc_gmc_sep(sc_coords, gmc_coords, nn_number=3)

	# convert separation to arcsecs
	nn3_sep = (nn3_sep*u.deg).to(u.arcsec).value

	# histogram for 1st nearest neighbor
	filename = data_dir + '%s/%s_sc_gmc_1nn_hist'%(gal_name, gal_name)
	sc_gmc_sep_hist(nn3_sep[:,0], nn3_dist[:,0], age=sc_age, filename=filename, age_split=10, nn_label='1st', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')
	# histogram for 2nd nearest neighbor
	filename = data_dir + '%s/%s_sc_gmc_2nn_hist'%(gal_name, gal_name)
	sc_gmc_sep_hist(nn3_sep[:,1], nn3_dist[:,1], age=sc_age, filename=filename, age_split=10, nn_label='2nd', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')
	# histogram for 3rd nearest neighbor
	filename = data_dir + '%s/%s_sc_gmc_3nn_hist'%(gal_name, gal_name)
	sc_gmc_sep_hist(nn3_sep[:,2], nn3_dist[:,2], age=sc_age, filename=filename, age_split=10, nn_label='3rd', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')
