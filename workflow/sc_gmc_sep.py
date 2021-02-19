"""
main script for running the code to get sc-gmc nearest neighbors

2021-01-04
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

import matplotlib
# non-interactive plots
matplotlib.use('agg')

# set path of the data dir
data_dir = '/cherokee1/turner/phangs/cf/data/'

# read in list of galaxies 
master_galaxy_list = ascii.read('master_galaxy.list')
galaxy_list = ascii.read('galaxy.list')

gal_id 		= galaxy_list['id']
gal_alt_id 	= galaxy_list['alt_id']
gal_dist 	= galaxy_list['dist']


mkhist = True

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
	gmc_cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_nativeres_nativenoise_props.fits'%(gal_name, gal_name))[1].data

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

	sc_cat_masked = np.delete(sc_cat, wfalse)

	# generate_gmc_cat_masked(galaxy_list, data_dir)

	# convert to SkyCoords
	sc_coords = SkyCoord(ra=sc_ra*u.deg, dec=sc_dec*u.deg, frame='icrs', distance=dist*u.Mpc)
	gmc_coords = SkyCoord(ra=gmc_ra*u.deg, dec=gmc_dec*u.deg, frame='fk5', distance=dist*u.Mpc)

	# get the nearest neighbor
	nn_sep, nn_dist, idx_sc, idx_gmc = sc_gmc_sep(sc_coords, gmc_coords, nn_number=1, return_idx=True)

	# convert separation to arcsecs
	nn_sep = (nn_sep*u.deg).to(u.arcsec).value

	# reshape the arrays
	nn_sep  = nn_sep.reshape(len(nn_sep))
	nn_dist = nn_dist.reshape(len(nn_dist))
	idx_sc  = idx_sc.reshape(len(idx_sc))
	idx_gmc = idx_gmc.reshape(len(idx_gmc))

	# get nearest neighbor gmc cloudnum and radius
	nn_gmc_cloudnum  = gmc_cat['CLOUDNUM'][idx_gmc]
	nn_gmc_radius_pc = gmc_cat['RAD_PC'][idx_gmc]

	# convert gmc radius from pc to arcsec
	nn_gmc_radius_radian = nn_gmc_radius_pc/(gal_dist[0]*1e6)
	nn_gmc_radius_asec = nn_gmc_radius_radian*u.rad.to(u.arcsec)

	# going to make a pandas dataframe for the star clusters within the overlap mask with:
	# ra, dec, hst_x, hst_y, alma_x, alma_y, sc_age + err, sc_mass + err, sc_ebv + err, environmental_mask_value, nn_gmc_cloudnum, nn_gmc_radius_pc + arcsec, nn_gmc_sep, nn_gmc_dist

	# need the star cluster pixel positions within the ALMA maps for easy placement within the environmental masks 
	alma_x, alma_y = get_alma_pixel_coords(sc_ra, sc_dec, gal_name, data_dir)

	# read in evironmental mask
	env_mask_hdulist = fits.open(data_dir + '%s/alma/%s_simple.fits'%(gal_name, gal_name))
	env_mask = env_mask_hdulist[0].data 
	env_mask_hdr = env_mask_hdulist[0].header

	# wcs for the environmenal mask
	w_env_mask = WCS(env_mask_hdr, fobj=env_mask_hdulist)
	# convert star cluster ra, dec to pixel locations in the environ mask
	env_mask_x, env_mask_y = w_env_mask.wcs_world2pix(sc_ra, sc_dec, 1.0)

	# convert to integers
	env_mask_x_int = np.array([np.round(x).astype(int) for x in env_mask_x])
	env_mask_y_int = np.array([np.round(y).astype(int) for y in env_mask_y])
	
	sc_env_mask_val = np.array([env_mask[y, x] for y,x in zip(env_mask_y_int, env_mask_x_int) ], dtype='int')

	# create dictionary with all the data which will go into the pandas dataframe
	df_data = {'id': sc_cat_masked['ID_PHANGS_ALLSOURCES_V0_9'], 'ra': sc_ra, 'dec': sc_dec, 'hst_x': sc_x, 'hst_y': sc_y, 'alma_x': alma_x, 'alma_y': alma_y,
			   'age': sc_cat_masked['PHANGS_AGE_MINCHISQ'], 'age_err': sc_cat_masked['PHANGS_AGE_MINCHISQ_ERR'], 
			   'mass': sc_cat_masked['PHANGS_MASS_MINCHISQ'], 'mass_err': sc_cat_masked['PHANGS_MASS_MINCHISQ_ERR'], 
			   'ebv': sc_cat_masked['PHANGS_EBV_MINCHISQ'], 'ebv_err': sc_cat_masked['PHANGS_EBV_MINCHISQ_ERR'],
			   'env_mask_val': sc_env_mask_val, 'nn_cloudnum': nn_gmc_cloudnum, 'nn_gmc_radius_pc': nn_gmc_radius_pc, 'nn_gmc_radius_asec': nn_gmc_radius_asec, 
			   'nn_sep_asec': nn_sep, 'nn_dist_pc': nn_dist }
	
	sc_df = Table(df_data).to_pandas()

	# output to a csv for easy manipulating later on
	sc_df.to_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12.csv'%(gal_name, gal_name), index=False)


	if mkhist:

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

