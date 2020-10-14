"""
main script for compiling the cluster catalog info and making class 1 and 2 catalogs and whatnot

2020-10-13
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
import astropy.units as u
import sys

sys.path.append('/cherokee1/turner/phangs/cf/utils')

from gmc_cat_utils import *

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

	# check if the catalog fits files exists
	try:
		cat = fits.open(data_dir + 'alma/%s/%s_12m+7m+tp_co21_native_props.fits'%(gal_name,gal_name))[1].data

	except FileNotFoundError:
		print(data_dir + 'alma/%s/%s_12m+7m+tp_co21_native_props.fits not found, skipping'%(gal_name,gal_name))
		continue

	n_total = len(cat)
	print('')
	print('%s total number of gmcs = %i'%(gal_name, n_total))
	print('')

	# pull out gmc coordinates 
	x			= cat['XCTR_PIX']
	y			= cat['YCTR_PIX']
	ra 			= cat['XCTR_DEG']
	dec			= cat['YCTR_DEG']
	pa_rad		= cat['POSANG'] * u.rad
	maj_pix		= cat['MOMMAJPIX']
	min_pix		= cat['MOMMINPIX']

	# will need to convert from pixels to arcsecs for major and minor axis
	# need image header
	mom0   = fits.open(data_dir + 'alma/%s/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))[0]
	header = mom0.header
	# degrees per pix
	cdelt = header['CDELT2']
	# convert to degrees
	maj_deg = maj_pix * cdelt * u.deg
	min_deg = min_pix * cdelt * u.deg
	# convert to arcsecs
	maj_asec = maj_deg.to(u.arcsec)
	min_asec = min_deg.to(u.arcsec)

	# convert postiion angle from radians to degrees
	pa_deg = pa_rad.to(u.deg)

	# dictionary
	coord_gmc = {'x': x, 'y': y, 'ra': ra, 'dec': dec, 'pa': pa_deg, 'maj_deg': maj_deg, 'min_deg': min_deg, 'maj_pix': maj_pix, 'min_pix': min_pix}

	# make ds9 region files for all the GMCs
	mk_gmc_ds9_regions(coord_gmc, data_dir + 'alma/%s/%s_gmc_cat'%(gal_name, gal_name), color='blue')
