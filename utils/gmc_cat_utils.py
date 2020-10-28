"""
utilities/functions for manipulating the gmc catalogs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u


def mk_gmc_ds9_regions(coord_dict, filename, color='blue'):
	""" function to make ds9 region files for the gmc catalog 

	"""

	x = coord_dict['x']
	y = coord_dict['y']
	pa = coord_dict['pa'].value
	maj_pix = coord_dict['maj_pix']
	min_pix = coord_dict['min_pix']
	
	f = open(filename + '.%s'%color + '.pix.reg', 'w')

	f.write('image\n')
	for i in range(len(coord_dict['x'])):
		f.write('ellipse(%.7f,%.7f,%.12f,%.12f,%.7f) # color=%s\n'%(x[i], y[i], maj_pix[i], min_pix[i], pa[i], color) )

	f.close()

	ra = coord_dict['ra']
	dec = coord_dict['dec']
	maj_deg = coord_dict['maj_deg'].value
	min_deg = coord_dict['min_deg'].value


	f = open(filename + '.%s'%color + '.deg.reg', 'w')

	f.write('fk5\n')
	for i in range(len(coord_dict['x'])):
		f.write('ellipse(%.7f,%.7f,%.12f,%.12f,%.7f) # color=%s\n'%(ra[i], dec[i], maj_deg[i], min_deg[i], pa[i], color) )

	f.close()

def all_galaxies_gmc_region_files(galaxy_list, data_dir):
	""" loop through all the given galaxies and make ds9 region files for the GMCs

	"""

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
			cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits'%(gal_name,gal_name))[1].data
	
		except FileNotFoundError:
			print(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits not found, skipping'%(gal_name,gal_name))
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
		mom0   = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))[0]
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
		mk_gmc_ds9_regions(coord_gmc, data_dir + '%s/alma/%s_gmc_cat'%(gal_name, gal_name), color='blue')
