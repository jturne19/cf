"""
utilities/functions for manipulating the images and image related things

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.wcs import WCS
import glob

def flc_to_ds9_region(flc):
	""" create ds9 region file for the given flc file 

	"""

	hdulist = fits.open(flc)

	hdr1 = fits.open(flc)[1].header
	hdr4 = fits.open(flc)[4].header

	chip_x1 = 0.0
	chip_x2 = hdr1['naxis1']
	chip_x3 = hdr1['naxis1']
	chip_x4 = 0.0
	chip_y1 = 0.0
	chip_y2 = 0.0
	chip_y3 = hdr1['naxis2']
	chip_y4 = hdr1['naxis2']

	chip_x = [chip_x1,chip_x2,chip_x3,chip_x4]
	chip_y = [chip_y1,chip_y2,chip_y3,chip_y4]

	# convert the pixel location to WCS 
	# chip 1
	w = WCS(hdr1, fobj=hdulist)
	chip1_ra, chip1_dec = w.wcs_pix2world(chip_x, chip_y, 1)

	# 	# chip 2
	w = WCS(hdr4, fobj=hdulist)
	chip2_ra, chip2_dec = w.wcs_pix2world(chip_x, chip_y, 1)

	# output to ds9 region file

	out1 = 'polygon(%.8f,%.8f,'%(chip1_ra[0], chip1_dec[0]) + '%.8f,%.8f,'%(chip1_ra[1], chip1_dec[1]) + '%.8f,%.8f,'%(chip1_ra[2], chip1_dec[2]) + '%.8f,%.8f)\n'%(chip1_ra[3], chip1_dec[3])
	out2 = 'polygon(%.8f,%.8f,'%(chip2_ra[0], chip2_dec[0]) + '%.8f,%.8f,'%(chip2_ra[1], chip2_dec[1]) + '%.8f,%.8f,'%(chip2_ra[2], chip2_dec[2]) + '%.8f,%.8f)\n'%(chip2_ra[3], chip2_dec[3])

	f = open('%s.deg.reg'%flc[:-9], 'w')
	f.write('fk5\n')
	f.write(out1)
	f.write(out2)
	f.close()

	return out1, out2


def all_galaxies_flc_region_files(galaxy_list, data_dir):
	""" loop through all the given galaxies and make ds9 region file for the 
	flc files which mark the footprints of the HST observations


	"""
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']

	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
	
		print(gal_name)

		flc_list = glob.glob(data_dir + '%s/hst/flc/*_flc.fits'%(gal_name))
		flc_list = np.sort(flc_list)

		f = open(data_dir + '%s/hst/flc/%s_all_flc.deg.reg'%(gal_name,gal_name), 'w')
		f.write('fk5\n')

		for flc in flc_list:

			out1, out2 = flc_to_ds9_region(flc)

			f.write(out1)
			f.write(out2)

		f.close()

