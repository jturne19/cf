"""
master file which contains all the utilities/functions 

"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import aplpy as ap

from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u

from scipy.optimize import curve_fit

from astroML.decorators import pickle_results
from astroML.correlation import bootstrap_two_point_angular

########################################################
# cluster catalog specific functions
########################################################

def new_clust_fits_from_classes(base_cat, classes, filename):
	""" takes input base catalog, grabs clusters with only the given
	classes, and outputs new fits file 
		
	Inputs:
	base_cat	astropy.io.fits.fitsrec.FITS_rec 	PHANGS fits cluster catalog file read in by astropy
	classes 	list								list with cluster classes wanted e.g. [1,2,3]
	filename 	str 								filename of the output fits file

	"""

	# grab just the cluster classes
	cc = base_cat['PHANGS_CLUSTER_CLASS']

	# grab index of clusters with the input classes regardless of the number of wanted classes
	w_unordered = np.array([np.where(cc == cl)[0] for cl in classes])
	w_unordered = np.concatenate(w_unordered).ravel()

	# sort
	w = np.sort(w_unordered)

	# make astropy table with the chosen clusters
	t = Table(base_cat[w])

	# output
	try:
		t.write(filename, overwrite=False)
	except:
		print(filename + ' already exists, exiting')


def all_galaxies_clust_cats(galaxy_list, data_dir, run_name, sc_base_cat_suffix='_phangshst_base_catalog',
							gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', sc_class='class12', mkfits=False):
	""" loop through all the given glaxies and print out the cluster numbers
	and make new fits tables with the chosen sc_class 

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	sc_base_cat_suffix	str 			suffix for the filename of the base star cluster catalog;
										included here to make it easier if new cluster catalogs are generated [and the name changes i guess?]
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	sc_class 			str 			which class of clusters to make the catalogs for; class12 or class123
	mkfits 				bool 			if true, then fits tables with the cleaned cluster catalogs will be saved;
										saves them as sc_cat_suffix + '.class12.fits' and sc_cat_suffix + '.class123.fits'
	"""
	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))

		# check if the base catalog fits files exists
		try:
			base_cat = fits.open(data_dir + '%s/hst/%s%s.fits'%(gal_name,gal_name,sc_base_cat_suffix))[1].data
	
		except FileNotFoundError:
			print(data_dir + '%s/hst/%s%s.fits not found, skipping'%(gal_name,gal_name,sc_base_cat_suffix))
			continue
		
		n_total = len(base_cat)
		print('total number clusters = %i'%(n_total))
		
		if sc_class == 'class123':
			classes = [1,2,3]
		else:
			classes = [1,2]

		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in base_cat.names:
			
			# make new fits files with the restricted classes
			if mkfits:
				new_clust_fits_from_classes(base_cat, classes, data_dir + '%s/%s/%s%s.%s.fits'%(gal_name,run_name,gal_name,sc_base_cat_suffix,sc_class))
	
			# read those in
			sccat = fits.open(data_dir + '%s/%s/%s%s.%s.fits'%(gal_name,run_name,gal_name,sc_base_cat_suffix,sc_class))[1].data
						
			print('number of %s = %i'%(sc_class, len(sccat)))
			print('')


def mk_clust_ds9_regions(coord_dict, radius_pix, radius_asec, filename, color='green'):
	""" function to make ds9 region files for the gmc catalog 

	"""
	x = coord_dict['x']
	y = coord_dict['y']

	# hst image pixel region file
	f = open(filename + '.%ipix.%s.pix.reg'%(radius_pix,color), 'w')

	f.write('image\n')
	for i in range(len(coord_dict['x'])):
		f.write('circle(%.7f,%.7f,%i) # color=%s\n'%(x[i], y[i], radius_pix, color))
	f.close()

	ra = coord_dict['ra']
	dec = coord_dict['dec']

	# ra, dec region file
	f = open(filename + '.%ipix.%s.deg.reg'%(radius_pix,color), 'w')

	f.write('fk5\n')
	for i in range(len(coord_dict['x'])):
		f.write('circle(%.7f,%.7f,%.3f") # color=%s\n'%(ra[i], dec[i], radius_asec, color) )

	f.close()

def all_galaxies_clust_region_files(galaxy_list, data_dir, run_name, sc_base_cat_suffix='_phangshst_base_catalog',
									gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', 
									sc_class='class12', radius_pix=10):
	""" loop through all the given galaxies and make ds9 region files for the star clusters
	makes separate region files for class 1,2,3 and class 1,2 if available
	makes hst image pixel region files and degree region files
	
	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	sc_base_cat_suffix	str 			suffix for the filename of the base star cluster catalog;
										included here to make it easier if new cluster catalogs are generated [and the name changes i guess?]
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	sc_class 			str 			which class of clusters to make the catalogs for; class12 or class123 
	radius_pix 			int 			sets the ds9 circle region radius in pixels; default is 10 pixels

	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))

		# check if the catalog fits files exists
		try:
			cat = fits.open(data_dir + '%s/hst/%s%s.fits'%(gal_name,gal_name,sc_base_cat_suffix))[1].data
		except FileNotFoundError:
			print(data_dir + '%s/hst/%s%s.fits not found, skipping'%(gal_name,gal_name,sc_base_cat_suffix))
			continue
	
		# pull out sc coordinates 
		x	= cat['PHANGS_X']
		y	= cat['PHANGS_Y']
		ra 	= cat['PHANGS_RA']
		dec	= cat['PHANGS_DEC']
	
		coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
		# read in image header to get cdelt so we can convert pixel to arcsec
		# use 275 because its UVIS for all of them
		image_hdr = fits.getheader(data_dir + '%s/hst/%s_uvis_f275w_exp_drc_sci.fits'%(gal_name, gal_name))
	
		# degrees per pix
		cdelt = image_hdr['CD2_2']
		# convert input pixel radius to degrees
		radius_deg = radius_pix * cdelt * u.deg
		# conver to arcsec
		radius_asec = radius_deg.to(u.arcsec).value
	
		mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + '%s/%s/%s_allclusters'%(gal_name,run_name,gal_name), color='green')
	
		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in cat.names:
			
			sccat = fits.open(data_dir + '%s/%s/%s%s.%s.fits'%(gal_name,run_name,gal_name,sc_base_cat_suffix,sc_class))[1].data
			
			x	= sccat['PHANGS_X']
			y	= sccat['PHANGS_Y']
			ra 	= sccat['PHANGS_RA']
			dec	= sccat['PHANGS_DEC']
			coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
			mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + '%s/%s/%s_%s'%(gal_name,run_name,gal_name,sc_class), color='green')

########################################################
# GMC catalog specific functions
########################################################

def mk_gmc_ds9_regions(coord_dict, filename, color='blue'):
	""" function to make ds9 region (both pixel and degree) files for the given gmc catalog 
	
	Inputs:
	coord_dict 		dictionary 		dictionary which holds the GMC x,y positions, position angle, major and minor axis in pixels
	filename 		str 			path and filename to name the output ds9 region files
	color 			str 			color for the ds9 regions; default is blue

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

def all_galaxies_gmc_region_files(galaxy_list, data_dir, run_name, gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props'):
	""" loop through all the given galaxies and make ds9 region files for the GMCs

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog

	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']
	
	# loop through all galaxies
	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))

		# check if the catalog fits files exists
		try:
			cat = fits.open(data_dir + '%s/alma/%s%s.fits'%(gal_name, gal_name, gmc_cat_suffix))[1].data
	
		except FileNotFoundError:
			print(data_dir + '%s/alma/%s%s.fits not found, skipping'%(gal_name, gal_name, gmc_cat_suffix))
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
		mk_gmc_ds9_regions(coord_gmc, data_dir + '%s/%s/%s_gmc_cat'%(gal_name, run_name, gal_name), color='blue')

def generate_gmc_cat_masked(galaxy_list, data_dir, run_name, gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props'):
	""" function to generate the gmc catalog for just the gmcs which lie within the
	HST-ALMA footprint overlap mask

	"""
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):

		gal_name = gal_id[i]
		print('')
		print(gal_name)

		# read in the full gmc cat
		gmc_cat = fits.open(data_dir + '%s/alma/%s%s.fits'%(gal_name, gal_name, gmc_cat_suffix))[1].data

		# grab gmc ra, dec
		ra, dec = gmc_cat['XCTR_DEG']*u.deg, gmc_cat['YCTR_DEG']*u.deg
		coords = SkyCoord(ra, dec, frame='fk5', distance=gal_dist[i])

		# read in the footprint overlap mask
		mask_hdulist = fits.open(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name))
		mask_header = mask_hdulist[0].header
		mask = mask_hdulist[0].data

		# mask/hst wcs
		w_mask = WCS(mask_header, fobj=mask_hdulist)

		# convert gmc ra, dec to mask pixel locations
		mask_x, mask_y = w_mask.wcs_world2pix(ra, dec, 1.0)

		mask_x_int = np.array([np.round(x).astype(int) for x in mask_x])
		mask_y_int = np.array([np.round(y).astype(int) for y in mask_y])

		in_mask =np.array([True if mask[y,x] == 1 else False for y,x in zip(mask_y_int, mask_x_int)])

		gmc_cat_in_mask = gmc_cat[in_mask]

		# output
		t = Table(gmc_cat_in_mask)
		t.write(data_dir + '%s/%s/%s_gmc_cat_masked.fits'%(gal_name, run_name, gal_name), format='fits', overwrite=True)




########################################################
# fits image specific functions
########################################################

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


def get_alma_pixel_coords(ra, dec, gal_name, data_dir, return_int=False):
	""" given the ra and dec of an object, return the pixel position in the ALMA maps

	"""
	# read in alma map hdulist
	alma_hdulist = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))
	alma_header =alma_hdulist[0].header

	# get alma map wcs 
	w_alma = WCS(alma_header, fobj=alma_hdulist, naxis=2)
	# convert ra, dec to alma map pixels
	alma_x, alma_y = w_alma.wcs_world2pix(ra, dec, 1.0)

	if return_int:
		# integer pixel location if requested
		alma_x_int = np.array([np.round(x).astype(int) for x in alma_x])
		alma_y_int = np.array([np.round(y).astype(int) for y in alma_y])

		return alma_x_int, alma_y_int

	else:
		
		return alma_x, alma_y


def generate_overlap_mask(galaxy_list, data_dir):
	""" function to go through the given galaxies and 
	create a 'mask' which identifies where the HST footprints
	overlap with the ALMA footprints
	which is needed for defining the region to populate the 
	random data points for the correlation functions

	"""
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']

	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
		print('')
		print(gal_name)

		# read in hst 555 (V-band) fits image
		hst_hdulist = fits.open(data_dir + '%s/hst/%s_uvis_f555w_exp_drc_sci.fits'%(gal_name,gal_name))
		hst_data = hst_hdulist[0].data
		hst_header = hst_hdulist[0].header

		# wcs coordinates
		w_hst = WCS(hst_header, fobj=hst_hdulist)

		# grid of pixel locations
		hst_pixel_grid = np.meshgrid(np.arange(0,len(hst_data[0,:]),1), np.arange(0,len(hst_data[:,0]),1))

		# convert pixels to wcs 
		ra, dec = w_hst.wcs_pix2world(hst_pixel_grid[0], hst_pixel_grid[1], 0.5)

		# read in moment 0 ALMA map
		alma_hdulist = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))
		alma_data = alma_hdulist[0].data
		alma_header = alma_hdulist[0].header

		# get the ALMA wcs info
		w_alma = WCS(alma_header, fobj=alma_hdulist, naxis=2)

		# convert the ra, dec to alma pixel locations
		alma_x, alma_y = w_alma.wcs_world2pix(ra, dec, 0.5)

		# round the float pixel locations to integers which will be used for indexing the ALMA data values
		alma_x_int = np.array([np.round(x).astype(int) for x in alma_x])
		alma_y_int = np.array([np.round(y).astype(int) for y in alma_y])

		# create array of nans the same shape as the hst data which will hold the ALMA data
		# corresponding to the HST pixels
		# i.e., mapping HST pixels to the larger ALMA ones
		alma_data_hst_pix = np.zeros(np.shape(hst_data)) + np.nan

		# find the maximum ALMA pixel locations
		alma_max_x = len(alma_data[0,:])
		alma_max_y = len(alma_data[:,0])

		# create dictionary with all our data which will go into the pandas data frame
		# HST pixel x, y locations, HST pixel ra, dec locations, HST pixel locations in the ALMA map,
		# HST pixel value, corresponding ALMA pixel value,
		# mask column is the whole point ---> 1 will mean the footprints overlap
		pixels_data = {'hst_x': hst_pixel_grid[0].flatten().astype('int16'), 'hst_y': hst_pixel_grid[1].flatten().astype('int16'), 
					   'ra': ra.flatten(), 'dec': dec.flatten(), 'alma_x': alma_x_int.flatten().astype('int16'), 'alma_y': alma_y_int.flatten().astype('int16'), 
					   'hst_val': hst_data.flatten(), 'alma_val': alma_data_hst_pix.flatten(), 'mask': np.zeros(len(hst_data.flatten()),dtype='int8')}

		# something's up with going from astropy fits read-in to pandas so have to make an astropy table first
		# and use the built-in function to convert to a pandas dataframe
		pixels_df = Table(pixels_data).to_pandas()

		# find the pixels which lie within the ALMA footprint because that's usually smaller than HST
		tmp_df = pixels_df.loc[(pixels_df['alma_x'] > 0) & (pixels_df['alma_x'] <= alma_max_x) & (pixels_df['alma_y'] > 0) & (pixels_df['alma_y'] <= alma_max_y)]
		# pull out the ALMA x, y pixel locations which act as idices for the ALMA data
		wy_alma = tmp_df['alma_y'].to_numpy() - 1
		wx_alma = tmp_df['alma_x'].to_numpy() - 1

		# pull out the ALMA values for the interested pixels
		alma_data_replace = alma_data[(wy_alma, wx_alma)]
		# replace this into the dataframe for the alma_val column
		pixels_df.loc[(pixels_df['alma_x'] > 0) & (pixels_df['alma_x'] <= alma_max_x) & (pixels_df['alma_y'] > 0) & (pixels_df['alma_y'] <= alma_max_y), ('alma_val')] = alma_data_replace
		# find where the HST value is not 0 and the ALMA value is not nan
		# this determines the overlap in footprint so mask is set to 1
		pixels_df.loc[(pixels_df['hst_val'] != 0) & (pixels_df['alma_val'].isnull() == False), ('mask')] = 1

		# pull out mask column
		mask = pixels_df['mask'].to_numpy()
		# reshape into the original HST image shape
		mask = np.reshape(mask, np.shape(hst_data)).astype(int)
		# write to a fits image using the HST header 
		hdu = fits.PrimaryHDU(mask, header=hst_header)
		hdu.writeto(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name), overwrite=True)


########################################################
# making figures/plots functions
########################################################

def outline_plot(gal_name, data_dir, run_name, gmc_cat_suffix, sc_base_cat_suffix, sc_coord_dict, center_deg, sc_class='class12', radius=0.04, bkgd=None, color_arr=[], color_code='' ):
	""" create an 'outline plot' for the given galaxy
	plot using the wcs info in the HST images
	shows the ellipses of the GMCs 
	shows the position of the star clusters 
	
	Inputs:
	gal_name		str 		name of the galaxy e.g., ngc0628
	data_dir 		str 		path to the directory containing the 'hst' and 'alma' directories
	run_name 			str 			name of the run/test; e.g., run01
	sc_base_cat_suffix	str 			suffix for the filename of the base star cluster catalog;
										included here to make it easier if new cluster catalogs are generated [and the name changes i guess?]
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	sc_coord_dict 	dict 		dictionary containing the star cluster coordinates: 'x', 'y', 'ra', 'dec'
	center_deg 		tuple 		(RA,DEC) coordinate of the center of the galaxy in degrees
	sc_class 		str 		chosen star cluster classes i.e., class12 or class123; needed for filename saving
	radius 			float		radius of the image view in degrees; default = 0.04
	bgkd 			bool		whether to include the 3-color image as the background
	color_arr 		list/arr	if included, the cluster data points will be colored by the given array i.e., cluster age
	color_code		str 		specifies the property given in color_arr; expects 'age' or 'mass'


	"""
	# read in GMC catalog
	# not sure if it's needed or not, probably just need the GMC ellipses
	gmc_cat = fits.open(data_dir + '%s/alma/%s%s.fits'%(gal_name,gal_name,gmc_cat_suffix))[1].data
	
	# path to GMC ellipses/region file
	gmc_region_file = data_dir + '%s/%s/%s_gmc_cat.blue.deg.reg'%(gal_name,run_name,gal_name)

	# path to HST image file 
	hst_image = data_dir + '%s/hst/%s_uvis_f555w_exp_drc_sci.fits'%(gal_name,gal_name)

	# path to full HST footprint file
	footprint_file = data_dir + '%s/hst/%s_full_footprint.pix.reg'%(gal_name,gal_name)


	if len(color_arr) > 0:
		fig = plt.figure(figsize=(8,6))
		f = ap.FITSFigure(hst_image, figure=fig)
	else:
		f = ap.FITSFigure(hst_image, figsize=(8,8))

	f.recenter(center_deg[0], center_deg[1], radius=radius)

	if bkgd != None:
		trilogy_image = data_dir + '%s/hst/png/%s_trilogy.png'%(gal_name, gal_name)
		f.show_rgb(trilogy_image)

	f.show_regions(gmc_region_file, layer='gmc', zorder=2)
	f.show_regions(footprint_file, layer='footprint', zorder=0)

	if len(color_arr) > 0:

		cm = plt.cm.get_cmap('spring_r')
		f.show_markers(sc_coord_dict['ra'], sc_coord_dict['dec'], layer='clusters', c=color_arr, cmap=cm, marker='o', s=30, edgecolor='black', lw=0.5, zorder=1)
		
		ax2= fig.add_axes([0.1,0.1,.8,.8], sharex=fig.get_axes()[0], sharey=fig.get_axes()[0], frameon=False)
		ax2.axis('off')
		scat = ax2.scatter(sc_coord_dict['ra'], sc_coord_dict['dec'], c=color_arr, cmap=cm, s=0) 
		cbar0 = fig.colorbar(scat, ax=fig.get_axes()[0], aspect=35, pad=0.03)

		cbar0.set_label(r'log(%s) [Myr]'%color_code, fontsize='large', labelpad=10)
	else:
		f.show_markers(sc_coord_dict['ra'], sc_coord_dict['dec'], layer='clusters', marker='o', s=15, facecolor='red')

	f.axis_labels.hide_y()
	f.axis_labels.hide_x()

	f.tick_labels.set_font(size='large', family='serif')
	f.set_theme('publication')
	f.tick_labels.set_xformat('hh:mm:ss')
	f.tick_labels.set_yformat('dd:mm:ss')

	if bkgd != None:
		
		f.save(data_dir + '%s/%s/%s_%s_outlineplot_%s.png'%(gal_name, run_name, gal_name, sc_class, bkgd))
		f.save(data_dir + '%s/%s/%s_%s_outlineplot_%s.pdf'%(gal_name, run_name, gal_name, sc_class, bkgd))
	
	elif len(color_arr) > 0:
		
		fig.savefig(data_dir + '%s/%s/%s_%s_outlineplot_%s.png'%(gal_name, run_name, gal_name, sc_class, color_code), bbox_inches='tight')
		fig.savefig(data_dir + '%s/%s/%s_%s_outlineplot_%s.pdf'%(gal_name, run_name, gal_name, sc_class, color_code), bbox_inches='tight')
	
	else:
		
		f.save(data_dir + '%s/%s/%s_%s_outlineplot.png'%(gal_name, run_name, gal_name, sc_class))
		f.save(data_dir + '%s/%s/%s_%s_outlineplot.pdf'%(gal_name, run_name, gal_name, sc_class))


	plt.close()


def all_galaxies_outline_plots(galaxy_list, data_dir, run_name, sc_base_cat_suffix='_phangshst_base_catalog',
							   gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', sc_class='class12', 
							   radius=[0.04], bkgd=None, color_code=''):
	""" loop through all the galaxies in the list and make the outline plots 

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	sc_base_cat_suffix	str 			suffix for the filename of the base star cluster catalog;
										included here to make it easier if new cluster catalogs are generated [and the name changes i guess?]
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	sc_class 			str 			which class clusters will be plotted? 'class12' or 'class123'; default is class12
	radius 				list 			list same length as galaxy_list which has radius of the image view in degrees; default is 0.04
	bkgd 				str 			if you want the trilogy 3-color image as background image, put 'trilogy' otherwise default is None
	color_code 			str 			if you want to color code the data points by age, put 'age', if by mass, put 'mass'; default is nothing so nothing color-coded

	"""
	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']
	gal_ra		= galaxy_list['center_ra']
	gal_dec 	= galaxy_list['center_dec']

	for i in range(len(galaxy_list)):

		gal_name = gal_id[i]
		center_deg = (gal_ra[i], gal_dec[i])
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))

		# read in the correct cluster catalog
		if 'class12' in sc_class:
			sc_cat = fits.open(data_dir + '%s/%s/%s%s.%s.fits'%(gal_name,run_name,gal_name,sc_base_cat_suffix,sc_class))[1].data
		else:
			# default to the class 1,2 catalog
			sc_cat = fits.open(data_dir + '%s/%s/%s%s.class12.fits'%(gal_name,run_name,gal_name,sc_base_cat_suffix))[1].data

		# create dictionary containing the star cluster coordinates to feed to the outline_plot function
		# positions
		x_sc   = sc_cat['PHANGS_X']
		y_sc   = sc_cat['PHANGS_Y']
		ra_sc  = sc_cat['PHANGS_RA']
		dec_sc = sc_cat['PHANGS_DEC']

		sc_coord_dict = {'x': x_sc, 'y': y_sc, 'ra': ra_sc, 'dec': dec_sc}

		gal_plot_radius = radius[i]

		if color_code == 'age':
			# read in cluster age
			color_input = np.log10(sc_cat['PHANGS_AGE_MINCHISQ'])
		elif color_code == 'mass':
			color_input = np.log10(sc_cat['PHANGS_MASS_MINCHISQ'])
		else:
			color_input = []

		outline_plot(gal_name, data_dir, run_name, gmc_cat_suffix, sc_base_cat_suffix, sc_coord_dict, center_deg, sc_class=sc_class, radius=gal_plot_radius, bkgd=bkgd, color_arr=color_input, color_code=color_code)


########################################################
# star cluster - gmc separation functions
########################################################

def sc_gmc_sep(sc_coords, gmc_coords, nn_number=3, return_idx=False):
	"""	takes the given star clusters and GMCs finds
	the nearest GMC to each star cluster i.e., nearest neighbor
	and returns the separation

	"""
	# search around basically full footprint to make sure we get at least the 3 nearest neighbors
	idx_sc, idx_gmc, sep, dist3d = search_around_sky(sc_coords, gmc_coords, 5*u.arcmin)

	nn_seps    = []
	nn_dists   = []
	nn_idx_sc  = []
	nn_idx_gmc = []

	# loop through all the star clusters
	for i in range(len(sc_coords)):

		wi = np.where(idx_sc == i)[0]

		sci   = idx_sc[wi]
		gmci  = idx_gmc[wi]
		sepi  = sep[wi]
		disti = dist3d[wi]

		wsort = np.argsort(sepi)
		# separation (deg) for the given nearest neighbor number
		nn_sep  = sepi[wsort][:nn_number].value
		nn_dist = disti[wsort][:nn_number].to(u.pc).value

		nn_seps.append(nn_sep)
		nn_dists.append(nn_dist)

		nn_idx_sc.append(sci[wsort][:nn_number])
		nn_idx_gmc.append(gmci[wsort][:nn_number])

	nn_seps = np.array(nn_seps)
	nn_dists = np.array(nn_dists)
		
	nn_idx_sc = np.array(nn_idx_sc)
	nn_idx_gmc = np.array(nn_idx_gmc)


	if return_idx:
		
		return nn_seps, nn_dists, nn_idx_sc, nn_idx_gmc

	else:

		return nn_seps, nn_dists


def sc_gmc_sep_hist(sep, dist, age, filename, age_split=10, nn_label='1st', sep_unit='arcsec', **kwargs):
	""" histograms for the given nearest neighbor separations and distances
	
	"""

	wyoung = np.where(age <= age_split)
	wold   = np.where(age > age_split)

	all_med = np.median(sep)
	young_med = np.median(sep[wyoung])
	old_med = np.median(sep[wold])

	fig, ax1 = plt.subplots(1,1, figsize=(5,5))

	ax1.hist(sep, label='All Clusters', histtype='step', **kwargs)
	ax1.axvline(all_med, ls='--', color='black', lw=1, zorder=0)
	
	ax1.hist(sep[wyoung], label=r'$\leq$ %i Myr'%age_split, histtype='step', edgecolor='#377eb8', lw=2.5, zorder=1)
	ax1.axvline(young_med, ls='--', color='#377eb8', lw=1, zorder=1)

	ax1.hist(sep[wold], label=r'$>$ %i Myr'%age_split, histtype='step', edgecolor='#e41a1c', lw=2.5, zorder=2)
	ax1.axvline(old_med, ls='--', color='#e41a1c', lw=1, zorder=2)

	ax2 = ax1.twiny()
	ax2.hist(dist, histtype='step', lw=0)

	ax1.set_xlabel('%s Nearest Neighbor Separation (%s)'%(nn_label, sep_unit))
	ax2.set_xlabel('Distance (pc)')

	ax1.set_ylabel('Number')

	ax1.legend(loc='upper right', fontsize='small')

	plt.savefig(filename + '.png', bbox_inches='tight')
	plt.savefig(filename + '.pdf', bbox_inches='tight')
	plt.close()

def all_galaxies_sc_gmc_sep(galaxy_list, data_dir, run_name, mkhist=True, sc_cat_suffix='_phangshst_base_catalog.class12', 
							gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', output_cat_suffix='_cluster_catalog_in_mask_class12'):
	""" nearest-neighbor analysis by finding the separations between star clusters and gmcs
	function version of sc_gmc_sep.py in the workflow directory

	will output a new cluster catalog as as csv which contains the clusters within the overlap mask and the information about the nearest neighbor
	gmcs (cloud number, radius, separations) and which environment the cluster lies in according to the environmental masks
	saves to the given run_name directory within each galaxy's data dir
	also can make histograms of the nearest neighbor separations - saves to the same directory as above

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	mkhist 				bool 			if true, then histograms of the 1st, 2nd, 3rd nearest-neighbor separations will be produced
	sc_cat_suffix 		str 			suffix for the filename of the star cluster catalog; defaults to the class 1 and 2 catalog;
										included here to make it easier to do the runs with class 1,2,3 and the associations
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	output_cat_suffix 	str 			suffix for the filename of the output star cluster catalog which has the new cluster catalog;
										defaults to the class 1,2 name
	"""
	# read in the galaxy list
	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']

	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):

		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))

		# read in star cluster catalog
		sc_cat = fits.open(data_dir + '%s/hst/%s%s.fits'%(gal_name, gal_name, sc_cat_suffix))[1].data

		# grab star cluster positions
		sc_x, sc_y,   = sc_cat['PHANGS_X'], sc_cat['PHANGS_Y']
		sc_ra, sc_dec = sc_cat['PHANGS_RA'], sc_cat['PHANGS_DEC']
		# grab star cluster ages
		sc_age = sc_cat['PHANGS_AGE_MINCHISQ']

		# read in GMC catalog 
		gmc_cat = fits.open(data_dir + '%s/alma/%s%s.fits'%(gal_name, gal_name, gmc_cat_suffix))[1].data

		# not sure why but some of the GMCs have nan radius so we drop those from the catalog
		wnanrad = np.where(np.isnan(gmc_cat['RAD3D_PC']))[0]
		gmc_cat = np.delete(gmc_cat, wnanrad)

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

		# keep only clusters within the mask
		sc_ra  = sc_ra[sc_in_mask]
		sc_dec = sc_dec[sc_in_mask]
		sc_x   = sc_x[sc_in_mask]
		sc_y   = sc_y[sc_in_mask]

		sc_age   = sc_age[sc_in_mask]
		
		sc_cat_masked   = sc_cat[sc_in_mask]

		# convert to SkyCoords
		sc_coords = SkyCoord(ra=sc_ra*u.deg, dec=sc_dec*u.deg, frame='icrs', distance=dist*u.Mpc)
		gmc_coords = SkyCoord(ra=gmc_ra*u.deg, dec=gmc_dec*u.deg, frame='fk5', distance=dist*u.Mpc)

		# get the nearest neighbor
		nn_sep, nn_dist, idx_sc, idx_gmc = sc_gmc_sep(sc_coords, gmc_coords, nn_number=1, return_idx=True)
	
		# convert separation from degrees to arcsecs
		nn_sep = (nn_sep*u.deg).to(u.arcsec).value
	
		# reshape the arrays
		nn_sep  = nn_sep.reshape(len(nn_sep))
		nn_dist = nn_dist.reshape(len(nn_dist))
		idx_sc  = idx_sc.reshape(len(idx_sc))
		idx_gmc = idx_gmc.reshape(len(idx_gmc))

		# get nearest neighbor gmc cloudnum, radius, and mass [luminosity mass]
		# R_3D for the cloud radius as this takes into account the galaxy scale height (assumes 100 pc) 
		# which restricts cloud sizes within the disk of the galaxy
		nn_gmc_cloudnum  = gmc_cat['CLOUDNUM'][idx_gmc]
		nn_gmc_radius_pc = gmc_cat['RAD3D_PC'][idx_gmc]
		nn_gmc_mlum 	 = gmc_cat['MLUM_MSUN'][idx_gmc]

		# convert gmc radius from pc to arcsec
		nn_gmc_radius_radian = nn_gmc_radius_pc/(gal_dist[0]*1e6)
		nn_gmc_radius_asec   = nn_gmc_radius_radian*u.rad.to(u.arcsec)

		# going to make a pandas dataframe for the star clusters within the overlap mask with:
		# ra, dec, hst_x, hst_y, alma_x, alma_y, sc_age + err, sc_mass + err, sc_ebv + err, environmental_mask_value, nn_gmc_cloudnum, nn_gmc_mlum, nn_gmc_radius_pc + arcsec, nn_gmc_sep, nn_gmc_dist

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
		df_data = {'id': sc_cat_masked['ID_PHANGS_ALLSOURCES_V0_9'], 'class': sc_cat_masked['PHANGS_CLUSTER_CLASS'], 'ra': sc_ra, 'dec': sc_dec, 'hst_x': sc_x, 'hst_y': sc_y, 'alma_x': alma_x, 'alma_y': alma_y,
				   'age': sc_cat_masked['PHANGS_AGE_MINCHISQ'], 'age_err': sc_cat_masked['PHANGS_AGE_MINCHISQ_ERR'], 
				   'mass': sc_cat_masked['PHANGS_MASS_MINCHISQ'], 'mass_err': sc_cat_masked['PHANGS_MASS_MINCHISQ_ERR'], 
				   'ebv': sc_cat_masked['PHANGS_EBV_MINCHISQ'], 'ebv_err': sc_cat_masked['PHANGS_EBV_MINCHISQ_ERR'],
				   'env_mask_val': sc_env_mask_val, 'nn_cloudnum': nn_gmc_cloudnum, 'nn_gmc_mlum': nn_gmc_mlum, 'nn_gmc_radius_pc': nn_gmc_radius_pc, 'nn_gmc_radius_asec': nn_gmc_radius_asec, 
				   'nn_sep_asec': nn_sep, 'nn_dist_pc': nn_dist }
		
		sc_df = Table(df_data).to_pandas()
	
		# output to a csv for easy manipulating later on
		sc_df.to_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, output_cat_suffix), index=False)

		if mkhist:
			
			# get the 3 nearest neighbor separations and distances 
			nn3_sep, nn3_dist = sc_gmc_sep(sc_coords, gmc_coords, nn_number=3)

			# convert separation to arcsecs
			nn3_sep = (nn3_sep*u.deg).to(u.arcsec).value

			# histogram for 1st nearest neighbor
			filename = data_dir + '%s/%s/%s_sc_gmc_1nn_hist'%(gal_name, run_name, gal_name)
			sc_gmc_sep_hist(nn3_sep[:,0], nn3_dist[:,0], age=sc_age, filename=filename, age_split=10, nn_label='1st', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')
			# histogram for 2nd nearest neighbor
			filename = data_dir + '%s/%s/%s_sc_gmc_2nn_hist'%(gal_name, run_name, gal_name)
			sc_gmc_sep_hist(nn3_sep[:,1], nn3_dist[:,1], age=sc_age, filename=filename, age_split=10, nn_label='2nd', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')
			# histogram for 3rd nearest neighbor
			filename = data_dir + '%s/%s/%s_sc_gmc_3nn_hist'%(gal_name, run_name, gal_name)
			sc_gmc_sep_hist(nn3_sep[:,2], nn3_dist[:,2], age=sc_age, filename=filename, age_split=10, nn_label='3rd', sep_unit='arcsec', bins=10, lw=2.5, edgecolor='black')

def bootstrap_median_error(data, sigma, nbootstraps=10000):
	"""	boostrap estimate of the error on the median of the distribution of input data and their errors/sigma
		
	"""

	medians = np.zeros(nbootstraps)

	for i in range(nbootstraps):

		rand = np.random.normal(data, sigma)

		medians[i] = np.median(rand)

	std = np.std(medians)

	return std

def generate_sc_gmc_assoc_df(galaxy_list, data_dir, run_name, sc_mask_cat_suffix='_cluster_catalog_in_mask_class12', 
							 gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', filename='_cluster_catalog_in_mask_class12_assoc_gmc' ):
	""" generates the dataframe (as a csv file) which adds to the dataframe/csv generated in all_galaxies_sc_gmc_sep
	adds columns: 
	cloud number of the gmc of which the star cluster is associated [nan if no association]
	the 'association number' i.e., 0 = unassociated, 1 = w/in 1 gmc rad, 2 = b/w 1 and 2 gmc rad, 3 = b/w 2 and 3 gmc rad

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	sc_mask_cat_suffix	str 			suffix for the filename of the star cluster within the overlap mask catalog; defaults to the class 1 and 2 catalog;
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	filename		 	str 			suffix for the filename of the output star cluster dataframe; defaults to the class 1,2 name

	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']

	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):
	
		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# check if the run_name directory exists and if not, create it
		if not os.path.exists(data_dir + '%s/%s'%(gal_name, run_name)):
			os.makedirs(data_dir + '%s/%s'%(gal_name, run_name))
	
		# read in the csv of the cluster catalog within the overlap mask
		sc_df = pd.read_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, sc_mask_cat_suffix))
	
		# read in the gmc cat
		gmc_cat = fits.open(data_dir + '%s/alma/%s%s.fits'%(gal_name, gal_name, gmc_cat_suffix))[1].data
	
		# get the star cluster ra, dec
		sc_ra, sc_dec = sc_df['ra'].to_numpy()*u.deg, sc_df['dec'].to_numpy()*u.deg
		sc_coords = SkyCoord(sc_ra, sc_dec, frame='icrs', distance=dist*u.Mpc)
	
		# initialize empty lists which will keep running lists of the gmcs and associated star clusters
		# assoc_track will keep track of the 'association number' which means:
		# 1 = w/in 1 gmc rad, 2 = b/w 1 and 2 gmc rad, 3 = b/w 2 and 3 gmc rad
		gmc_track = []
		sc_track  = []
		assoc_track = []
	
		# loop through all the gmcs
		for gmc in gmc_cat:
	
			gmc_coord = SkyCoord(gmc['XCTR_DEG']*u.deg, gmc['YCTR_DEG']*u.deg, frame='fk5', distance=dist*u.Mpc)
			
			# gmc radius		
			gmc_rad_pc = gmc['RAD_PC']
			# convert gmc radius in pc to arcsec
			gmc_rad_rad  = gmc_rad_pc/(dist*1e6)
			gmc_rad_asec = gmc_rad_rad*u.rad.to(u.arcsec)
	
			# get separation between the gmc and all the star clusters
			sep =  gmc_coord.separation(sc_coords).to(u.arcsec)
			
			# find clusters within 1 gmc radius
			w1 = np.where( sep.value <= gmc_rad_asec)[0]
			# find clusters between 1 and 2 gmc radius
			w2 = np.where((sep.value <= 2*gmc_rad_asec) & (sep.value > gmc_rad_asec))[0]
			# find clusters between 2 and 3 gmc radius
			w3 = np.where((sep.value <= 3*gmc_rad_asec) & (sep.value > 2*gmc_rad_asec))[0]
			
			# check if w1 isn't empty
			if len(w1) > 0: 
				# loop through its contents
				for j in range(len(w1)):
					# keep track of the *index* of the gmc currently being used
					gmc_track.append(gmc['CLOUDNUM'] - 1)
					# keep track of the index of the star cluster
					sc_track.append(w1[j])
					# place 1 in assoc_track to signify that this cluster is within 1 gmc radius
					assoc_track.append(1)
			
			# same agian but for between 1 and 2 gmc radius
			if len(w2) > 0:
				for j in range(len(w2)):
					gmc_track.append(gmc['CLOUDNUM'] - 1)
					sc_track.append(w2[j])
					assoc_track.append(2)
			
			# same agian but for between 2 and 3 gmc radius
			if len(w3) > 0: 
				for j in range(len(w3)):
					gmc_track.append(gmc['CLOUDNUM'] - 1)
					sc_track.append(w3[j])
					assoc_track.append(3)
	
		# make the lists numpy arrays
		gmc_track = np.array(gmc_track)
		sc_track = np.array(sc_track)
		assoc_track = np.array(assoc_track)
	
		# star cluster can only be associated with a single gmc so
		# we look for duplicates of the star clusters found associated with the GMCs
		# using the counts given by np.unique
		uniq, c = np.unique(sc_track, return_counts=True)
	
		# duplicate star clusters are where the counts are greater than 1
		wdup = np.where(c > 1)[0]
		
		# loop through the dupilcate star clusters
		for sc_idx in uniq[wdup]:
	
			# get the index within sc_track where the duplicate is found
			w_sc_dup = np.where(sc_track == sc_idx)[0]
			# get the indices of the gmcs (within gmc_cat) of the multiple gmcs associated with the single star cluster
			wgmc = gmc_track[w_sc_dup]
			# find which of these gmcs is the most massive
			wmassmax = np.argmax(gmc_cat['MLUM_MSUN'][wgmc])
	
			# drop the most massive gmc from the list of duplicates
			w_sc_dup = np.delete(w_sc_dup, wmassmax)
	
			# then delete the duplicates from the arrays
			gmc_track = np.delete(gmc_track, w_sc_dup)
			sc_track = np.delete(sc_track, w_sc_dup)
			assoc_track = np.delete(assoc_track, w_sc_dup)
		
		# make a copy of the dataframe
		sc_df2 = sc_df.copy()
	
		# column with the gmc cloudnum of which the star cluster is associated
		# star clusters unassociated with a gmc will be set to nan
		assoc_gmc_cloudnum = np.zeros(len(sc_df2)) + np.nan
		assoc_gmc_cloudnum[sc_track] = gmc_cat['CLOUDNUM'][gmc_track]
	
		# column with the association number 
		# 0 = unassociated, 1 = w/in 1 gmc rad, 2 = b/w 1 and 2 gmc rad, 3 = b/w 2 and 3 gmc rad
		assoc_num = np.zeros(len(sc_df2))
		assoc_num[sc_track] = assoc_track

		# associated gmc mass and radius
		assoc_gmc_mlum = np.zeros(len(sc_df2)) + np.nan
		assoc_gmc_mlum[sc_track] = gmc_cat['MLUM_MSUN'][gmc_track]
		
		assoc_gmc_radius = np.zeros(len(sc_df2)) + np.nan
		assoc_gmc_radius[sc_track] = gmc_cat['RAD3D_PC'][gmc_track]
	
		# insert those new columns
		sc_df2['assoc_gmc_cloudnum']  = assoc_gmc_cloudnum
		sc_df2['assoc_gmc_mlum'] 	  = assoc_gmc_mlum
		sc_df2['assoc_gmc_radius_pc'] = assoc_gmc_radius
		sc_df2['assoc_num'] 		  = assoc_num
	
		# output to csv file
		sc_df2.to_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, filename), index=False)


def make_sc_gmc_assoc_hists(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_class12_assoc_gmc', 
							plot_errorbars=False):
	""" loop through all the galaxies and make the histograms of the cluster ages with their gmc associations

	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):

		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# read in the invidual galaxy dataframe rather than using the mega
		df = pd.read_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, assoc_cat_suffix))

		w0 = df['assoc_num'] == 0
		w1 = df['assoc_num'] == 1
		w2 = df['assoc_num'] == 2
		w3 = df['assoc_num'] == 3

		sc_gmc_assoc_hist(df, filename=data_dir+'%s/%s/%s_sc_gmc_assoc_hist'%(gal_name, run_name, gal_name), errorbars=plot_errorbars)

		# star cluster ages
		age_all = df['age'].to_numpy()
		lage_all = np.log10(age_all)
		age_err_all = df['age_err'].to_numpy()
		lage_err_all = age_err_all/age_all/np.log(10)

		# errors on the median ages 
		med_age_sigma_all = bootstrap_median_error(age_all, age_err_all)
		med_age_sigma1    = bootstrap_median_error(age_all[w1], age_err_all[w1])
		med_age_sigma2    = bootstrap_median_error(age_all[w2], age_err_all[w2])
		med_age_sigma3    = bootstrap_median_error(age_all[w3], age_err_all[w3])
		med_age_sigma0    = bootstrap_median_error(age_all[w0], age_err_all[w0])

		# log the stats for each galaxy
		if i == 0:
			f = open(data_dir + 'sc_gmc_assoc_stats.%s.txt'%run_name, 'w')
			f.write(gal_name + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), med_age_sigma_all))
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), med_age_sigma1))
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), med_age_sigma2))
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), med_age_sigma3))
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), med_age_sigma0))
			f.write('\n')
			f.close()

		else:
			f = open(data_dir + 'sc_gmc_assoc_stats.%s.txt'%run_name, 'a')
			f.write(gal_name + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), med_age_sigma_all))
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), med_age_sigma1))
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), med_age_sigma2))
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), med_age_sigma3))
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), med_age_sigma0))
			f.write('\n')

			f.close()


def generate_mega_df(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_class12_assoc_gmc', output=True):
	""" generates the 'mega df' which is a dataframe which has the star clusters (and all their nearest neighbor GMCs and GMC association info) 
	for all the galaxies in the galaxy_list


	"""
	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']

	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):

		gal_name = gal_id[i]

		# read in the csv
		df = pd.read_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, assoc_cat_suffix))

		# need column for the galaxy name
		gal_name_col = np.array([gal_name for i in range(len(df))])

		df.insert(loc=0, column='gal_name', value=gal_name_col)

		if i == 0:
			# initialize the mega df using the first galaxy
			mega_df = df.copy()
		else:
			mega_df = mega_df.append(df)

	if output:
		mega_df.to_csv(data_dir + 'sc_gmc_assoc_mega.%s.csv'%run_name, index=False)

	return mega_df

def sc_gmc_assoc_hist(df, filename, errorbars=False, **kwargs):
	""" histograms of the star cluster ages based on their association with gmcs

	"""

	w0 = df['assoc_num'] == 0
	w1 = df['assoc_num'] == 1
	w2 = df['assoc_num'] == 2
	w3 = df['assoc_num'] == 3

	age_all = df['age'].to_numpy()
	lage_all = np.log10(age_all)
	age_err_all = df['age_err'].to_numpy()
	lage_err_all = age_err_all/age_all/np.log(10)

	med_age_sigma_all = bootstrap_median_error(age_all, age_err_all)
	med_age_sigma1    = bootstrap_median_error(age_all[w1], age_err_all[w1])
	med_age_sigma2    = bootstrap_median_error(age_all[w2], age_err_all[w2])
	med_age_sigma3    = bootstrap_median_error(age_all[w3], age_err_all[w3])
	med_age_sigma0    = bootstrap_median_error(age_all[w0], age_err_all[w0])

	# gets bins for all clusters so we can have the same bins for all the hists
	hist, bin_edges = np.histogram(lage_all, bins=10)

	fig, ax1 = plt.subplots(1,1, figsize=(6,4))

	ax1.hist(lage_all, bins=bin_edges,     density=True, histtype='step', edgecolor='black',   lw=2.5, label='All SCs (%i)'%(len(df)))
	ax1.hist(lage_all[w1], bins=bin_edges, density=True, histtype='step', edgecolor='#377eb8', lw=2.5, label=r'$\leq 1$ $R_{\rm GMC}$ (%i)'%(len(df[w1])))
	ax1.hist(lage_all[w2], bins=bin_edges, density=True, histtype='step', edgecolor='#ff7f00', lw=2.5, label=r'$1 < R_{\rm GMC} \leq 2 $ (%i)'%(len(df[w2])))
	ax1.hist(lage_all[w3], bins=bin_edges, density=True, histtype='step', edgecolor='#e41a1c', lw=2.5, label=r'$2 < R_{\rm GMC} \leq 3 $ (%i)'%(len(df[w3])))
	ax1.hist(lage_all[w0], bins=bin_edges, density=True, histtype='step', edgecolor='#984ea3', lw=2.5, label=r'Unassociated (%i)'%(len(df[w0])))

	ax1.axvline(np.median(lage_all),     ls='--', color='black',   lw=1.5)
	ax1.axvline(np.median(lage_all[w1]), ls='--', color='#377eb8', lw=1.5)
	ax1.axvline(np.median(lage_all[w2]), ls='--', color='#ff7f00', lw=1.5)
	ax1.axvline(np.median(lage_all[w3]), ls='--', color='#e41a1c', lw=1.5)
	ax1.axvline(np.median(lage_all[w0]), ls='--', color='#984ea3', lw=1.5)

	xmi, xma, ymi, yma = ax1.axis()

	# medians
	ax1.scatter(np.median(lage_all), yma, marker='D', color='black', s=30)
	ax1.scatter(np.median(lage_all[w1]), yma, marker='D', color='#377eb8', s=30)
	ax1.scatter(np.median(lage_all[w2]), yma, marker='D', color='#ff7f00', s=30)
	ax1.scatter(np.median(lage_all[w3]), yma, marker='D', color='#e41a1c', s=30)
	ax1.scatter(np.median(lage_all[w0]), yma, marker='D', color='#984ea3', s=30)

	if errorbars:
		ax1.errorbar(np.median(lage_all), yma, xerr=med_age_sigma_all, ecolor='black', zorder=0)
		ax1.errorbar(np.median(lage_all[w1]), yma, xerr=med_age_sigma1, ecolor='#377eb8', zorder=0)
		ax1.errorbar(np.median(lage_all[w2]), yma, xerr=med_age_sigma2, ecolor='#ff7f00', zorder=0)
		ax1.errorbar(np.median(lage_all[w3]), yma, xerr=med_age_sigma3, ecolor='#e41a1c', zorder=0)
		ax1.errorbar(np.median(lage_all[w0]), yma, xerr=med_age_sigma0, ecolor='#984ea3', zorder=0)


	ax1.set_xlabel('log(Age) [Myr]')
	ax1.set_ylabel('Number')

	ax1.set_xlim(xmi, xma+0.6)

	plt.legend(loc='upper right', fontsize='x-small')


	plt.savefig(filename + '.png', bbox_inches='tight')
	plt.savefig(filename + '.pdf', bbox_inches='tight')

	plt.close()


def all_galaxies_sc_gmc_assoc(galaxy_list, data_dir, run_name, sc_mask_cat_suffix='_cluster_catalog_in_mask_class12', 
							  gmc_cat_suffix='_12m+7m+tp_co21_nativeres_nativenoise_props', sc_class='class12'):
	""" function form of sc_gmc_assoc.py - loops through all the galaxies and find the gmcs associated with each star cluster

	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	sc_cat_suffix 		str 			suffix for the filename of the star cluster catalog; defaults to the class 1 and 2 catalog;
										included here to make it easier to do the runs with class 1,2,3 and the associations
	gmc_cat_suffix  	str 			suffix for the filename of the gmc catalog; defaults to the latest gmc catalog
	sc_class 			str 			which class of clusters to make the catalogs for; class12 or class123

	"""
	# do the star cluster - GMC assocation and output to a csv dataframe
	generate_sc_gmc_assoc_df(galaxy_list, data_dir, run_name, sc_mask_cat_suffix=sc_mask_cat_suffix, gmc_cat_suffix=gmc_cat_suffix, 
							 filename='_cluster_catalog_in_mask_%s_assoc_gmc'%sc_class)

	# create the histograms of the cluster ages split by how they are associated with a gmc 
	make_sc_gmc_assoc_hists(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_%s_assoc_gmc'%sc_class, 
							plot_errorbars=False)

	# generate the 'mega df' 
	mega_df = generate_mega_df(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_%s_assoc_gmc'%sc_class, output=True)

	# grab indices of the clustes for each association number
	w0 = mega_df['assoc_num'] == 0 	# unassociated
	w1 = mega_df['assoc_num'] == 1 	# w/in 1 gmc radius
	w2 = mega_df['assoc_num'] == 2 	# b/w 1 and 2 gmc radii
	w3 = mega_df['assoc_num'] == 3 	# b/w 2 and 3 gmc radii

	# bootstrap estimate of the uncertainty on the median ages
	med_age_sigma_all = bootstrap_median_error(mega_df['age'],     mega_df['age_err'] )
	med_age_sigma0    = bootstrap_median_error(mega_df['age'][w0], mega_df['age_err'][w0] )
	med_age_sigma1    = bootstrap_median_error(mega_df['age'][w1], mega_df['age_err'][w1] )
	med_age_sigma2    = bootstrap_median_error(mega_df['age'][w2], mega_df['age_err'][w2] )
	med_age_sigma3    = bootstrap_median_error(mega_df['age'][w3], mega_df['age_err'][w3] )

	# print out stats
	print('All star clusters median age: %.2f +/- %.2f Myr'%(np.median(mega_df.age), med_age_sigma_all))
	print('1 Rgmc median age: %.2f +/- %.2f'%(mega_df.loc[w1].median()['age'], med_age_sigma1))
	print('1 - 2 Rgmc median age: %.2f +/- %.2f'%(mega_df.loc[w2].median()['age'], med_age_sigma2))
	print('2 - 3 Rgmc median age: %.2f +/- %.2f'%(mega_df.loc[w3].median()['age'], med_age_sigma3))
	print('Unassoc median age: %.2f+/- %.2f'%(mega_df.loc[w0].median()['age'], med_age_sigma0))

	# append the stats across all galaxies 
	f = open(data_dir + 'sc_gmc_assoc_stats.%s.txt'%run_name, 'a')
	f.write('All Galaxies\n')
	f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(mega_df.age), med_age_sigma_all))
	f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(mega_df[w1].age), med_age_sigma1))
	f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(mega_df[w2].age), med_age_sigma2))
	f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(mega_df[w3].age), med_age_sigma3))
	f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(mega_df[w0].age), med_age_sigma0))
	f.write('\n')
	f.close()

	"""	now do things by environmental mask locations 

	simple environmental masks cheatsheet
	1 = center (small bulge, nuclear ring & disk)
	2 = bar (excluding bar ends)
	3 = bar ends (overlap of bar and spiral)
	4 =​ interbar (R_gal < R_bar, but outside bar footprint)
	5 = ​spiral arms inside interbar (R_gal < R_bar)
	6 = ​spiral arms (R_gal > R_bar)
	7 =​ interarm (only the R_gal spanned by spiral arms, and R_gal > R_bar)
	8 = ​outer disc (R_gal > spiral arm ends, only for galaxies with identified spirals)
	9 = ​disc (R_gal > R_bar) where no spiral arms were identified (e.g. flocculent spirals)

	simplified further
	1 =​ center
	2 + 3 =​ bar
	4 + 7 + 8 =​ interarm
	5 + 6 =​ spiral arms
	9 =​ disc in galaxies without spirals
	
	"""
	# get indices for the clusters of each enviro - need np.where so we can get mulitple conditions and can us iloc later
	wcenter			= np.where(mega_df['env_mask_val'] == 1)
	wbar_idx		= np.where((mega_df['env_mask_val'] == 2) | (mega_df['env_mask_val'] == 3) ) 
	winterarm_idx	= np.where((mega_df['env_mask_val'] == 4) | (mega_df['env_mask_val'] == 7) | (mega_df['env_mask_val'] == 8))
	wspiral_idx		= np.where((mega_df['env_mask_val'] == 5) | (mega_df['env_mask_val'] == 6))
	wdisk			= np.where(mega_df['env_mask_val'] == 9)

	# list with all the enviro indices
	wall = [wcenter, wbar_idx[0], winterarm_idx[0], wspiral_idx[0], wdisk]
	# list of the enviro names
	names = ['center', 'bar', 'interarm', 'spiralarm', 'disk']

	# loop through to each enviro
	for i in range(len(wall)):

		# make a temp dataframe with just the clusters of the current enviro
		df = mega_df.iloc[wall[i]]

		# make histogram of cluster ages split by association number
		sc_gmc_assoc_hist(df, filename=data_dir+'sc_gmc_assoc_hist_%s.%s'%(names[i], run_name))

		# star cluster ages and errors
		age_all = df['age'].to_numpy()
		lage_all = np.log10(age_all)
		age_err_all = df['age_err'].to_numpy()
		lage_err_all = age_err_all/age_all/np.log(10)

		# indices for each association number
		w0 = df['assoc_num'] == 0
		w1 = df['assoc_num'] == 1
		w2 = df['assoc_num'] == 2
		w3 = df['assoc_num'] == 3

		# bootstrap errors on the median ages
		med_age_sigma_all = bootstrap_median_error(age_all, age_err_all)
		med_age_sigma1    = bootstrap_median_error(age_all[w1], age_err_all[w1])
		med_age_sigma2    = bootstrap_median_error(age_all[w2], age_err_all[w2])
		med_age_sigma3    = bootstrap_median_error(age_all[w3], age_err_all[w3])
		med_age_sigma0    = bootstrap_median_error(age_all[w0], age_err_all[w0])

		# log the stats for each env
		if i == 0:
			f = open(data_dir + 'sc_gmc_assoc_stats_env.%s.txt'%run_name, 'w')
			f.write(names[i] + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), med_age_sigma_all) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), med_age_sigma1 ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), med_age_sigma2 ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), med_age_sigma3 ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), med_age_sigma0 ) )
			f.write('\n')
			f.close()
		else:
			f = open(data_dir + 'sc_gmc_assoc_stats_env.%s.txt'%run_name, 'a')
			f.write(names[i] + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), med_age_sigma_all) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), med_age_sigma1 ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), med_age_sigma2 ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), med_age_sigma3 ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), med_age_sigma0 ) )
			f.write('\n')
			f.close()

def auto_corr(df, min_bin=1.1e-5, nbins=10, nbootstraps=50, method='landy-szalay', rseed=222, gmc=False):
	""" function to calculate the auto-correlation for the given dataframe
	uses the astroML function bootstrap_two_point_angular 

	Inputs:
	df 			pandas DataFrame 	dataframe which holds the objects to do the correlation function for
	min_bin 	float 				the angular location of the first/minimum bin
	nbins 		int 				the number of radial bins over which to do the correlation
	nbootstraps int 				number of bootstraps to perform for the error estimation; default is 50
	method 		str 				estimator method to use for correlation function; landy-szalay or standard; default is landy-szalay
	rseed 		int 				the seed value which gets used for the numpy.random
	gmc 		bool 				set to true if the df used is the GMC catalog since it has different keywords for ra,dec

	Outputs:
	results
		results[0] == bins 			list of the bin edges; len(bins) == nbins + 1
		results[1] == corr 			list of the correlation values for each bin
		results[2] == corr_err 		list of the bootstrap estimated errors on the correlation values
		results[3] == bootstraps 	list of lists of bootstrapped correlation values in each bin; len(bootstraps) == nbootstraps

	""" 
	np.random.seed(rseed)

	bins = 10 ** np.linspace(np.log10(min_bin), np.log10(0.1), nbins+1)
	results = [bins]

	if gmc:
		results += bootstrap_two_point_angular(df['XCTR_DEG'], df['YCTR_DEG'], bins=bins, method=method, Nbootstraps=nbootstraps)

	else:
		results += bootstrap_two_point_angular(df['ra'], df['dec'], bins=bins, method=method, Nbootstraps=nbootstraps)

	return results

def powerlaw_func(theta, Aw, alpha):
	""" a powerlaw function of the form 
	f(theta) = Aw * theta^alpha

	"""
	return Aw * theta**alpha

def tpcf(df, dist, **kwargs):
	""" runs the bootstrap two point corrrelation function and the power law fit

	Inputs:
	df 		pandas DataFrame 	dataframe which holds the objects to do the correlation function for
	dist 	float 				distance to galaxy in Mpc
	kwargs 	dictionary			keyword arguments to pass on to the auto_corr function

	Outputs:
	bins_centers_pc 	list 	center positions of the bins in parsecs 
	corr 				list 	correlation values for each bin; 1 + omega(theta)
	corr_err 			list   	bootstrap estimated errors on the correlation values
	power_law_fits 		list 	the best-fit for powerlaws; [A_w (deg), error, A_w (pc), error, alpha, error ]


	"""
	# perform the auto-correlation
	bins, corr, corr_err, bootstraps = auto_corr(df, **kwargs)
	
	# find bin centers [degrees]
	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# bin centers as in pc
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad) 

	# add 1 so the correlation is 1 + omega(theta)
	corr = corr + 1

	# need to drop nans for the power law fitting
	wnnan = np.where(np.isnan(corr)==False)

	# power-law fit
	popt_ang, pcov = curve_fit(powerlaw_func, bin_centers[wnnan], corr[wnnan])
	perr_ang 	   = np.sqrt(np.diag(pcov))
	popt_pc, pcov  = curve_fit(powerlaw_func, bin_centers_pc[wnnan], corr[wnnan])
	perr_pc 	   = np.sqrt(np.diag(pcov))

	# sometimes the error doesn't converge so replace those with 0 (instead of inf)
	winf = np.where(np.isinf(perr_ang))[0]
	if len(winf) > 0:
		perr_ang[winf] = 0
		perr_pc[winf] = 0

	return bin_centers_pc, corr, corr_err, [popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]]

def all_galaxies_tpcf(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_class12_assoc_gmc', sc_class='class12', nbins=10 ):
	""" function form of tpcf.py - loop through all the galaxies and do the two-point correlation function analysis
	
	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	assoc_cat_suffix 	str 			suffix of the filename for the csv which holds the star cluster - gmc association dataframe
	sc_class 			str 			which class of clusters to make the catalogs for; class12 or class123
	nbins 				int; list		the number of radial bins over which to do the correlation; if a list, it'll loop through all the given nbsins


	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']

	for i in range(len(galaxy_list)):

		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# read in the star cluster cat in the hst-alma footprint overlap mask
		sc_df = pd.read_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, assoc_cat_suffix))

		# read in the gmc cat in the hst-alma footprint overlap mask
		gmc_cat = fits.open(data_dir + '%s/%s/%s_gmc_cat_masked.fits'%(gal_name, run_name, gal_name))[1].data
		gmc_df = Table(gmc_cat).to_pandas()

		# check if nbins is a list or int; if int, make it a list of len 1
		if type(nbins) == int:
			nbins = [nbins]
		
		# loop through the nbins in the list
		for j in range(len(nbins)):

			# two-point correlation function on the all the star clusters
			bin_centers_pc_all, corr_all, corr_err_all, pl_fit_all = tpcf(sc_df, dist, nbins=nbins[j])

			# now for clusters <= 10 Myr
			wleq10 = sc_df['age'] <= 10
			bin_centers_pc_young, corr_young, corr_err_young, pl_fit_young = tpcf(sc_df.loc[wleq10], dist, nbins=nbins[j])

			# now for clusters > 10 Myr
			w10 = sc_df['age'] > 10
			bin_centers_pc_old, corr_old, corr_err_old, pl_fit_old = tpcf(sc_df.loc[w10], dist, nbins=nbins[j])

			# now gmcs
			bin_centers_pc_gmc, corr_gmc, corr_err_gmc, pl_fit_gmc = tpcf(gmc_df, dist, nbins=nbins[j], min_bin=3e-4, gmc=True)

			# write out the power-law best fit parameters for all, young, old
			f = open(data_dir + '%s/%s/%s_tpcf_fits.nbins%02d.dat'%(gal_name, run_name, gal_name, nbins[j]), 'w')

			f.write('{:<6}  '.format('# bin'))
			f.write('{:<6}  '.format('Aw_deg'))
			f.write('{:<5}  '.format('error'))
			f.write('{:<6}  '.format('Aw_pc'))
			f.write('{:<6}  '.format('error'))
			f.write('{:<6}  '.format('alpha'))
			f.write('{:<5}  '.format('error'))
			f.write('\n')

			f.write('{:<6}  '.format('all'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_all[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_all[5])))
			f.write('\n')

			f.write('{:<6}  '.format('<= 10'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_young[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_young[5])))
			f.write('\n')

			f.write('{:<6}  '.format('> 10'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_old[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_old[5])))
			f.write('\n')

			f.write('{:<6}  '.format('gmc'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_gmc[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_gmc[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_gmc[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_gmc[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_gmc[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_gmc[5])))
			f.close()

			# create figure

			fig, ax = plt.subplots(1,1, figsize=(5,5))
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xlabel(r'$r$ [pc]')
			ax.set_ylabel(r'$1 + \omega(\theta)$')

			# all clusters
			ax.errorbar(bin_centers_pc_all, corr_all, yerr=corr_err_all, fmt='k-o', ecolor='black', markersize=5, lw=1.5, 
					label=r'All SCs $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_all[4], pl_fit_all[5], len(sc_df)))
			# clusters <= 10 Myr
			ax.errorbar(bin_centers_pc_young, corr_young, yerr=corr_err_young, fmt='-o', color='#377eb8', ecolor='#377eb8', markersize=5, lw=1.5, 
					label=r'$\leq 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_young[4], pl_fit_young[5], len(sc_df.loc[wleq10])))
			# clusters > 10 Myr
			ax.errorbar(bin_centers_pc_old, corr_old, yerr=corr_err_old, fmt='-o', color='#e41a1c', ecolor='#e41a1c', markersize=5, lw=1.5, 
					label=r'$> 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_old[4], pl_fit_old[5], len(sc_df.loc[w10])))
			# gmcs
			ax.errorbar(bin_centers_pc_gmc, corr_gmc, yerr=corr_err_gmc, fmt='-o', color='#E68310', ecolor='#E68310', markersize=5, lw=1.5, 
					label=r'GMCs $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_gmc[4], pl_fit_gmc[5], len(gmc_df)))
			# plot vertical line at mean GMC radius
			ax.axvline(gmc_df.mean()['RAD3D_PC'], lw=1.1, c='#999999', zorder=0)

			plt.legend(loc='upper right', fontsize='x-small')
			plt.savefig(data_dir + '%s/%s/%s_tpcf.nbins%02d.png'%(gal_name, run_name, gal_name, nbins[j]), bbox_inches='tight')
			plt.savefig(data_dir + '%s/%s/%s_tpcf.nbins%02d.pdf'%(gal_name, run_name, gal_name, nbins[j]), bbox_inches='tight')
			plt.close()

			logbins_sc  = np.log10(bin_centers_pc_all)
			logbins_gmc = np.log10(bin_centers_pc_gmc)


			# write out the bin centers and correlation values
			f = open(data_dir + '%s/%s/%s_tpcf.nbins%02d.dat'%(gal_name, run_name, gal_name, nbins[j]), 'w')
			f.write('# two-point correlation function values (1 + omega(theta)); bin centers in are given in log(pc)\n')

			f.write('{:<8}  '.format('nbins%02d'%nbins[j]))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(logbins_sc[k])))
			
			f.write('\n')
			f.write('{:<8}  '.format('corr_all'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_all[k])))

			f.write('\n')
			f.write('{:<8}  '.format('corr_yng'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_young[k])))

			f.write('\n')
			f.write('{:<8}  '.format('corr_old'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_old[k])))
	
			f.write('\n')
			f.write('{:<8}  '.format('bins_gmc'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(logbins_gmc[k])))
			
			f.write('\n')
			f.write('{:<8}  '.format('corr_gmc'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_gmc[k])))

			f.close()

def cross_correlate(sc_df, gmc_df, rand_sc_df, rand_gmc_df, dist, nbins=16, min_bin=1.1e-5):
	""" calculate the cross-correlation between the given star clusters and GMC with the Landy & Szalay (1993) estimator
	makes use of astropy search_around_sky which employs kdtrees to find separations b/w objects 
	follows the same methodology as astroML two point correlation function (https://www.astroml.org/user_guide/correlation_functions.html)
	but modified to just use search_around_sky (so i didn't have to write the kdtree stuff since it's implemented already in search_around_sky)
	and follows the math for doing cross-correlation rather than auto-correlation

	Inputs:
	sc_df 		pandas DataFrame 	dataframe which contains the star cluster catalog; needs to include the cluster ra, dec
	gmc_df 		pandas DataFrame 	dataframe which contains the gmc catalog; needs to include the gmc cluster ra, dec
	rand_sc_df 	pandas DataFrame 	dataframe which contains the random star cluster catalog over which the correlation will be compared to; needs the random cluster ra, dec
	rand_gmc_df pandas DataFrame 	dataframe which contains the random gmc catalog over which the correlation will be compared to; needs the random gmc ra, dec
	dist 		float 				distance to the galaxy in Mpc; needed for calculating the separations b/w clusters and gmcs in pc
	nbins 		int 				the number of radial bins over which to do the correlation
	min_bin 	float 				the angular location of the first/minimum bin

	Outputs:
	bins 		array 				array of the angular bin edges (degrees)
	corr 		array 				correlation values --> 1 + omega(theta)

	"""

	bins = 10 ** np.linspace(np.log10(min_bin), np.log10(0.1), nbins+1)

	# total numbers in the two catalogs
	N_Dsc  = len(sc_df)
	N_Dgmc = len(gmc_df)

	N_Rsc  = len(rand_sc_df)
	N_Rgmc = len(rand_gmc_df)

	# make sky coords for all the cats
	sc_coords  = SkyCoord(ra=sc_df['ra']*u.deg, dec=sc_df['dec']*u.deg, frame='icrs', distance=dist*u.Mpc)
	gmc_coords = SkyCoord(ra=gmc_df['XCTR_DEG']*u.deg, dec=gmc_df['YCTR_DEG']*u.deg, frame='fk5', distance=dist*u.Mpc)
	
	rand_sc_coords  = SkyCoord(ra=rand_sc_df['ra']*u.deg, dec=rand_sc_df['dec']*u.deg, frame='icrs', distance=dist*u.Mpc)
	rand_gmc_coords = SkyCoord(ra=rand_gmc_df['ra']*u.deg, dec=rand_gmc_df['dec']*u.deg, frame='icrs', distance=dist*u.Mpc)

	# get separations for all pairs 
	DscDgmc_sas = search_around_sky(sc_coords, gmc_coords, seplimit=5*u.deg)
	DscRgmc_sas = search_around_sky(sc_coords, rand_gmc_coords, seplimit=5*u.deg)
	RscDgmc_sas = search_around_sky(rand_sc_coords, gmc_coords, seplimit=5*u.deg)
	RscRgmc_sas = search_around_sky(rand_sc_coords, rand_gmc_coords, seplimit=5*u.deg)

	DscDgmc_sep = DscDgmc_sas[2]
	DscRgmc_sep = DscRgmc_sas[2]
	RscDgmc_sep = RscDgmc_sas[2]
	RscRgmc_sep = RscRgmc_sas[2]

	# loop through bins [drop first bin/min_bin]
	lbins = bins[1:] * u.deg

	corr = []
	for j in range(len(lbins)):

		if j == 0:
			wDDr = np.where((DscDgmc_sep >= min_bin*u.deg) & (DscDgmc_sep < lbins[j]) )[0]
			wDRr = np.where((DscRgmc_sep >= min_bin*u.deg) & (DscRgmc_sep < lbins[j]) )[0]
			wRDr = np.where((RscDgmc_sep >= min_bin*u.deg) & (RscDgmc_sep < lbins[j]) )[0]
			wRRr = np.where((RscRgmc_sep >= min_bin*u.deg) & (RscRgmc_sep < lbins[j]) )[0]

		else:
			wDDr = np.where((DscDgmc_sep >= lbins[j-1]) & (DscDgmc_sep < lbins[j]) )[0]
			wDRr = np.where((DscRgmc_sep >= lbins[j-1]) & (DscRgmc_sep < lbins[j]) )[0]
			wRDr = np.where((RscDgmc_sep >= lbins[j-1]) & (RscDgmc_sep < lbins[j]) )[0]
			wRRr = np.where((RscRgmc_sep >= lbins[j-1]) & (RscRgmc_sep < lbins[j]) )[0]

		# check for empty RR, and replace with nan if needed
		if len(wRRr) == 0:
			corr.append(np.nan)
		else:
			factor1 = (N_Rsc * N_Rgmc * len(wDDr))/(N_Dsc * N_Dgmc * len(wRRr))
			factor2 = (N_Rsc * len(wDRr))/(N_Dsc * len(wRRr))
			factor3 = (N_Rgmc * len(wRDr))/(N_Dgmc * len(wRRr))
	
			xi_r = factor1 - factor2 - factor3 + 1
	
			corr.append(xi_r)


	corr = np.array(corr) + 1

	return bins, corr

def cross_correlate_bootstrap(sc_df, gmc_df, wcs_hst, xmax_hst, ymax_hst, wcs_alma, xmax_alma, ymax_alma, mask, dist, nbootstraps=100, **kwargs):
	""" run the cross-correlation with a non-parametric bootstrap estimation of the errors on the correlation values
	calls on the cross_correlate function above but runs through the given number of bootstraps and returns an estimation of the error

	bootstrap estimation is done by randomly resampling the star cluster and gmc catalogs and then recalculating the correlation over
	a new random cluster and gmc catalog
	the correlation value (1 + omega(theta)) is saved for each bin for each bootstrap run
	and error is estimated as the standard deviation (with ddof = 1) of the correlation values in each bin

	Inputs:
	sc_df 			pandas DataFrame 		dataframe which contains the star cluster catalog; needs to include the cluster ra, dec
	gmc_df 			pandas DataFrame 		dataframe which contains the gmc catalog; needs to include the gmc cluster ra, dec
	wcs_hst 		astropy.wcs.wcs.WCS 	astropy wcs object of the hst image
	xmax_hst 		int 					the maximum pixel location in the x-direction of the hst image; i.e., if hst image is 13000x14000 pixels, xmax_hst = 13000
	ymax_hst 		int 					the maximum pixel location in the y-direction of the hst image; i.e., if hst image is 13000x14000 pixels, ymax_hst = 14000
	wcs_alma 		astropy.wcs.wcs.WCS 	astropy wcs object of the alma image
	xmax_alma 		int 					the maximum pixel location in the x-direction of the alma image; i.e., if alma image is 1599x1598 pixels, xmax_alma = 1599
	ymax_alma		int 					the maximum pixel location in the y-direction of the alma image; i.e., if alma image is 1599x1598 pixels, ymax_alma = 1598
	mask 			array 					the data array of the hst-alma overlap mask image
	dist 			float 					distance to the galaxy in Mpc; needed for calculating the separations b/w clusters and gmcs in pc
	nbootstraps 	int 					number of bootstraps to run through
	**kwargs 		dict 					keyword arguments to pass on to the cross_correlate funciton; really just nbins and min_bin

	Outputs:
	bins_centers_pc	array 					centers of the bins in parsecs
	corr 			array 					correlation values in each bin --> 1 + omega(theta)
	corr_err 		array 					1 sigma error on the correlation values; if correlation is nan, error will be 0
	power_law_fits 	list 				 	the best-fit for powerlaws; [A_w (deg), error, A_w (pc), error, alpha, error ]
	
	"""

	bootstraps = []

	for i in range(nbootstraps):

		# generate random star cluster and gmc catalogs
		rand_sc_df  = generate_random_sc(sc_df, xmax_hst, ymax_hst, wcs_hst, mask)
		rand_gmc_df = generate_random_gmc(gmc_df, xmax_alma, ymax_alma, wcs_alma, wcs_hst, mask)

		# random resampling of the data (sc and gmcs) unless its the first time through
		# then just run the original data
		if i > 0:
		
			rand_sc_ind  = np.random.randint(0, len(sc_df), len(sc_df) )
			rand_gmc_ind = np.random.randint(0, len(gmc_df), len(gmc_df) )

			sc_boot  = sc_df.iloc[rand_sc_ind]
			gmc_boot = gmc_df.iloc[rand_gmc_ind]

		else:

			sc_boot = sc_df
			gmc_boot = gmc_df

		# run cross-correlation; bins won't change through out the bootstraps so it's ok to overwrite but we do want to return it later
		bins, corr = cross_correlate(sc_boot, gmc_boot, rand_sc_df, rand_gmc_df, dist, **kwargs)

		# save the correlation array from each bootstrap
		bootstraps.append(corr)

	# since the first one of the bootstraps was the original data, this is the correlation to return
	corr = bootstraps[0]

	# since there are nans in the correlation results, we mask them out and comput the standard deviations in each bin 
	# delta degree of freedom = 1 because the bootstraps are computed from a random sample of the population
	# bins with nan correlation will be given 0 for the error on the correlation
	corr_err = np.asarray(np.ma.masked_invalid(bootstraps).std(0, ddof=1))

	# get centers of bins [degrees]
	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# bin centers in parsec
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad)

	# will need to drop nans for power law fitting
	wnnan = np.where(np.isnan(corr) == False)

	try:
		# power-law fit
		popt_ang, pcov = curve_fit(powerlaw_func, bin_centers[wnnan], corr[wnnan])
		perr_ang 	   = np.sqrt(np.diag(pcov))
		popt_pc, pcov  = curve_fit(powerlaw_func, bin_centers_pc[wnnan], corr[wnnan])
		perr_pc 	   = np.sqrt(np.diag(pcov))
	except:
		print('\ncross correlation power law fit failed for nbins = %i'%(len(bin_centers)))

		popt_ang	= [0,0]
		perr_ang	= [0,0]
		popt_pc		= [0,0]
		perr_pc		= [0,0]
		popt_ang 	= [0,0]
		perr_ang	= [0,0]

	# sometimes the error doesn't converge so replace those with 0 (instead of inf)
	winf = np.where(np.isinf(perr_ang))[0]
	if len(winf) > 0:
		perr_ang[winf] = 0
		perr_pc[winf] = 0

	return bin_centers_pc, corr, corr_err, [popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]]

def generate_random_sc(sc_df, xmax, ymax, wcs_hst, mask, rseed=222):
	""" generates a catalog of randomly placed 'star clusters'

	Inputs:
	sc_df 		pandas DataFrame 		dataframe which contains the star cluster catalog; really just to get the number of actual star clusters
	xmax 	 	int 					the maximum pixel location in the x-direction of the hst image; i.e., if hst image is 13000x14000 pixels, xmax_hst = 13000
	ymax 	 	int 					the maximum pixel location in the y-direction of the hst image; i.e., if hst image is 13000x14000 pixels, ymax_hst = 14000
	wcs_hst 	astropy.wcs.wcs.WCS 	astropy wcs object of the hst image
	mask 		array 					the data array of the hst-alma overlap mask image
	rseed 		int 					the seed value which gets used for the numpy.random

	Output:
	rand_sc_df 	pandas DataFrame 		dataframe of the random star cluster catalog; includes x, y, ra, dec

	"""
	
	np.random.seed(rseed)

	# generate random uniform distribution of star clusters over the entire HST image
	rand_x_full = np.random.uniform(0, xmax, size=4*len(sc_df))
	rand_y_full = np.random.uniform(0, ymax, size=4*len(sc_df))

	# convert x,y values to integers b/c np.random.uniform can only produce floats
	rand_x_full_int = np.array([int(np.round(x)) for x in rand_x_full ])
	rand_y_full_int = np.array([int(np.round(y)) for y in rand_y_full ])

	# limit to those in the hst-alma overlap mask
	rand_in_mask = np.array([True if mask[y,x] == 1 else False for y,x in zip(rand_y_full_int, rand_x_full_int)])

	# convert to ra, dec using the hst image wcs info
	rand_ra, rand_dec = wcs_hst.wcs_pix2world(rand_x_full[rand_in_mask], rand_y_full[rand_in_mask], 1.0)

	# create dictionary with the random clusters' x,y,ra,dec positions
	rand_sc_data = {'x': rand_x_full[rand_in_mask], 'y': rand_y_full[rand_in_mask], 'ra': rand_ra, 'dec': rand_dec }

	# create a dataframe from the dictionary
	rand_sc_df = pd.DataFrame(rand_sc_data)

	return rand_sc_df

def generate_random_gmc(gmc_df, xmax, ymax, wcs_alma, wcs_hst, mask, rseed=222):
	""" generates a catalog of randomly placed 'gmcs'

	Inputs:
	gmc_df 		pandas DataFrame 		dataframe which contains the gmc catalog; really just to get the number of actual gmcs
	xmax 	 	int 					the maximum pixel location in the x-direction of the alma image; i.e., if alma image is 1599x1598 pixels, xmax_alma = 1599
	ymax 	 	int 					the maximum pixel location in the y-direction of the alma image; i.e., if alma image is 1599x1598 pixels, ymax_alma = 1598
	wcs_alma 	astropy.wcs.wcs.WCS 	astropy wcs object of the alma image
	wcs_hst 	astropy.wcs.wcs.WCS 	astropy wcs object of the hst image
	mask 		array 					the data array of the hst-alma overlap mask image
	rseed 		int 					the seed value which gets used for the numpy.random

	Output:
	rand_gmc_df pandas DataFrame 		dataframe of the random gmc catalog; includes x, y, ra, dec

	"""

	# generate random distribution of gmcs over the entire ALMA image
	rand_x_full = np.random.uniform(0, xmax, size=2*len(gmc_df))
	rand_y_full = np.random.uniform(0, ymax, size=2*len(gmc_df))

	# convert pixels to ra, dec
	rand_ra_full, rand_dec_full = wcs_alma.wcs_pix2world(rand_x_full, rand_y_full, 1.0)

	# conver ra, dec to HST pixels
	rand_x_full, rand_y_full = wcs_hst.wcs_world2pix(rand_ra_full, rand_dec_full, 1.0)

	# convert x,y values to integers b/c np.random.uniform can only produce floats
	rand_x_full_int = np.array([int(np.round(x)) for x in rand_x_full ])
	rand_y_full_int = np.array([int(np.round(y)) for y in rand_y_full ])

	# limit to those in the hst-alma overlap mask
	rand_in_mask = np.array([True if mask[y,x] == 1 else False for y,x in zip(rand_y_full_int, rand_x_full_int)])
	
	# convert to ra, dec using the hst image wcs info
	rand_ra, rand_dec = wcs_hst.wcs_pix2world(rand_x_full[rand_in_mask], rand_y_full[rand_in_mask], 1.0)
	
	# create dictionary with the random gmcs' x,y,ra,dec positions
	rand_gmc_data = {'x': rand_x_full[rand_in_mask], 'y': rand_y_full[rand_in_mask], 'ra': rand_ra, 'dec': rand_dec }
	
	# create a dataframe from the dictionary
	rand_gmc_df = pd.DataFrame(rand_gmc_data)

	return rand_gmc_df


def all_galaxies_crosscf(galaxy_list, data_dir, run_name, assoc_cat_suffix='_cluster_catalog_in_mask_class12_assoc_gmc', sc_class='class12', nbins=10 ):
	""" function form of cf.py - loop through all the galaxies and do the cross correlation function analysis
	
	Inputs:
	galaxy_list			astropy Table 	table that holds the list of galaxies to perform the analysis on
	data_dir 			str 			path to the data directory; e.g., /cherokee1/turner/phangs/cf/data/
	run_name 			str 			name of the run/test; e.g., run01
	assoc_cat_suffix 	str 			suffix of the filename for the csv which holds the star cluster - gmc association dataframe
	sc_class 			str 			which class of clusters to make the catalogs for; class12 or class123
	nbins 				int; list		the number of radial bins over which to do the correlation; if a list, it'll loop through all the given nbsins


	"""

	gal_id 		= galaxy_list['id']
	gal_dist 	= galaxy_list['dist']

	for i in range(len(galaxy_list)):

		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# read in the star cluster cat in the hst-alma footprint overlap mask
		sc_df = pd.read_csv(data_dir + '%s/%s/%s%s.csv'%(gal_name, run_name, gal_name, assoc_cat_suffix))

		# read in the gmc cat in the hst-alma footprint overlap mask
		gmc_cat = fits.open(data_dir + '%s/%s/%s_gmc_cat_masked.fits'%(gal_name, run_name, gal_name))[1].data
		gmc_df = Table(gmc_cat).to_pandas()

		# read in the HST-ALMA overlap mask (uses same wcs info as the HST image)
		mask_hdu = fits.open(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name))
		mask = mask_hdu[0].data
		mask_header = mask_hdu[0].header
		wcs_hst = WCS(mask_header, fobj=mask_hdu)
		# max pixel numbers
		ymax_hst, xmax_hst = np.shape(mask)

		# read in an ALMA image
		alma_hdulist = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))
		alma_header = alma_hdulist[0].header
		alma_data = alma_hdulist[0].data
		wcs_alma = WCS(alma_header, fobj=alma_hdulist, naxis=2)
		ymax_alma, xmax_alma = np.shape(alma_data)

		# check if nbins is a list or int; if int, make it a list of len 1
		if type(nbins) == int:
			nbins = [nbins]

		for j in range(len(nbins)):

			# cross-correlation function using all the star clusters
			bin_centers_pc_all, corr_all, corr_err_all, pl_fit_all = cross_correlate_bootstrap(sc_df, gmc_df, wcs_hst, xmax_hst, ymax_hst, 
																							  wcs_alma, xmax_alma, ymax_alma, mask, dist, 
																							  nbootstraps=50, nbins=nbins[j])
			
			# now with clusters <= 10 Myr
			wleq10 = sc_df['age'] <= 10
			bin_centers_pc_young, corr_young, corr_err_young, pl_fit_young = cross_correlate_bootstrap(sc_df.loc[wleq10], gmc_df, wcs_hst, xmax_hst, ymax_hst, 
																									   wcs_alma, xmax_alma, ymax_alma, mask, dist, 
																									   nbootstraps=50, nbins=nbins[j])

			# now with clusters > 10 Myr
			w10 = sc_df['age'] > 10
			bin_centers_pc_old, corr_old, corr_err_old, pl_fit_old = cross_correlate_bootstrap(sc_df.loc[w10], gmc_df, wcs_hst, xmax_hst, ymax_hst, 
																							  wcs_alma, xmax_alma, ymax_alma, mask, dist, 
																							  nbootstraps=50, nbins=nbins[j])

			# create figure
			fig, ax = plt.subplots(1,1, figsize=(5,5))
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xlabel(r'$r$ [pc]')
			ax.set_ylabel(r'$1 + \omega(\theta)$')

			# all clusters
			ax.errorbar(bin_centers_pc_all, corr_all, yerr=corr_err_all, fmt='k-o', ecolor='black', markersize=5, lw=1.5, 
					label=r'All SCs $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_all[4], pl_fit_all[5], len(sc_df)))
			# clusters <= 10 Myr
			ax.errorbar(bin_centers_pc_young, corr_young, yerr=corr_err_young, fmt='-o', color='#377eb8', ecolor='#377eb8', markersize=5, lw=1.5, 
					label=r'$\leq 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_young[4], pl_fit_young[5], len(sc_df.loc[wleq10])))
			# clusters > 10 Myr
			ax.errorbar(bin_centers_pc_old, corr_old, yerr=corr_err_old, fmt='-o', color='#e41a1c', ecolor='#e41a1c', markersize=5, lw=1.5, 
					label=r'$> 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(pl_fit_old[4], pl_fit_old[5], len(sc_df.loc[w10])))
			# plot vertical line at mean GMC radius
			ax.axvline(gmc_df.mean()['RAD3D_PC'], lw=1.1, c='#999999', zorder=0)

			plt.legend(loc='upper right', fontsize='x-small')
			plt.savefig(data_dir + '%s/%s/%s_crosscf.nbins%02d.png'%(gal_name, run_name, gal_name, nbins[j]), bbox_inches='tight')
			plt.savefig(data_dir + '%s/%s/%s_crosscf.nbins%02d.pdf'%(gal_name, run_name, gal_name, nbins[j]), bbox_inches='tight')
			plt.close()


			# write out the power law best fits for all, young, old
			f = open(data_dir + '%s/%s/%s_crosscf_fits.nbins%02d.dat'%(gal_name, run_name, gal_name, nbins[j]), 'w')

			f.write('{:<6}  '.format('# bin'))
			f.write('{:<6}  '.format('Aw_deg'))
			f.write('{:<5}  '.format('error'))
			f.write('{:<6}  '.format('Aw_pc'))
			f.write('{:<6}  '.format('error'))
			f.write('{:<6}  '.format('alpha'))
			f.write('{:<5}  '.format('error'))
			f.write('\n')

			f.write('{:<6}  '.format('all'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_all[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_all[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_all[5])))
			f.write('\n')

			f.write('{:<6}  '.format('<= 10'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_young[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_young[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_young[5])))
			f.write('\n')

			f.write('{:<6}  '.format('> 10'))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[0])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_old[1])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[2])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[3])))
			f.write('{:>6}  '.format('%.3f'%(pl_fit_old[4])))
			f.write('{:>5}  '.format('%.3f'%(pl_fit_old[5])))
			f.close()

			logbins  = np.log10(bin_centers_pc_all)

			# write out the bin centers and correlation values
			f = open(data_dir + '%s/%s/%s_crosscf.nbins%02d.dat'%(gal_name, run_name, gal_name, nbins[j]), 'w')
			f.write('# cross correlation function values (1 + omega(theta)); bin centers in are given in log(pc)\n')
			
			f.write('{:<8}  '.format('nbins%02d'%nbins[j]))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(logbins[k])))
			
			f.write('\n')
			f.write('{:<8}  '.format('corr_all'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_all[k])))

			f.write('\n')
			f.write('{:<8}  '.format('corr_yng'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_young[k])))

			f.write('\n')
			f.write('{:<8}  '.format('corr_old'))
			for k in range(nbins[j]):
				f.write('{:>6}  '.format('%.3f'%(corr_old[k])))

			f.close()