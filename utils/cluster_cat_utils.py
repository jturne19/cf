"""
utilities/functions for manipulating the star cluster catalogs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u

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


def all_galaxies_clust_cats(galaxy_list, data_dir, mkfits=False):
	""" loop through all the given glaxies and print out the cluster numbers
	and make new fits tables with class 1,2,3 and class 1,2

	"""
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	
	
	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
		gal_alt_name = gal_alt_id[i]
		
		print('')
		print(gal_name)
	
		# check if the base catalog fits files exists
		try:
			base_cat = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits'%(gal_name,gal_name))[1].data
	
		except FileNotFoundError:
			print(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits not found, skipping'%(gal_name,gal_name))
			continue
		
		n_total    = len(base_cat)
		print('total number clusters = %i'%(n_total))
		
		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in base_cat.names:
			
			# make new fits files with just class 1,2,3 and 1,2
			if mkfits:
				new_clust_fits_from_classes(base_cat, [1,2,3], data_dir + 'hst/%s/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))
				new_clust_fits_from_classes(base_cat, [1,2], data_dir + 'hst/%s/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))
	
			# read those in
			cat123 = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))[1].data
			cat12  = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))[1].data
			
			n_class123 = len(cat123)
			n_class12  = len(cat12)
		
			print('number of class 1,2,3 = %i'%(n_class123))
			print('number of class 1,2   = %i'%(n_class12))
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

def all_galaxies_clust_region_files(galaxy_list, data_dir, radius_pix=10):
	""" loop through all the given galaxies and make ds9 region files for the star clusters
	makes separate region files for class 1,2,3 and class 1,2 if available
	makes hst image pixel region files and degree region files

	radius_pix sets the circle radius in pixels, default is 10 pixels
	
	"""

	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):
	
		gal_name = gal_id[i]
	
		print(gal_name)
	
		# check if the catalog fits files exists
		try:
			cat = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits'%(gal_name,gal_name))[1].data
		except FileNotFoundError:
			print(data_dir + 'hst/%s/%s_phangshst_base_catalog.fits not found, skipping'%(gal_name,gal_name))
			continue
	
		# pull out sc coordinates 
		x	= cat['PHANGS_X']
		y	= cat['PHANGS_Y']
		ra 	= cat['PHANGS_RA']
		dec	= cat['PHANGS_DEC']
	
		coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
		# read in image header to get cdelt so we can convert pixel to arcsec
		# use 275 because its UVIS for all of them
		image_hdr = fits.getheader(data_dir + 'hst/%s/%s_uvis_f275w_exp_drc_sci.fits'%(gal_name, gal_name))
	
		# degrees per pix
		cdelt = image_hdr['CD2_2']
		# convert input pixel radius to degrees
		radius_deg = radius_pix * cdelt * u.deg
		# conver to arcsec
		radius_asec = radius_deg.to(u.arcsec).value
	
		mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + 'hst/%s/%s_allclusters'%(gal_name,gal_name), color='green')
	
		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in cat.names:
			
			cat123 = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))[1].data
			
			x	= cat123['PHANGS_X']
			y	= cat123['PHANGS_Y']
			ra 	= cat123['PHANGS_RA']
			dec	= cat123['PHANGS_DEC']
			coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
			mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + 'hst/%s/%s_class123'%(gal_name,gal_name), color='green')
	
			cat12  = fits.open(data_dir + 'hst/%s/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))[1].data
	
			x	= cat123['PHANGS_X']
			y	= cat123['PHANGS_Y']
			ra 	= cat123['PHANGS_RA']
			dec	= cat123['PHANGS_DEC']
			coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
			mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + 'hst/%s/%s_class12'%(gal_name,gal_name), color='green')
