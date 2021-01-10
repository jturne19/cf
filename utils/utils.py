"""
master file which contains all the utilities/functions 

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
			base_cat = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.fits'%(gal_name,gal_name))[1].data
	
		except FileNotFoundError:
			print(data_dir + '%s/hst/%s_phangshst_base_catalog.fits not found, skipping'%(gal_name,gal_name))
			continue
		
		n_total    = len(base_cat)
		print('total number clusters = %i'%(n_total))
		
		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in base_cat.names:
			
			# make new fits files with just class 1,2,3 and 1,2
			if mkfits:
				new_clust_fits_from_classes(base_cat, [1,2,3], data_dir + '%s/hst/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))
				new_clust_fits_from_classes(base_cat, [1,2], data_dir + '%s/hst/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))
	
			# read those in
			cat123 = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))[1].data
			cat12  = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))[1].data
			
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
			cat = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.fits'%(gal_name,gal_name))[1].data
		except FileNotFoundError:
			print(data_dir + '%s/hst/%s_phangshst_base_catalog.fits not found, skipping'%(gal_name,gal_name))
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
	
		mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + '%s/hst/%s_allclusters'%(gal_name,gal_name), color='green')
	
		# check if the catalog has the cluster classifications
		# and skip galaxies missing the classifications
		if 'PHANGS_CLUSTER_CLASS' in cat.names:
			
			cat123 = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class123.fits'%(gal_name,gal_name))[1].data
			
			x	= cat123['PHANGS_X']
			y	= cat123['PHANGS_Y']
			ra 	= cat123['PHANGS_RA']
			dec	= cat123['PHANGS_DEC']
			coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
			mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + '%s/hst/%s_class123'%(gal_name,gal_name), color='green')
	
			cat12  = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class12.fits'%(gal_name,gal_name))[1].data
	
			x	= cat123['PHANGS_X']
			y	= cat123['PHANGS_Y']
			ra 	= cat123['PHANGS_RA']
			dec	= cat123['PHANGS_DEC']
			coord_sc = {'x': x, 'y': y, 'ra': ra, 'dec': dec}
	
			mk_clust_ds9_regions(coord_sc, radius_pix, radius_asec, data_dir + '%s/hst/%s_class12'%(gal_name,gal_name), color='green')


########################################################
# GMC catalog specific functions
########################################################

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

def generate_gmc_cat_masked(galaxy_list, data_dir):
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
		gmc_cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits'%(gal_name, gal_name))[1].data

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
		t.write(data_dir + '%s/alma/%s_gmc_cat_masked.fits'%(gal_name, gal_name), format='fits', overwrite=True)




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

def outline_plot(gal_name, data_dir, sc_coord_dict, center_deg, sc_class='class12', radius=0.04, bkgd=None, color_arr=[], color_code='' ):
	""" create an 'outline plot' for the given galaxy
	plot using the wcs info in the HST images
	shows the ellipses of the GMCs 
	shows the position of the star clusters 
	
	Inputs:
	gal_name		str 		name of the galaxy e.g., ngc0628
	data_dir 		str 		path to the directory containing the 'hst' and 'alma' directories
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
	gmc_cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits'%(gal_name,gal_name))[1].data
	
	# path to GMC ellipses/region file
	gmc_region_file = data_dir + '%s/alma/%s_gmc_cat.blue.deg.reg'%(gal_name, gal_name)

	# path to HST image file 
	hst_image = data_dir + '%s/hst/%s_uvis_f555w_exp_drc_sci.fits'%(gal_name, gal_name)

	# path to full HST footprint file
	footprint_file = data_dir + '%s/hst/%s_full_footprint.pix.reg'%(gal_name, gal_name)



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
		f.save(data_dir + '%s/%s_%s_outlineplot_%s.png'%(gal_name, gal_name, sc_class, bkgd))
		f.save(data_dir + '%s/%s_%s_outlineplot_%s.pdf'%(gal_name, gal_name, sc_class, bkgd))
	elif len(color_arr) > 0:
		fig.savefig(data_dir + '%s/%s_%s_outlineplot_%s.png'%(gal_name, gal_name, sc_class, color_code), bbox_inches='tight')
		fig.savefig(data_dir + '%s/%s_%s_outlineplot_%s.pdf'%(gal_name, gal_name, sc_class, color_code), bbox_inches='tight')

	else:
		f.save(data_dir + '%s/%s_%s_outlineplot.png'%(gal_name, gal_name, sc_class))
		f.save(data_dir + '%s/%s_%s_outlineplot.pdf'%(gal_name, gal_name, sc_class))


	plt.close()


def all_galaxies_outline_plots(galaxy_list, data_dir, sc_class='class12', radius=[0.04], bkgd=None, color_code=''):
	"""

	"""
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	gal_ra		= galaxy_list['center_ra']
	gal_dec 	= galaxy_list['center_dec']

	for i in range(len(galaxy_list)):

		gal_name = gal_id[i]
		gal_alt_name = gal_alt_id[i]
		center_deg = (gal_ra[i], gal_dec[i])

		print('')
		print(gal_name)

		# read in the correct cluster catalog
		if 'class12' in sc_class:
			sc_cat = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.%s.fits'%(gal_name, gal_name, sc_class))[1].data
		else:
			# default to the class 1,2 catalog
			sc_cat = fits.open(data_dir + '%s/hst/%s_phangshst_base_catalog.class12.fits'%(gal_name, gal_name))[1].data

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


		outline_plot(gal_name, data_dir, sc_coord_dict, center_deg, sc_class=sc_class, radius=gal_plot_radius, bkgd=bkgd, color_arr=color_input, color_code=color_code )


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


def bootstrap_median_error(data, sigma, nbootstraps=10000):
	"""	boostrap estimate of the error on the median of the distribution of input data and their errors/sigma
		
	"""

	medians = np.zeros(nbootstraps)

	for i in range(nbootstraps):

		rand = np.random.normal(data, sigma)

		medians[i] = np.median(rand)

	std = np.std(medians)

	return std


