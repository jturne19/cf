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

########################################################
# making figures/plots functions
########################################################

def outline_plot(gal_name, data_dir, sc_coord_dict, center_deg, radius=0.04, bkgd=None, color_arr=[], color_code='' ):
	""" create an 'outline plot' for the given galaxy
	plot using the wcs info in the HST images
	shows the ellipses of the GMCs 
	shows the position of the star clusters 
	
	Inputs:
	gal_name		str 		name of the galaxy e.g., ngc0628
	data_dir 		str 		path to the directory containing the 'hst' and 'alma' directories
	sc_coord_dict 	dict 		dictionary containing the star cluster coordinates: 'x', 'y', 'ra', 'dec'
	center_deg 		tuple 		(RA,DEC) coordinate of the center of the galaxy in degrees
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
		f.save(data_dir + '%s/%s_outlineplot_%s.png'%(gal_name, gal_name, bkgd))
		f.save(data_dir + '%s/%s_outlineplot_%s.pdf'%(gal_name, gal_name, bkgd))
	elif len(color_arr) > 0:
		fig.savefig(data_dir + '%s/%s_outlineplot_%s.png'%(gal_name, gal_name, color_code), bbox_inches='tight')
		fig.savefig(data_dir + '%s/%s_outlineplot_%s.pdf'%(gal_name, gal_name, color_code), bbox_inches='tight')

	else:
		f.save(data_dir + '%s/%s_outlineplot.png'%(gal_name, gal_name))
		f.save(data_dir + '%s/%s_outlineplot.pdf'%(gal_name, gal_name))


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


		outline_plot(gal_name, data_dir, sc_coord_dict, center_deg, radius=gal_plot_radius, bkgd=bkgd, color_arr=color_input, color_code=color_code )