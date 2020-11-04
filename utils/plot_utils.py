"""
utilities/functions for creating plots/figures


"""
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
import aplpy as ap

def outline_plot(gal_name, data_dir, sc_coord_dict, center_ra, center_dec, radius=0.04):
	""" create an 'outline plot' for the given galaxy
	plot using the wcs info in the HST images
	shows the ellipses of the GMCs 
	shows the position of the star clusters 
	
	Inputs:
	gal_name		str 	name of the galaxy e.g., ngc0628
	data_dir 		str 	path to the directory containing the 'hst' and 'alma' directories
	sc_coord_dict 	dict 	dictionary containing the star cluster coordinates: 'x', 'y', 'ra', 'dec'


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


	fig = ap.FITSFigure(hst_image)
	fig.recenter(center_ra, center_dec, radius=radius)

	fig.show_regions(gmc_region_file, layer='gmc')
	fig.show_regions(footprint_file, layer='footprint')
	fig.show_markers(sc_coord_dict['ra'], sc_coord_dict['dec'], layer='clusters', marker='o', s=15, facecolor='red')

	fig.axis_labels.hide_y()
	fig.axis_labels.hide_x()

	fig.tick_labels.set_font(size='xx-large', family='serif')
	fig.set_theme('publication')
	fig.tick_labels.set_xformat('hh:mm:ss')
	fig.tick_labels.set_yformat('dd:mm:ss')

	fig.save(data_dir + '%s/%s_outlineplot.png'%(gal_name, gal_name))
	fig.save(data_dir + '%s/%s_outlineplot.pdf'%(gal_name, gal_name))

	fig.close()


def all_galaxies_outline_plots(galaxy_list, data_dir, sc_class='class12', radius=[0.04]):
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
		center_ra = gal_ra[i]
		center_dec = gal_dec[i]

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

		outline_plot(gal_name, data_dir, sc_coord_dict, center_ra=center_ra, center_dec=center_dec, radius=gal_plot_radius)