"""
utilities/functions for creating plots/figures


"""
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
import aplpy as ap


from palettable.scientific.sequential import Tokyo_20
from palettable.scientific.sequential import LaJolla_20
from palettable.scientific.sequential import Turku_20_r

from palettable.cartocolors.sequential import SunsetDark_7
from palettable.cartocolors.sequential import agSunset_7
from palettable.cmocean.sequential import Speed_20
from palettable.cmocean.sequential import Thermal_20
from palettable.cmocean.sequential import Turbid_20
from palettable.cmocean.sequential import Matter_20
from palettable.cmocean.sequential import Oxy_20
from palettable.colorbrewer.sequential import Oranges_9
from palettable.colorbrewer.sequential import Reds_9
from palettable.colorbrewer.sequential import YlGn_9_r
from palettable.colorbrewer.sequential import YlOrRd_9


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