"""
script for cross-correlation of the star cluster and gmc catalogs

loops through all galaxies given in the galaxy_list
reads in the cluster catalog and the GMC catalog
reads in the HST-ALMA overlap mask image which is needed for generating the random catalogs
generates random SC catalog/positions within the overlap mask
generates random GMC catalog/positions within the overlap mask
does cross-correlation functions for star clusters - GMCs
does cross-correlation functions for just young (<= 10 Myr) star clusters - GMCs
does cross-correlation functions for just old (> 10 Myr) star clusters - GMCs
does powerlaw-fits to the cross-correlation functions
outputs the best-fit powerlaw params to a text file
outputs figures of the cross-correlation functions

2021-02-11
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
from scipy.optimize import curve_fit

import sys
sys.path.append('/cherokee1/turner/phangs/cf/utils')
from utils import *

import matplotlib
# non-interactive plots
matplotlib.use('agg')
# interactive plots
# matplotlib.use('Qt5agg')

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

	bins = 10 ** np.linspace(np.log10(min_bin), np.log10(10), nbins)

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
	sc_df 		pandas DataFrame 		dataframe which contains the star cluster catalog; needs to include the cluster ra, dec
	gmc_df 		pandas DataFrame 		dataframe which contains the gmc catalog; needs to include the gmc cluster ra, dec
	wcs_hst 	astropy.wcs.wcs.WCS 	astropy wcs object of the hst image
	xmax_hst 	int 					the maximum pixel location in the x-direction of the hst image; i.e., if hst image is 13000x14000 pixels, xmax_hst = 13000
	ymax_hst 	int 					the maximum pixel location in the y-direction of the hst image; i.e., if hst image is 13000x14000 pixels, ymax_hst = 14000
	wcs_alma 	astropy.wcs.wcs.WCS 	astropy wcs object of the alma image
	xmax_alma 	int 					the maximum pixel location in the x-direction of the alma image; i.e., if alma image is 1599x1598 pixels, xmax_alma = 1599
	ymax_alma	int 					the maximum pixel location in the y-direction of the alma image; i.e., if alma image is 1599x1598 pixels, ymax_alma = 1598
	mask 		array 					the data array of the hst-alma overlap mask image
	dist 		float 					distance to the galaxy in Mpc; needed for calculating the separations b/w clusters and gmcs in pc
	nbootstraps int 					number of bootstraps to run through
	**kwargs 	dict 					keyword arguments to pass on to the cross_correlate funciton; really just nbins and min_bin

	Outputs:
	bins 		array 					array of the angular bin edges (degrees)
	corr 		array 					correlation values in each bin --> 1 + omega(theta)
	corr_err 	array 					1 sigma error on the correlation values; if correlation is nan, error will be 0

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

	return bins, corr, corr_err	



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

def powerlaw_func(theta, Aw, alpha):
	""" a powerlaw function of the form 
	f(theta) = Aw * theta^alpha

	"""
	return Aw * theta**alpha

if __name__ == '__main__':

	data_dir = '/cherokee1/turner/phangs/cf/data/'
	# read in list of galaxies 
	master_galaxy_list = ascii.read('master_galaxy.list')
	galaxy_list = ascii.read('galaxy.list')
	
	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	gal_cen_ra  = galaxy_list['center_ra']
	gal_cen_dec = galaxy_list['center_dec']
	
	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):
	
		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)
	
		# read in the star cluster dataframe
		sc_df = pd.read_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12_assoc_gmc.csv'%(gal_name, gal_name))
	
		# read in the gmc cat
		gmc_cat = fits.open(data_dir + '%s/alma/%s_gmc_cat_masked.fits'%(gal_name, gal_name))[1].data
		gmc_df = Table(gmc_cat).to_pandas()
	
		# read in the HST-ALMA overlap mask (uses same wcs info as the HST image)
		mask_hdu = fits.open(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name))
		mask = mask_hdu[0].data
		mask_header = mask_hdu[0].header
		wcs_mask = WCS(mask_header, fobj=mask_hdu)
		# max pixel numbers
		ymax_hst, xmax_hst = np.shape(mask)
	
		# read in an ALMA image
		alma_hdulist = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))
		alma_header = alma_hdulist[0].header
		alma_data = alma_hdulist[0].data
		wcs_alma = WCS(alma_header, fobj=alma_hdulist, naxis=2)
		ymax_alma, xmax_alma = np.shape(alma_data)
		
		# cross correlation (bootstrap error estimatation) for all star clusters/gmcs
		bins, corr_all, corr_all_err = cross_correlate_bootstrap(sc_df, gmc_df, wcs_mask, xmax_hst, ymax_hst, wcs_alma, xmax_alma, ymax_alma, mask, dist, nbootstraps=100, nbins=16, min_bin=1.1e-5 )
		
		# get centers of bins
		bin_centers = 0.5 * (bins[1:] + bins[:-1])
		# drop nans
		wnnan = np.where(np.isnan(corr_all)==False)
		bin_centers_all = bin_centers[wnnan]
		corr_all = corr_all[wnnan]
		corr_all_err = corr_all_err[wnnan]
		# bin centers in pc
		bin_centers_pc_all = dist*1e6 * bin_centers_all*u.deg.to(u.rad)
		
		# cross correlation with only the young star clusters (<= 10 Myr)
		wyoung = sc_df['age'] <= 10
		bins, corr_young, corr_young_err = cross_correlate_bootstrap(sc_df[wyoung], gmc_df, wcs_mask, xmax_hst, ymax_hst, wcs_alma, xmax_alma, ymax_alma, mask, dist, nbootstraps=100, nbins=16, min_bin=1.1e-5 )
		
		# drop nans
		wnnan = np.where(np.isnan(corr_young)==False)
		bin_centers_young = bin_centers[wnnan]
		corr_young = corr_young[wnnan]
		corr_young_err = corr_young_err[wnnan]
		bin_centers_pc_young = dist*1e6 * bin_centers_young*u.deg.to(u.rad)
	
		# cross correlation with only the old star clusters (> 10 Myr)
		wold = sc_df['age'] > 10
		bins, corr_old, corr_old_err = cross_correlate_bootstrap(sc_df[wold], gmc_df, wcs_mask, xmax_hst, ymax_hst, wcs_alma, xmax_alma, ymax_alma, mask, dist, nbootstraps=100, nbins=16, min_bin=1.1e-5 )
		
		# drop nans
		wnnan = np.where(np.isnan(corr_old)==False)
		bin_centers_old = bin_centers[wnnan]
		corr_old = corr_old[wnnan]
		corr_old_err = corr_old_err[wnnan]
		bin_centers_pc_old = dist*1e6 * bin_centers_old*u.deg.to(u.rad)
		
		# fit the power law
		popt_ang_all, pcov = curve_fit(powerlaw_func, bin_centers_all, corr_all)
		popt_pc_all, pcov = curve_fit(powerlaw_func, bin_centers_pc_all, corr_all)
		perr_ang_all = np.sqrt(np.diag(pcov))
		perr_pc_all = np.sqrt(np.diag(pcov))

		popt_ang_young, pcov = curve_fit(powerlaw_func, bin_centers_young, corr_young)
		popt_pc_young, pcov = curve_fit(powerlaw_func, bin_centers_pc_young, corr_young)
		perr_ang_young = np.sqrt(np.diag(pcov))
		perr_pc_young = np.sqrt(np.diag(pcov))

		popt_ang_old, pcov = curve_fit(powerlaw_func, bin_centers_old, corr_old)
		popt_pc_old, pcov = curve_fit(powerlaw_func, bin_centers_pc_old, corr_old)
		perr_ang_old = np.sqrt(np.diag(pcov))
		perr_pc_old = np.sqrt(np.diag(pcov))

		# output text file with all the power law best-fits
		f = open(data_dir + '%s/%s_crosscf_fits.dat'%(gal_name, gal_name), 'w')
		f.write('#bin\t Aw\t\t  err\t\t Aw_pc\t\t err\t\t alpha\t err\n')
		f.write('all\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_all[0], perr_ang_all[0], popt_pc_all[0], perr_pc_all[0], popt_ang_all[1], perr_ang_all[1]))
		f.write('young\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_young[0], perr_ang_young[0], popt_pc_young[0], perr_pc_young[0], popt_ang_young[1], perr_ang_young[1]))
		f.write('old\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_old[0], perr_ang_old[0], popt_pc_old[0], perr_pc_old[0], popt_ang_old[1], perr_ang_old[1]))
		f.close()

		# generate figures
		fig, ax = plt.subplots(1,1, figsize=(5,5))
		ax.set_xscale('log')
		ax.set_yscale('log')
		# ax.set_xlabel(r'$\theta$ (deg)')
		ax.set_xlabel(r'$r$ [pc]')
		ax.set_ylabel(r'$1 + \omega(\theta)$')
	
		ax.errorbar(bin_centers_pc_all, corr_all, yerr=corr_all_err, fmt='k-o', markersize=5, lw=1.5, label='All SCs (%i) \& GMCs (%i)'%(len(sc_df), len(gmc_df)))
		ax.errorbar(bin_centers_pc_young, corr_young, yerr=corr_young_err, fmt='-o', color='#377eb8', markersize=5, lw=1.5, label=r'$\leq 10$ Myr (%i) \& GMCs (%i)'%(len(sc_df.loc[wyoung]), len(gmc_df)))

		ax.errorbar(bin_centers_pc_old, corr_old, yerr=corr_old_err, fmt='-o', color='#e41a1c', markersize=5, lw=1.5, label=r'$> 10$ Myr (%i) \& GMCs (%i)'%(len(sc_df.loc[wold]), len(gmc_df)))

		plt.legend(loc='upper right', fontsize='x-small')

		plt.savefig(data_dir + '%s/%s_crosscf.png'%(gal_name, gal_name), bbox_inches='tight')
		plt.savefig(data_dir + '%s/%s_crosscf.pdf'%(gal_name, gal_name), bbox_inches='tight')
		plt.close()
