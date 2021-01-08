"""
script for cross-correlation of the star cluster and gmc catalogs

2021-01-08
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
# matplotlib.use('agg')
# interactive plots
# matplotlib.use('Qt5agg')

def cross_correlate(sc_df, gmc_df, rand_sc_df, rand_gmc_df, dist, nbins=16, min_bin=1.1e-5, ):
	"""
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

	# looop through bins [drop first bin/min_bin]
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

def generate_random_sc(sc_df, xmax, ymax, wcs_hst, rseed=222):
	
	np.random.seed(rseed)

	# generate random distribution of star clusters over the entire HST image
	rand_x_full = np.random.uniform(0, xmax, size=4*len(sc_df))
	rand_y_full = np.random.uniform(0, ymax, size=4*len(sc_df))

	# limit to those in the hst-alma overlap mask
	rand_x_full_int = np.array([int(np.round(x)) for x in rand_x_full ])
	rand_y_full_int = np.array([int(np.round(y)) for y in rand_y_full ])

	rand_in_mask = np.array([True if mask[y,x] == 1 else False for y,x in zip(rand_y_full_int, rand_x_full_int)])

	rand_ra, rand_dec = wcs_hst.wcs_pix2world(rand_x_full[rand_in_mask], rand_y_full[rand_in_mask], 1.0)

	rand_sc_data = {'x': rand_x_full[rand_in_mask], 'y': rand_y_full[rand_in_mask], 'ra': rand_ra, 'dec': rand_dec }

	rand_sc_df = pd.DataFrame(rand_sc_data)

	return rand_sc_df

def generate_random_gmc(gmc_df, xmax, ymax, wcs_alma, wcs_hst, rseed=222):

	# generate random distribution of gmcs over the entire ALMA image
	rand_x_full = np.random.uniform(0, xmax, size=2*len(gmc_df))
	rand_y_full = np.random.uniform(0, ymax, size=2*len(gmc_df))

	# convert pixels to ra, dec
	rand_ra_full, rand_dec_full = wcs_alma.wcs_pix2world(rand_x_full, rand_y_full, 1.0)

	# conver ra, dec to HST pixels
	rand_x_full, rand_y_full = wcs_hst.wcs_world2pix(rand_ra_full, rand_dec_full, 1.0)

	# limit to those in the hst-alma overlap mask
	rand_x_full_int = np.array([int(np.round(x)) for x in rand_x_full ])
	rand_y_full_int = np.array([int(np.round(y)) for y in rand_y_full ])

	rand_in_mask = np.array([True if mask[y,x] == 1 else False for y,x in zip(rand_y_full_int, rand_x_full_int)])

	rand_ra, rand_dec = wcs_hst.wcs_pix2world(rand_x_full[rand_in_mask], rand_y_full[rand_in_mask], 1.0)

	rand_gmc_data = {'x': rand_x_full[rand_in_mask], 'y': rand_y_full[rand_in_mask], 'ra': rand_ra, 'dec': rand_dec }

	rand_gmc_df = pd.DataFrame(rand_gmc_data)

	return rand_gmc_df

def powerlaw_func(theta, Aw, alpha):
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
	
		# read in the HST-ALMA overlap mask 
		mask_hdu = fits.open(data_dir + '%s/%s_hst_alma_overlap_mask.fits'%(gal_name, gal_name))
		mask = mask_hdu[0].data
		mask_header = mask_hdu[0].header
		w_mask = WCS(mask_header, fobj=mask_hdu)
	
		ymax, xmax = np.shape(mask)
	
		# read in an ALMA image
		alma_hdulist = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_broad_mom0.fits'%(gal_name, gal_name))
		alma_header = alma_hdulist[0].header
		alma_data = alma_hdulist[0].data
		w_alma = WCS(alma_header, fobj=alma_hdulist, naxis=2)
	
		aymax, axmax = np.shape(alma_data)
		
		rand_sc_df = generate_random_sc(sc_df, xmax, ymax, w_mask)
		rand_gmc_df = generate_random_gmc(gmc_df, axmax, aymax, w_alma, w_mask)
	
		# # check up on our distributions
		# f = ap.FITSFigure(data_dir + '%s/hst/%s_uvis_f555w_exp_drc_sci.fits'%(gal_name, gal_name), figsize=(8,8))
		# f.recenter(gal_cen_ra[i], gal_cen_dec[i], radius=0.024)
	
		# f.show_markers(sc_df['ra'], sc_df['dec'], layer='Dsc', c='blue', marker='o', s=15)
		# f.show_markers(rand_sc_df['ra'], rand_sc_df['dec'], layer='Rsc', c='grey', marker='o', s=15)
		# f.show_markers(gmc_df['XCTR_DEG'], gmc_df['YCTR_DEG'], layer='Dgmc', c='red', marker='o', s=15)
		# f.show_markers(rand_gmc_df['ra'], rand_gmc_df['dec'], layer='Rgmc', c='orange', marker='o', s=15)
	
	
		bins, corr_all = cross_correlate(sc_df, gmc_df, rand_sc_df, rand_gmc_df, dist=dist, min_bin=1.1e-5)
	
		# get centers of bins
		bin_centers = 0.5 * (bins[1:] + bins[:-1])
		# drop nans
		wnnan = np.where(np.isnan(corr_all)==False)
		bin_centers_all = bin_centers[wnnan]
		corr_all = corr_all[wnnan]
		# bin centers as in pc
		bin_centers_pc_all = dist*1e6 * bin_centers_all*u.deg.to(u.rad)
	
		wyoung = sc_df['age'] <= 10
		rand_sc_df = generate_random_sc(sc_df.loc[wyoung], xmax, ymax, w_mask)
		rand_gmc_df = generate_random_gmc(gmc_df, axmax, aymax, w_alma, w_mask)
	
		bins, corr_young = cross_correlate(sc_df.loc[wyoung], gmc_df, rand_sc_df, rand_gmc_df, dist=dist, min_bin=1.1e-5)
		
		# drop nans
		wnnan = np.where(np.isnan(corr_young)==False)
		bin_centers_young = bin_centers[wnnan]
		corr_young = corr_young[wnnan]
		bin_centers_pc_young = dist*1e6 * bin_centers_young*u.deg.to(u.rad)
	
		wold = sc_df['age'] > 10
		rand_sc_df = generate_random_sc(sc_df.loc[wold], xmax, ymax, w_mask)
		rand_gmc_df = generate_random_gmc(gmc_df, axmax, aymax, w_alma, w_mask)
	
		bins, corr_old = cross_correlate(sc_df.loc[wold], gmc_df, rand_sc_df, rand_gmc_df, dist=dist, min_bin=1.1e-5)
		
		# drop nans
		wnnan = np.where(np.isnan(corr_old)==False)
		bin_centers_old = bin_centers[wnnan]
		corr_old = corr_old[wnnan]
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

		# output
		f = open(data_dir + '%s/%s_crosscf_fits.dat'%(gal_name, gal_name), 'w')
		f.write('#bin\t Aw\t\t  err\t\t Aw_pc\t\t err\t\t alpha\t err\n')
		f.write('all\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_all[0], perr_ang_all[0], popt_pc_all[0], perr_pc_all[0], popt_ang_all[1], perr_ang_all[1]))
		f.write('young\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_young[0], perr_ang_young[0], popt_pc_young[0], perr_pc_young[0], popt_ang_young[1], perr_ang_young[1]))
		f.write('old\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang_old[0], perr_ang_old[0], popt_pc_old[0], perr_pc_old[0], popt_ang_old[1], perr_ang_old[1]))
		f.close()

		fig, ax = plt.subplots(1,1, figsize=(5,5))
		ax.set_xscale('log')
		ax.set_yscale('log')
		# ax.set_xlabel(r'$\theta$ (deg)')
		ax.set_xlabel(r'$r$ [pc]')
		ax.set_ylabel(r'$1 + \omega(\theta)$')
	
		ax.plot(bin_centers_pc_all, corr_all, 'k-o', markersize=5, lw=1.5, label='All SCs (%i) \& GMCs (%i)'%(len(sc_df), len(gmc_df)))
		ax.plot(bin_centers_pc_young, corr_young, '-o', color='#377eb8', markersize=5, lw=1.5, label=r'$\leq 10$ Myr (%i) \& GMCs (%i)'%(len(sc_df.loc[wyoung]), len(gmc_df)))
		ax.plot(bin_centers_pc_old, corr_old, '-o', color='#e41a1c', markersize=5, lw=1.5, label=r'$> 10$ Myr (%i) \& GMCs (%i)'%(len(sc_df.loc[wold]), len(gmc_df)))

		plt.legend(loc='upper right', fontsize='x-small')

		plt.savefig(data_dir + '%s/%s_crosscf.png'%(gal_name, gal_name), bbox_inches='tight')
		plt.savefig(data_dir + '%s/%s_crosscf.pdf'%(gal_name, gal_name), bbox_inches='tight')
		plt.close()
