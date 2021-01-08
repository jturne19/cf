"""
script for running the cross-correlation function estimation

2021-01-07
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

from astroML.decorators import pickle_results
from astroML.correlation import bootstrap_two_point_angular

import sys
sys.path.append('/cherokee1/turner/phangs/cf/utils')
from utils import *

import matplotlib
# non-interactive plots
matplotlib.use('agg')
# interactive plots
# matplotlib.use('Qt5agg')


def auto_corr(df, min_bin=1.1e-5, nbins=10, nbootstraps=50, method='landy-szalay', rseed=222, gmc=False):
	
	np.random.seed(rseed)

	bins = 10 ** np.linspace(np.log10(min_bin), np.log10(10), nbins)
	results = [bins]

	if gmc:
		results += bootstrap_two_point_angular(df['XCTR_DEG'], df['YCTR_DEG'], bins=bins, method=method, Nbootstraps=nbootstraps)

	else:
		results += bootstrap_two_point_angular(df['ra'], df['dec'], bins=bins, method=method, Nbootstraps=nbootstraps)

	return results

def powerlaw_func(theta, Aw, alpha):
	return Aw * theta**alpha


data_dir = '/cherokee1/turner/phangs/cf/data/'

# read in list of galaxies 
master_galaxy_list = ascii.read('master_galaxy.list')
galaxy_list = ascii.read('galaxy.list')

gal_id 		= galaxy_list['id']
gal_alt_id 	= galaxy_list['alt_id']
gal_dist 	= galaxy_list['dist']

for i in range(len(galaxy_list)):

	# galaxy props
	gal_name = gal_id[i]
	dist = gal_dist[i]
	print('')
	print(gal_name)


	# read in the star cluster cat in the hst-alma footprint overlap mask
	sc_df = pd.read_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12_assoc_gmc.csv'%(gal_name, gal_name))

	# read in the gmc cat in the hst-alma footprint overlap mask
	gmc_cat = fits.open(data_dir + '%s/alma/%s_gmc_cat_masked.fits'%(gal_name, gal_name))[1].data
	gmc_df = Table(gmc_cat).to_pandas()

	f = open(data_dir + '%s/%s_tpcf_fits.dat'%(gal_name, gal_name), 'w')
	f.write('#bin\t Aw\t\t  err\t\t Aw_pc\t\t err\t\t alpha\t err\n')


	fig, ax = plt.subplots(1,1, figsize=(5,5))
	ax.set_xscale('log')
	ax.set_yscale('log')
	# ax.set_xlabel(r'$\theta$ (deg)')
	ax.set_xlabel(r'$r$ [pc]')
	ax.set_ylabel(r'$1 + \omega(\theta)$')

	# all clusters auto-correlation
	bins, corr, corr_err, bootstraps = auto_corr(sc_df, )

	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# drop nans
	wnnan = np.where(np.isnan(corr)==False)
	bin_centers = bin_centers[wnnan]
	corr = corr[wnnan] + 1
	corr_err = corr_err[wnnan]

	# bin centers as in pc
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad) 

	# power-law fit
	popt_ang, pcov = curve_fit(powerlaw_func, bin_centers, corr)
	popt_pc, pcov = curve_fit(powerlaw_func, bin_centers_pc, corr)
	perr_ang = np.sqrt(np.diag(pcov))
	perr_pc = np.sqrt(np.diag(pcov))
	f.write('all\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]))

	ax.errorbar(bin_centers_pc, corr, yerr=corr_err, fmt='k-o', ecolor='black', markersize=5, lw=1.5, label=r'All SCs $\alpha=%.2f\pm%.2f$ (%i) '%(popt_ang[1], perr_ang[1], len(sc_df)))

	# cluters <= 10 Myr
	wleq10 = sc_df['age'] <= 10

	bins, corr, corr_err, bootstraps = auto_corr(sc_df.loc[wleq10], )

	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# drop nans
	wnnan = np.where(np.isnan(corr)==False)
	bin_centers = bin_centers[wnnan]
	corr = corr[wnnan] + 1
	corr_err = corr_err[wnnan]

	# bin centers as in pc
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad) 

	popt_ang, pcov = curve_fit(powerlaw_func, bin_centers, corr)
	popt_pc, pcov = curve_fit(powerlaw_func, bin_centers_pc, corr)
	perr_ang = np.sqrt(np.diag(pcov))
	perr_pc = np.sqrt(np.diag(pcov))
	f.write('<=10\t %.3f\t %.3f\t %.3f\t%.3f\t\t %.3f\t %.3f\n'%(popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]))

	ax.errorbar(bin_centers_pc, corr, yerr=corr_err, fmt='-o', color='#377eb8', ecolor='#377eb8', markersize=5, lw=1.5, label=r'$\leq 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(popt_ang[1], perr_ang[1], len(sc_df.loc[wleq10])))

	# cluters > 100 Myr

	w10 = sc_df['age'] > 10

	bins, corr, corr_err, bootstraps = auto_corr(sc_df.loc[w10], )

	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# drop nans
	wnnan = np.where(np.isnan(corr)==False)
	bin_centers = bin_centers[wnnan]
	corr = corr[wnnan] + 1
	corr_err = corr_err[wnnan]

	# bin centers as in pc
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad) 

	popt_ang, pcov = curve_fit(powerlaw_func, bin_centers, corr)
	popt_pc, pcov = curve_fit(powerlaw_func, bin_centers_pc, corr)
	perr_ang = np.sqrt(np.diag(pcov))
	perr_pc = np.sqrt(np.diag(pcov))
	f.write('>10\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]))

	ax.errorbar(bin_centers_pc, corr, yerr=corr_err, fmt='-o', color='#e41a1c', ecolor='#e41a1c', markersize=5, lw=1.5, label=r'$> 10$ Myr $\alpha=%.2f\pm%.2f$ (%i) '%(popt_ang[1], perr_ang[1], len(sc_df.loc[w10])))


	# gmcs
	bins, corr, corr_err, bootstraps = auto_corr(gmc_df, min_bin=3e-4, gmc=True)

	bin_centers = 0.5 * (bins[1:] + bins[:-1])
	# drop nans
	wnnan = np.where(np.isnan(corr)==False)
	bin_centers = bin_centers[wnnan]
	corr = corr[wnnan] + 1
	corr_err = corr_err[wnnan]

	# bin centers as in pc
	bin_centers_pc = dist*1e6 * bin_centers*u.deg.to(u.rad) 

	popt_ang, pcov = curve_fit(powerlaw_func, bin_centers, corr)
	popt_pc, pcov = curve_fit(powerlaw_func, bin_centers_pc, corr)
	perr_ang = np.sqrt(np.diag(pcov))
	perr_pc = np.sqrt(np.diag(pcov))
	f.write('gmc\t\t %.3f\t %.3f\t\t %.3f\t\t %.3f\t\t %.3f\t %.3f\n'%(popt_ang[0], perr_ang[0], popt_pc[0], perr_pc[0], popt_ang[1], perr_ang[1]))
	f.close()

	ax.errorbar(bin_centers_pc, corr, yerr=corr_err, fmt='-o', color='#E68310', ecolor='#E68310', markersize=5, lw=1.5, label=r'GMCs $\alpha=%.2f\pm%.2f$ (%i) '%(popt_ang[1], perr_ang[1], len(gmc_df)))

	ax.axvline(gmc_df.mean()['RAD_PC'], lw=1.1, c='#999999', zorder=0)

	plt.legend(loc='upper right', fontsize='x-small')
	plt.savefig(data_dir + '%s/%s_tpcf.png'%(gal_name, gal_name), bbox_inches='tight')
	plt.savefig(data_dir + '%s/%s_tpcf.pdf'%(gal_name, gal_name), bbox_inches='tight')
	plt.close()

