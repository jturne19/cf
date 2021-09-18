"""
script for making the corner plots
taken out from inside the pipeline for NGC3351 so it might be a little weird

To run, you need:

line 39: directory = path to the output directory e.g., 'cigale/out_run01'; will need to grab the observations.fits from it

line  43: bc03 track in color-color space -- I made it from running cigale in the model-ouput mode
line  50: header file that has the photfnu and vega zero point
line  61: original phangs-hst catalog fits table
line  65: 'models' array is shaped like: models[:,0] = ages, models[:,1] = masses, models[:,2] = E(B-V) reddenings
	      each row in models is a particular model on the model-grid used in the SED fitting	
line  67: 'models_prior' is the same as above but with the krumholz fiducial prior applied
line 148: 'age_chi2min' is the best-fit (chi2 minimized) age
line 149: 'age_mean_prior' is the best-estimate (mean) from the prior-applied models
line 150: 'mass_chi2min', 'mass_mean_prior' same as above but for masses
line 152: 'ebv_chi2min', 'ebv_mean_prior' same as above but for ebv reddenings
line 184: postage stamp of the cluster as a png. comment out lines 184-187 if you don't have postage stamp

"""
# not sure if all these imports are necessary but just in case 
import numpy as np
from astropy.io import fits
from statsmodels.stats.weightstats import DescrStatsW
import glob
import os.path
from astropy.table import Table
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patheffects as pe
import matplotlib as mpl
from matplotlib import colors
import seaborn as sns
from scipy import stats

# read in the observations to make  the color-color plot
fname = '%s/observations.fits'%(directory)
obsdata = fits.open(fname)[1].data

# read in the bc03 track
fname = 'cigale/sfh_results/out_sfh_ssp/models-block-0.fits'
track = fits.open(fname)[1].data
# track fluxes [mJy]
trackflux = np.array([track['F336W_UVIS_CHIP2'],track['F438W_UVIS_CHIP2'], track['F555W_UVIS_CHIP2'],track['F814W_UVIS_CHIP2']])
# convert to [Jy]
trackflux *= 1e-3
# header info for photfnu and zero points
photfnu, zpvega = np.loadtxt('photometry/header_info_ngc3351_prime.txt', usecols=[4,6], unpack=True)
# keep the right filters
photfnu = photfnu[[0,1,2,4,6]]
zpvega = zpvega[[0,1,2,4,6]]
# convert to vega mags
trackmag = np.array([ -2.5 * np.log10(f/pfn) + zp for f,pfn,zp in zip(trackflux, photfnu, zpvega) ])
# track colors
trackred  = trackmag[2] - trackmag[3]
trackblue = trackmag[0] - trackmag[1]

# read in original catalog to get x, y, ra, dec positions
fname = 'n3351/ngc3351_phangshst_base_catalog_v1.class12.fits'
catalog = fits.open(fname)[1].data

# log the ages and masses
models[:,0] = np.log10(models[:,0])
models[:,1] = np.log10(models[:,1])
models_prior[:,0] = np.log10(models_prior[:,0])
models_prior[:,1] = np.log10(models_prior[:,1])

# create the figure with subplots
fig, axs = plt.subplots(3,3, figsize=(12,12))
# adjust subplot spacing
plt.subplots_adjust(wspace=0.2, hspace=0.05)

# create color bars
cm2 = plt.cm.get_cmap('Blues')

zmax = []

# looping through the 3x3 subplots
for i in range(3):
	for j in range(3):
		# 2d histograms
		if j < i :
			plot_bins = 50

			kde2d = stats.gaussian_kde([models_prior[:,j], models_prior[:,i]], bw_method=None, weights=plike)
			xi, yi = np.mgrid[models_prior[:,j].min():models_prior[:,j].max():plot_bins*1j, models_prior[:,i].min():models_prior[:,i].max():plot_bins*1j]
			zi = kde2d(np.vstack([xi.flatten(), yi.flatten()]))
			zmax.append(zi.max())
			norm = colors.Normalize(vmin=1e-2, vmax=zi.max())

			h2 = axs[i,j].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cm2, vmin=1e-2, norm=norm)

		# 1d histograms along the diagonals
		elif i == j :

			kde = stats.gaussian_kde(models_prior[:,i], bw_method=None, weights=plike)
			x = np.linspace(models_prior[:,j].min(), models_prior[:,j].max(), 500)
			y = kde(x)
			w_nonzero = np.where(y > 1e-6)[0]
			x, y = x[w_nonzero], y[w_nonzero]

			axs[i,j].plot(x, y, linewidth=2.5, color='#377eb8')
			# axs[i,j].hist(models_prior[:,i], weights=plike, bins=bins, density=True, color='#969696', linewidth=2.3)


# set age axis limits
axs[0,0].set_xlim(-0.3, 4.3)
axs[1,0].set_xlim(-0.3, 4.3)
axs[2,0].set_xlim(-0.3, 4.3)

# set mass axis limits
axs[1,1].set_xlim(models_prior[:,1].min() - 0.2, models_prior[:,1].max() + 0.2)
axs[2,1].set_xlim(models_prior[:,1].min() - 0.2, models_prior[:,1].max() + 0.2)
axs[1,0].set_ylim(models_prior[:,1].min() - 0.2, models_prior[:,1].max() + 0.2)

# set ebv axis limits
axs[2,2].set_xlim(-0.1, 1.6)
axs[2,0].set_ylim(-0.1, 1.6)
axs[2,1].set_ylim(-0.1, 1.6)

# subplot panel labels
axs[2,0].set_xlabel(r'log(Age) [Myr]')
axs[2,1].set_xlabel(r'log(Mass) [M$_{\odot}$]')
axs[2,2].set_xlabel(r'E(B-V)')
axs[0,0].set_ylabel(r'PDF')
axs[1,0].set_ylabel(r'log(Mass) [M$_{\odot}$]')
axs[2,0].set_ylabel(r'E(B-V)')

# turn off tick labels for the appropriate panels
axs[0,0].set_xticklabels('')
axs[1,0].set_xticklabels('')
axs[1,1].set_xticklabels('')
axs[2,1].set_yticklabels('')
# turn the top right and middle right panels off
axs[1,2].axis('off')
axs[0,2].axis('off')

# place the colorbars for the 2d histograms in the upper right panel
norm = colors.Normalize(vmin=0, vmax=zmax[0])
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cm2), ax=axs[0,2], orientation='horizontal', pad=0.1, panchor=(0.5, 0.5))
cbar.ax.xaxis.set_label_position('top')
cbar.ax.set_xlabel('2D PDF')

# plot lines for the best-fits
# chi2min for the default, bayes likelihood-weighted mean for the prior
axs[0,0].axvline(np.log10(age_chi2min),     color='black',   linestyle='--', linewidth=1.5)
axs[0,0].axvline(np.log10(age_mean_prior),  color='#e41a1c', linestyle='--', linewidth=1.5)
axs[1,1].axvline(np.log10(mass_chi2min),    color='black',   linestyle='--', linewidth=1.5)
axs[1,1].axvline(np.log10(mass_mean_prior), color='#e41a1c', linestyle='--', linewidth=1.5)
axs[2,2].axvline(ebv_chi2min,     			color='black',   linestyle='--', linewidth=1.5)
axs[2,2].axvline(ebv_mean_prior,  			color='#e41a1c', linestyle='--', linewidth=1.5)

# text for displaying the cluster's position
wclust = np.where(catalog['ID_PHANGS_ALLSOURCES_V0_9'] == cluster_id)[0]
clust = catalog[wclust][0]

plt.figtext(0.68, 0.86, 'Cluster %i'%cluster_id, fontsize='large')
plt.figtext(0.68, 0.84, r'%.6f$^\circ$, %.6f$^\circ$'%(clust['PHANGS_RA'], clust['PHANGS_DEC']), fontsize='large')
plt.figtext(0.68, 0.82, '%.2f, %.2f'%(clust['PHANGS_X'], clust['PHANGS_Y']), fontsize='large')

# color-color diagram in top middle panel
# grab just the one cluster's data
clust_obs = obsdata[wclust][0]
flux = np.array([clust_obs['F336W_UVIS_CHIP2'],clust_obs['F438W_UVIS_CHIP2'],clust_obs['F555W_UVIS_CHIP2'],clust_obs['F814W_UVIS_CHIP2']])
# convert to mags
mag  = np.array([ -2.5 * np.log10(f/pfn) + zp for f,pfn,zp in zip(flux, photfnu, zpvega) ])
red  = mag[2] - mag[3]
blue = mag[0] - mag[1]

# plot the track
axs[0,1].plot(trackred, -1*trackblue, 'k-', linewidth=1.8)
# plot cluster color
axs[0,1].scatter(red, -1*blue, marker='D', color='#377eb8', s=30, edgecolor='black', lw=0.5)
axs[0,1].xaxis.set_label_position('top')
axs[0,1].yaxis.set_label_position('right')
axs[0,1].tick_params(labeltop=True, labelbottom=False, labelleft=False, labelright=True)
axs[0,1].set_xlabel('F555 - F814')
axs[0,1].set_ylabel('F336 - F438', rotation=-90, labelpad=15)
axs[0,1].set_yticklabels(['2.0','1.5','1.0','0.5','0.0','-0.5'])

# postage stamp in the middle right panel
fname = 'n3351/phangs/maps/postage/v1/noframe/ps_noframe{:0<6}.png'.format('%06d'%cluster_id)

img = mpimg.imread(fname)
axs[1,2].imshow(img)

plt.savefig('%s'%figurename + '{:<06}.png'.format('%06d'%cluster_id), bbox_inches='tight')
plt.savefig('%s'%figurename + '{:<06}.pdf'.format('%06d'%cluster_id), bbox_inches='tight')
plt.close()