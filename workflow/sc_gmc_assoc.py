"""
script for finding the star clusters associated with each gmc

2020-01-05
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

import sys
sys.path.append('/cherokee1/turner/phangs/cf/utils')
from utils import *

import matplotlib
# non-interactive plots
matplotlib.use('agg')
# interactive plots
# matplotlib.use('Qt5agg')


def generate_sc_gmc_assoc_df(galaxy_list, data_dir):

	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']

	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):
	
		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)
	
		# read in the csv of the cluster catalog within the overlap mask
		sc_df = pd.read_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12.csv'%(gal_name, gal_name))
	
		# read in the gmc cat
		gmc_cat = fits.open(data_dir + '%s/alma/%s_12m+7m+tp_co21_native_props.fits'%(gal_name, gal_name))[1].data
	
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
	
		# insert those new columns
		sc_df2['assoc_gmc_cloudnum'] = assoc_gmc_cloudnum
		sc_df2['assoc_num'] = assoc_num
	
		# output to csv file
		sc_df2.to_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12_assoc_gmc.csv'%(gal_name, gal_name), index=False)

def make_sc_gmc_assoc_hists(galaxy_list, data_dir):
	""" loop through all the galaxies and make the histograms of the cluster ages with their gmc associations

	"""

	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	
	for i in range(len(galaxy_list)):

		# galaxy props
		gal_name = gal_id[i]
		dist = gal_dist[i]
		print('')
		print(gal_name)

		# read in the invidual galaxy dataframe rather than using the mega
		df = pd.read_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12_assoc_gmc.csv'%(gal_name, gal_name))

		w0 = df['assoc_num'] == 0
		w1 = df['assoc_num'] == 1
		w2 = df['assoc_num'] == 2
		w3 = df['assoc_num'] == 3

		sc_gmc_assoc_hist(df, filename=data_dir+'%s/%s_sc_gmc_assoc_hist'%(gal_name, gal_name))

		# star cluster ages
		age_all = df['age'].to_numpy()
		lage_all = np.log10(age_all)
		age_err_all = df['age_err'].to_numpy()
		lage_err_all = age_err_all/age_all/np.log(10)

		# log the stats for each galaxy
		if i == 0:
			f = open(data_dir + 'sc_gmc_assoc_stats.txt', 'w')
			f.write(gal_name + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), 1.2533*np.std(age_all)/np.sqrt(len(age_all)) ) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), 1.2533*np.std(age_all[w1])/np.sqrt(len(age_all[w1])) ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), 1.2533*np.std(age_all[w2])/np.sqrt(len(age_all[w2])) ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), 1.2533*np.std(age_all[w3])/np.sqrt(len(age_all[w3])) ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), 1.2533*np.std(age_all[w0])/np.sqrt(len(age_all[w0])) ) )
			f.write('\n')
			f.close()

		else:
			f = open(data_dir + 'sc_gmc_assoc_stats.txt', 'a')
			f.write(gal_name + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), 1.2533*np.std(age_all)/np.sqrt(len(age_all)) ) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), 1.2533*np.std(age_all[w1])/np.sqrt(len(age_all[w1])) ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), 1.2533*np.std(age_all[w2])/np.sqrt(len(age_all[w2])) ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), 1.2533*np.std(age_all[w3])/np.sqrt(len(age_all[w3])) ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), 1.2533*np.std(age_all[w0])/np.sqrt(len(age_all[w0])) ) )
			f.write('\n')

			f.close()


def generate_mega_df(galaxy_list, data_dir, output=True):

	gal_id 		= galaxy_list['id']
	gal_alt_id 	= galaxy_list['alt_id']
	gal_dist 	= galaxy_list['dist']
	

	# loop through all the galaxies in the list
	for i in range(len(galaxy_list)):

		gal_name = gal_id[i]

		# read in the csv
		df = pd.read_csv(data_dir + '%s/%s_cluster_catalog_in_mask_class12_assoc_gmc.csv'%(gal_name, gal_name))

		# need column for the galaxy name
		gal_name_col = np.array([gal_name for i in range(len(df))])

		df.insert(loc=0, column='gal_name', value=gal_name_col)


		if i == 0:
			# initialize the mega df using the first galaxy
			mega_df = df.copy()
		else:
			mega_df = mega_df.append(df)

	if output:
		mega_df.to_csv(data_dir + 'sc_gmc_assoc_mega.csv', index=False)

	return mega_df




if __name__ == '__main__':

	# set path of the data dir
	data_dir = '/cherokee1/turner/phangs/cf/data/'
	
	# read in list of galaxies 
	master_galaxy_list = ascii.read('master_galaxy.list')
	galaxy_list = ascii.read('galaxy.list')
	
	# generate_sc_gmc_assoc_df(galaxy_list, data_dir)

	# make_sc_gmc_assoc_hists(galaxy_list, data_dir)

	# mega_df = generate_mega_df(galaxy_list, data_dir, output=True)

	mega_df = pd.read_csv(data_dir + 'sc_gmc_assoc_mega.csv')
	
	"""	simple environmental masks cheatsheet
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
	wcenter       = mega_df['env_mask_val'] == 1
	wbar_idx      = np.where((mega_df['env_mask_val'] == 2) | (mega_df['env_mask_val'] == 3) ) 
	winterarm_idx = np.where((mega_df['env_mask_val'] == 4) | (mega_df['env_mask_val'] == 7) | (mega_df['env_mask_val'] == 8))
	wspiral_idx   = np.where((mega_df['env_mask_val'] == 5) | (mega_df['env_mask_val'] == 6))
	wdisk         = mega_df['env_mask_val'] == 9

	wall = [wcenter, wbar_idx, winterarm_idx, wspiral_idx, wdisk]
	names = ['center', 'bar', 'interarm', 'spiralarm', 'disk']

	for i in range(len(wall)):

		df = mega_df.loc[wall[i]]

		# sc_gmc_assoc_hist(df, filename=data_dir+'sc_gmc_assoc_hist_%s'%(names[i]))

		# star cluster ages
		age_all = df['age'].to_numpy()
		lage_all = np.log10(age_all)
		age_err_all = df['age_err'].to_numpy()
		lage_err_all = age_err_all/age_all/np.log(10)

		w0 = df['assoc_num'] == 0
		w1 = df['assoc_num'] == 1
		w2 = df['assoc_num'] == 2
		w3 = df['assoc_num'] == 3

		# log the stats for each env
		if i == 0:
			f = open(data_dir + 'sc_gmc_assoc_stats_env.txt', 'w')
			f.write(names[i] + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), 1.2533*np.std(age_all)/np.sqrt(len(age_all)) ) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), 1.2533*np.std(age_all[w1])/np.sqrt(len(age_all[w1])) ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), 1.2533*np.std(age_all[w2])/np.sqrt(len(age_all[w2])) ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), 1.2533*np.std(age_all[w3])/np.sqrt(len(age_all[w3])) ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), 1.2533*np.std(age_all[w0])/np.sqrt(len(age_all[w0])) ) )
			f.write('\n')
			f.close()
		else:
			f = open(data_dir + 'sc_gmc_assoc_stats_env.txt', 'a')
			f.write(names[i] + '\n')
			f.write('All star clusters median age: %.2f +/- %.2f Myr \n'%(np.median(age_all), 1.2533*np.std(age_all)/np.sqrt(len(age_all)) ) )
			f.write('Within 1 R_gmc median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w1]), 1.2533*np.std(age_all[w1])/np.sqrt(len(age_all[w1])) ) )
			f.write('1 < R_gmc <= 2 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w2]), 1.2533*np.std(age_all[w2])/np.sqrt(len(age_all[w2])) ) )
			f.write('2 < R_gmc <= 3 median age:    %.2f +/- %.2f Myr \n'%(np.median(age_all[w3]), 1.2533*np.std(age_all[w3])/np.sqrt(len(age_all[w3])) ) )
			f.write('Unassociated median age:      %.2f +/- %.2f Myr \n'%(np.median(age_all[w0]), 1.2533*np.std(age_all[w0])/np.sqrt(len(age_all[w0])) ) )
			f.write('\n')
			f.close()