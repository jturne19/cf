"""
utilities/functions for manipulating the star cluster catalogs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table

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



