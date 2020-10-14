"""
utilities/functions for manipulating the gmc catalogs

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky


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