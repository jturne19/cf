# cf
PHANGS correlation function work with ALMA GMC catalog and HST cluster catalog

browse the figures showing the location of star clusters and GMCs - [http://physics.uwyo.edu/~turner/research/correlation/plots/](http://physics.uwyo.edu/~turner/research/correlation/plots/)


## Data

### Current Galaxies
| Galaxy   |Number GMCs | Number Star Clusters | Number Class 1,2,3   | Number Class 1,2     | Notes                           |
|----------|------------|----------------------|----------------------|----------------------|---------------------------------|
| NGC0628  |    729     | 10498 (8501+1997)    |     772 (661+111)    |       567 (476+91)   | (Center + East pointings)       |
| NGC1433  |    177     |       2155           |         287          |          188         |                                 |
| NGC1559  |    727     |       8363           |         927          |          715         |                                 |
| NGC1566  |    1111    |       8679           |         835          |          672         |                                 |
| NGC1792  |    388     |       4295           |         615          |          521         |                                 |
| NGC3351  |    315     |       4310           |         468          |          302         |                                 |
| NGC3627  |    1064    |       10153          |         941          |          761         |                                 |
| NGC4535  |    603     |       2167           |         452          |          345         |                                 |
| NGC4548  |    207     |       637            |                      |                      | Missing cluster classifications |
| NGC4569  |    325     |       1159           |                      |                      | Missing cluster classifications |
| NGC4571  |    140     |       832            |                      |                      | Missing cluster classifications |
| NGC4689  |    473     |       1358           |                      |                      | Missing cluster classifications |
| NGC5248  |    642     |       3212           |                      |                      | Missing cluster classifications |
| Total    |    6,901   |       57,818         |         5,297        |          4,071       |                                 |


### cluster catalogs 

currently using the v1.0 internal catalogs

manually downloaded from `PHANGS-HST box > Cluster Detection > v1.0 Internal Catalogs`

`data_dir/ngc????/hst/ngc????_phangshst_base_catalog.fits`

PHANGS-HST images images downloaded from [PHANGS-HST STScI website](www.phangs.stsci.edu)

`data_dir/ngc????/hst/ngc????_uvis_f275w_exp_drc_sci.fits`


### GMC catalogs

currently using the native resolution v3.4_ST1.5 GMC catalogs

manually downloaded from `PHANGS google drive > scratch > gmccats > v3p4_ST1p5 > native`

`data_dir/ngc????/alma/ngc????_12m+7m+tp_co21_native_props.fits`

PHANGS-ALMA moment 0 broad maps also downloaded from `PHANGS google drive > Archive > ALMA > delivery_v3p4 > broad_maps`

`data_dir/ngc????/alma/ngc????_12m+7m+tp_co21_broad_mom0.fits`

## Scripts

| Name                 | Location  | Description |
|----------------------|-----------|-------------|
|`cluster_cat_utils.py`|`utils`    | contains functions for manipulating the cluster catalogs              |
|`gmc_cat_utils.py`    |`utils`    | contains functions for manipulating the gmc catalogs 				   |
|`image_utils.py`	   |`utils`    | contains functions for manipulating the images   	 				   |
|`plot_utils.py` 	   |`utils`    | contains functions for creating figures/plots 						   |
|`catalog_setup.py`	   |`workflow` | calling the functions for creating ds9 region files for everything and other 'setup' type things   |


### process notes
- generate ds9 region files for all the GMCs using `all_galaxies_gmc_region_files`
- generate ds9 region files for all the star clusters using `all_galaxies_clust_region_files`
- generate ds9 region files for the HST footprints using FLC files with `all_galaxies_flc_region_files`
- manually made HST full footprint ds9 region files since i couldn't think of a way to do this automatically
- combined the NGC0628C and NGC0628E cluster catalogs to make NGC0628. already have HST maps combined so this makes sense
- generate the outline plots which overlay the GMC ellipses with the star cluster positions and HST footprints with `all_galaxies_outline_plots`
