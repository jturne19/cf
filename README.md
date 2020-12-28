# cf
PHANGS correlation function work with ALMA GMC catalog and HST cluster catalog

browse the figures showing the location of star clusters and GMCs - [http://physics.uwyo.edu/~turner/research/correlation/plots/](http://physics.uwyo.edu/~turner/research/correlation/plots/)


## Data

### Current Galaxies
| Galaxy   |Number GMCs | Number Star Clusters | Number Class 1,2,3   | Number Class 1,2     | Notes                           |
|----------|------------|----------------------|----------------------|----------------------|---------------------------------|
| NGC0628  |    729     | 11354 (9033+2321)    |     789 (678+111)    |       580 (489+91)   | (Center + East pointings)       |
| NGC1365  |    951     |       3291           |         789          |          635         |                                 |
| NGC1433  |    177     |       2475           |         293          |          194         |                                 |
| NGC1559  |    727     |       8363           |         927          |          715         | Still v1.0                      |
| NGC1566  |    1111    |       9898           |         851          |          685         |                                 |
| NGC1672  |    480     |       8754           |                      |                      | Missing cluster classifications |
| NGC1792  |    388     |       4660           |         675          |          567         |                                 |
| NGC2775  |            |       628            |                      |                      | Missing GMC catalog & cluster classificaitons |
| NGC3351  |    315     |       4329           |         468          |          302         |                                 |
| NGC3627  |    1064    |       10673          |         958          |          774         |                                 |
| NGC4535  |    603     |       2167           |         452          |          345         | Still v1.0                      |
| NGC4548  |    207     |       788            |         271          |          195         |                                 |
| NGC4569  |    325     |       1309           |                      |                      | Missing cluster classifications |
| NGC4571  |    140     |       1085           |         262          |          162         |                                 |
| NGC4654  |            |       2812           |                      |                      | Missing GMC catalog & cluster classifications |
| NGC4689  |    473     |       1358           |                      |                      | Missing cluster classifications |
| NGC4826  |    169     |       1935           |                      |                      | Missing cluster classifications |
| NGC5248  |    642     |       3434           |                      |                      | Missing cluster classifications |
| Total    |    8,501   |       79,313         |         6,735        |          5,154       |                                 |


### cluster catalogs 

updated to the v1.1 internal catalogs

manually downloaded from `PHANGS-HST box > Cluster Detection > v1.1 Internal Catalogs`

`data_dir/ngc????/hst/ngc????_phangshst_base_catalog.fits`

PHANGS-HST images images downloaded from [PHANGS-HST STScI website](https://phangs.stsci.edu/)

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
|`utils.py`			   |`utils`    | contains all the functions/utilities 		             |					 
|`catalog_setup.py`	   |`workflow` | calling the functions for creating ds9 region files for everything and other 'setup' type things   |
|`sc_gmc_sep.py`	   |`workflow` | calling the functions for finding the star cluster - gmc nearest neighbors and making histograms   |


### process notes
- generate ds9 region files for all the GMCs using `all_galaxies_gmc_region_files`
- generate ds9 region files for all the star clusters using `all_galaxies_clust_region_files`
- generate ds9 region files for the HST footprints using FLC files with `all_galaxies_flc_region_files`
- manually made HST full footprint ds9 region files since i couldn't think of a way to do this automatically
- combined the NGC0628C and NGC0628E cluster catalogs to make NGC0628. already have HST maps combined so this makes sense
- generate the outline plots which overlay the GMC ellipses with the star cluster positions and HST footprints with `all_galaxies_outline_plots`
- generate masks which define the overlap of the HST and ALMA footprints with `generate_overlap_mask`. need to use uwpa with more memory for the larger hst images (628, 3351, 3627)
- find star cluster - gmc nearest neighbors and makes histograms split at 10 Myr in `sc_gmc_sep.py`