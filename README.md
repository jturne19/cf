# cf
PHANGS correlation function work with ALMA GMC catalog and HST cluster catalog


## Data

### Current Galaxies
| Galaxy   |Number GMCs | Number Star Clusters | Number Class 1,2,3   | Number Class 1,2     | Notes                           |
|----------|------------|----------------------|----------------------|----------------------|---------------------------------|
| NGC0628C |            |       8501           |         661          |          476         |                                 |
| NGC0628E |            |       1997           |         111          |          91          |                                 |
| NGC1433  |            |       2155           |         287          |          188         |                                 |
| NGC1559  |            |       8363           |         927          |          715         |                                 |
| NGC1566  |            |       8679           |         835          |          672         |                                 |
| NGC1792  |            |       4295           |         615          |          521         |                                 |
| NGC3351  |            |       4310           |         468          |          302         |                                 |
| NGC3627  |            |       10153          |         941          |          761         |                                 |
| NGC4535  |            |       2167           |         452          |          345         |                                 |
| NGC4548  |            |       637            |                      |                      | Missing cluster classifications |
| NGC4569  |            |       1159           |                      |                      | Missing cluster classifications |
| NGC4571  |            |       832            |                      |                      | Missing cluster classifications |
| NGC4689  |            |       1358           |                      |                      | Missing cluster classifications |
| NGC5248  |            |       3212           |                      |                      | Missing cluster classifications |
|----------|------------|----------------------|----------------------|----------------------|---------------------------------|
| Total    |            |       57,818         |         5,297        |          4,071       |                                 |

### cluster catalogs 

currently using the v1.0 internal catalogs

manually downloaded from PHANGS-HST box > Cluster Detection > v1.0 Internal Catalogs 

`data_dir/hst/ngc????/ngc????_phangshst_base_catalog.fits`


### GMC catalogs

currently using the native resolution v3.4_ST1.5 GMC catalogs

manually downloaded from PHANGS google drive > scratch > gmccats > v3p4_ST1p5 > native

`data_dir/alma/ngc????/ngc????_12m+7m+tp_co21_native_props.fits`


## Scripts

| Name                 | Location  | Description |
|----------------------|-----------|-------------|
|`cluster_cat.py`      |`workflow` | create separate fits tables with class 1,2,3 and 1,2 for each galaxy  |
|`cluster_cat_utils.py`|`utils`    | contains functions for manipulating the cluster catalogs              |
