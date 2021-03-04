## using sqlite3 to build a database using the class123 mega df
useful tutorial:
[https://www.tutorialspoint.com/sqlite/index.htm](https://www.tutorialspoint.com/sqlite/index.htm)

```bash
sudo apt install sqlite3

sqlite3 clustersDB.db
```

set to csv mode 
```
.mode csv
```

import the csv mega df to the 'clusters' table
```
.import /cherokee1/turner/phangs/cf/data/sc_gmc_assoc_mega.class123.csv clusters
```

verify import
```
.schema clusters
```

change back to column mode for easier reading
```
.mode column
```

thinks everything is text so change that
rename clusters table to clusters_old
and create new clusters table with the data type set for each column
```sql
ALTER TABLE clusters RENAME TO clusters_old;

CREATE TABLE clusters (
	gal_name TEXT NOT NULL,
	id INTEGER NOT NULL, 
	class INTEGER NOT NULL,
	ra REAL NOT NULL,
	dec REAL NOT NULL,
	hst_x REAL NOT NULL,
	hst_y REAL NOT NULL,
	alma_x REAL NOT NULL,
	alma_y REAL NOT NULL,
	age REAL NOT NULL,
	age_err REAL NOT NULL,
	mass REAL NOT NULL,
	mass_err REAL NOT NULL,
	ebv REAL NOT NULL,
	ebv_err REAL NOT NULL,
	env_mask_val INTEGER NOT NULL,
	nn_cloudnum INTEGER NOT NULL,
	nn_gmc_mlum REAL NOT NULL,
	nn_gmc_radius_pc REAL NOT NULL,
	nn_gmc_radius_asec REAL NOT NULL,
	nn_sep_asec REAL NOT NULL,
	nn_dist_pc REAL NOT NULL,
	assoc_gmc_cloudnum INTEGER NULL,
	assoc_gmc_mlum REAL NULL,
	assoc_gmc_radius_pc REAL NULL,
	assoc_num INTEGER NOT NULL
);

INSERT INTO clusters (gal_name,id,class,ra,dec,hst_x,hst_y,alma_x,alma_y,age,age_err,mass,mass_err,ebv,ebv_err,env_mask_val,nn_cloudnum,nn_gmc_mlum,nn_gmc_radius_pc,nn_gmc_radius_asec,nn_sep_asec,nn_dist_pc,assoc_gmc_cloudnum,assoc_gmc_mlum,assoc_gmc_radius_pc,assoc_num)
	SELECT gal_name,id,class,ra,dec,hst_x,hst_y,alma_x,alma_y,age,age_err,mass,mass_err,ebv,ebv_err,env_mask_val,nn_cloudnum,nn_gmc_mlum,nn_gmc_radius_pc,nn_gmc_radius_asec,nn_sep_asec,nn_dist_pc,assoc_gmc_cloudnum,assoc_gmc_mlum,assoc_gmc_radius_pc,assoc_num
	FROM clusters_old;
```

## couple of useful examples

show everything 
```
SELECT * FROM clusters;
```

pull out specific galaxies
```
SELECT * 
FROM clusters
WHERE gal_name = 'ngc0628';
```

select just class 1 and 2
```
SELECT * 
FROM clusters
WHERE class IN (1, 2);
```

select just class 3
```
SELECT * 
FROM clusters
WHERE class = 3;
```

select class 1 and 2 in ngc3351
```
SELECT * 
FROM clusters
WHERE class IN (1, 2) AND gal_name = 'ngc3351';
```

get ages for all spiral arm clusters (assoc_num = )
```
SELECT age
FROM clusters 
WHERE env_mask_val IN (5,6);
```

select multiple columns in certain galaxy 
```
SELECT age, age_err, mass, mass_err
FROM clusters
WHERE gal_name = 'ngc3351';
```

clusters which are associated with gmcs 
```
SELECT *
FROM clusters 
WHERE assoc_num IS NOT 0;
```

## calling sqlite from python
```python
import numpy as np
import sqlite3

conn = sqlite3.connect('clustersDB.db')
print('Opened database successfully')

cursor = conn.execute('SELECT age, age_err, mass, mass_err FROM clusters')

# return a list of all the rows returned by the command executed in cursor
sc = cursor.fetchall()
# convert to numpy array
sc = np.array(sc)

# close connection
conn.close()
```