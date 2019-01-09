## `count_matrix.txt`

* First row: sample IDs
* First column: feature names
* Row name format (gene-level): `<gene_id>|<gene_name>|<gene_type>`
* Row name format (domain-level): `<gene_id>|<gene_name>|<gene_type>|<domain_id>`

## `feature_importances.txt`

* First row: None
* Column 1: feature name
* Column 2: feature importance
  
## `metrics.txt`

* First row: None
* Column 1: run index
* Column 2: test ROAUC
* Column 3: train ROAUC