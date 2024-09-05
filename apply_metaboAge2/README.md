# MetaboAge 2.0
This subfolder contains the scripts and data necessary to reproduce the MetaboAge 2.0 model, as described in the paper "Technical Report: A Comprehensive Comparison Between Quantifications of Nightingale Health's 1H-NMR Metabolomics Platform."
You can access the paper on https://pubmed.ncbi.nlm.nih.gov/38132863/.

## Applying the MetaboAge 2.0 Model
To apply the MetaboAge 2.0 model to your data, follow these steps:

Prepare Your Data:
1) Ensure to save your raw Nightingale data as a CSV file with the metabolites value (not the xls file).
	The CSV should have:
		Row names: Subject IDs.
		Column names: Metabolite names.
2) Open the script apply_metaboAge2.R.
3) Replace the placeholder metabo_mat with the path to your CSV file.
4) Run the script: Use the apply_metaboAge2.R script to apply the model.

By following these steps, you will be able to apply the MetaboAge 2.0 model to your data and reproduce the results as described in the publication.


## Authors
Daniele Bizzarri, 
Erik van den Akker