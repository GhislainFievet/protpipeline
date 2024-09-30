# Using protpipeline

The protpipeline package permits to build proteome from maxquant files. It handles the following steps:
- data cleaning, and sample => renaming max_quant_process.r
- filtering of the data related to missing values => prot_filter.r
- detection of the MNAR proteins => mnar_filter.r
- normalization of the data => prot_normalize.r
- imputation of the missing data => prot_impute.r
- rename of the protein from uniprot IDs to HGNC IDs => prot_rename.r
- generation of plots for QC purpose => prot_plots.r

All steps are integrated into a pipeline and can be called from the function protPipeline().
