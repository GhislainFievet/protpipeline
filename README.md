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

![Pipeline](https://github.com/GhislainFievet/protpipeline/blob/main/im/pipeline.png?raw=true)

## Installation

You can install `protPipeline` directly from GitHub using `devtools`:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install protPipeline from GitHub
devtools::install_github("GhislainFievet/protpipeline")
```

## Usage

```r
# output_dir is the output pipeline folder 
# max_quant_dir is the math of the folder containing the maxquant files
# yaml_config_file is the path of the YAML configuration file
protPipeline(output_dir, max_quant_dir, yaml_config_file, dir_name="my_dir_name")
```

Toy example:
```r
# MaxQuant folder
max_quant_dir <- system.file("extdata", "max_quant_toy_example", package="protpipeline")
# Conditions file: file containing samples id, the studied conditions, and the new sample names
toy_conditions_path <- system.file("extdata", "conditions_toy.txt", package="protpipeline")
# .yaml config file, edit and modify conditions_path variable to `toy_conditions_path`
yaml_config_file <- system.file("extdata", "config_toy_example.yaml", package="protpipeline")

# Run the pipeline
protPipeline("output_dir", max_quant_dir, yaml_config_file, dir_name="toy_example")
```


## Parameter file

```yaml
# Seed for reproductibility
seed: 123

# A file associating
# - a condition to each sample
# - a new sample id to each sample
#
# 2 available formats for conditions_path file
# format 1: tab separated
# label  condition   replicate
# NCY38 CR      1
# NCY54 CR      2
# NCY70 CR      3
# NCY200        CR      4
# AVI01 R       1
# CLM01 R       2
# CLM07 R       3
#
# format 2: tab separated
# label condition       replicate       new_name
# NCY38 CR      1       newname_1
# NCY54 CR      2       newname_2
# NCY70 CR      3       newname_3
# AVI01 R       1       newname_4
# CLM01 R       2       newname_5
# CLM07 R       3       newname_6
conditions_path : "data/conditions.txt"

# "prot" to construct proteome, "pep" to construct peptidome
pipeline_mode: "prot"
# Number of occurences of the same protein for peptides to be kept
peptide_occurence_filter: 2
# Name of directories to create
output_dir_new_names: "01_new_names"
output_dir_filtered: "02_filtered"
output_dir_normalized: "03_normalized"
output_dir_imputed: "04_imputed"
output_dir_prot_rename: "05_prot_rename"
output_dir_plots: "plots"
# Name of files to create
conditions_new_sample_names: "conditions_new_sample_names.txt"
prot_annotations_path: "prot_annotations.txt"
pep_annotations_path: "pep_annotations.txt"
output_proteome_data: "proteome_data.txt"
output_peptidome_data: "peptidome_data.txt"
output_both_data: "prot_pep_data.txt"
# proteinGroups_path and peptides_path are files from Max Quant
proteinGroups_path : "proteinGroups.txt"
peptides_path : "peptides.txt"

# Step 1: clean and new names
# 
# Remove contaminants
remove_reverse_identified_contaminant: TRUE

# Step 2: filtering and MNARs
#
# 4 available group_threshold_mode:
# "file" (look for thresholds in a file)
# "by_group" (use group_threshold on groups)
# "simple" (use global_threshold on all samples)
# "simple+group" (use group_threshold on at least 1 group and global_threshold on all samples)
group_threshold_mode : "simple+group"
# format for thresholds_path file
# tab separated
# R   CR   CS   S
# 16  7 15  5
thresholds_path : "data/thresholds.txt"
group_threshold : 0.6667
global_threshold: 0.5
MNAR_filter: TRUE
MNAR_threshold : "1/3"
atLeastOne : FALSE

# Step 3: normalization
#
# Normalisation methods: "none", "vsn", "quantiles"
norm_method: "quantiles"

# Step 4: imputation
#
# Imputation_method methods: "rf", "knn", "none"
imputation_method: "rf"
# if partial_imputation == TRUE, then it imputes the MAR proteins and the MNAR groups containing at most partial_imputation_threshold missing.
# if partial_imputation == FALSE, then it imputes only the MAR proteins.
partial_imputation: TRUE
partial_imputation_threshold: 0.3334
# knn.impute parameters for "knn" imputation method
# knn_k: 10
# knn_rowmax: 0.4
# knn_colmax: 0.8
knn_k: 10
knn_rowmax: 0.4
knn_colmax: 0.8
# knn.impute parameters for "rf" imputation method
# rf_knn_k: 10
# rf_knn_rowmax: 0.4
# rf_knn_colmax: 0.8
rf_knn_k: 10
rf_knn_rowmax: 0.4
rf_knn_colmax: 0.8

# Step 5: convert uniprot IDs to HGNC IDs
#
# rename methods can be "file", "biomart", or "none"
rename_method: "file"
# protein id should be: Protein.IDs column
# hgnc should be: HGNC.symbol_result column
# prot_rename_db_path : "data/FIXED_ProteinDB_annots4127.txt"
prot_rename_db_path : "data/simple_rename_db.txt"

# Step 6: QC plots

```
