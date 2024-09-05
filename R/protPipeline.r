.onLoad <- function(libname, pkgname) {
    # Vérifie si BiocManager est installé
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    # Installe les packages Bioconductor nécessaires
    BiocManager::install("DEP")
    BiocManager::install("SummarizedExperiment")
    BiocManager::install("preprocessCore")
    BiocManager::install("imputeGenomicAlignments")
}


protPipeline <- function(output_dir, max_quant_dir, yaml_config_file) {
    # args <- commandArgs(trailingOnly = TRUE)

    # if (length(args) < 3 ){
    #     print("You need to give an output_dir, the max quant dir and a yaml config file")
    #     quit()
    # }

    library(yaml)
    args = read_yaml(yaml_config_file)

    conditions_path = args$conditions_path
    thresholds_path = args$thresholds_path
    proteinGroups_path = file.path(max_quant_dir, args$proteinGroups_path)
    peptides_path = file.path(max_quant_dir, args$peptides_path)

    # pipeline_mode
    pipeline_mode = args$pipeline_mode

    filter = args$filter
    atLeastOne = args$atLeastOne
    remove_reverse_identified_contaminant = args$remove_reverse_identified_contaminant
    peptide_occurence_filter = args$peptide_occurence_filter
    imputation_method = args$imputation_method

    # lib path
    max_quant_process_path = args$max_quant_process_path
    prot_filter_path = args$prot_filter_path
    prot_normalize_path = args$prot_normalize_path
    prot_impute_path = args$prot_impute_path
    rf_impute_lib_path = args$rf_impute_lib_path
    prot_plots_path = args$prot_plots_path

    # output directories path
    output_dir_new_names = args$output_dir_new_names
    output_dir_filtered = args$output_dir_filtered
    output_dir_normalized = args$output_dir_normalized
    output_dir_imputed = args$output_dir_imputed
    output_dir_plots = args$output_dir_plots

    # output files
    conditions_new_sample_names = args$conditions_new_sample_names
    prot_annotations_path = args$prot_annotations_path
    pep_annotations_path = args$pep_annotations_path
    output_proteome_data = args$output_proteome_data
    output_peptidome_data = args$output_peptidome_data
    output_both_data = args$output_both_data

    # Normalisation method
    prot_norm_method = args$prot_norm_method
    pep_norm_method = args$pep_norm_method
    both_norm_method = args$both_norm_method


    # knn.impute parameters for "knn" imputation method
    knn_k = args$knn_k
    knn_rowmax = args$knn_rowmax
    knn_colmax = args$knn_colmax

    # knn.impute parameters for "rf" imputation method
    rf_knn_k = args$rf_knn_k
    rf_knn_rowmax = args$rf_knn_rowmax
    rf_knn_colmax = args$rf_knn_colmax


    ###### Load libraries ######

    library("DEP") # for Proteomics analysis
    library("dplyr") # for data manipulation
    library("proteus") # for normalisation function here and is suitable for
    library("SummarizedExperiment")
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(preprocessCore)
    source(max_quant_process_path)
    source(prot_filter_path)
    source(prot_normalize_path)
    source(prot_impute_path)
    source(prot_plots_path)

    ###### Create directories if needed ######

    if (!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_filtered))){
        dir.create(file.path(output_dir,output_dir_filtered), recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_normalized))){
        dir.create(file.path(output_dir,output_dir_normalized), recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_imputed))){
        dir.create(file.path(output_dir,output_dir_imputed), recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_imputed))){
        dir.create(file.path(output_dir,output_dir_imputed), recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_plots))){
        dir.create(file.path(output_dir,output_dir_plots), recursive = TRUE)
    }
    if (!dir.exists(file.path(output_dir,output_dir_new_names))){
        dir.create(file.path(output_dir,output_dir_new_names), recursive = TRUE)
    }

    ###### Only keep LFQ.intensity.XXX samples needed ######
    # For proteins
    if ( pipeline_mode == "both" || pipeline_mode == "prot"){
        data_new_sample_names_path <- file.path(output_dir,output_dir_new_names,output_proteome_data)
        conditions_new_sample_names_path <- file.path(output_dir,conditions_new_sample_names)
        prot_annotations_path <- file.path(output_dir,prot_annotations_path)
        df_prot_results <- prot_max_quant_process( proteinGroups_path, conditions_path, data_new_sample_names_path, conditions_new_sample_names_path, prot_annotations_path, remove_reverse_identified_contaminant)
    }
    # For peptides
    if ( pipeline_mode == "both" || pipeline_mode == "pep"){
        pep_data_new_sample_names_path <- file.path(output_dir,output_dir_new_names,output_peptidome_data)
        conditions_new_sample_names_path <- file.path(output_dir,conditions_new_sample_names)
        pep_annotations_path <- file.path(output_dir,pep_annotations_path)
        df_pep_results <- pep_max_quant_process( peptides_path, conditions_path, pep_data_new_sample_names_path, conditions_new_sample_names_path, pep_annotations_path, remove_reverse_identified_contaminant )
    }


    ###### Filter proteins ######
    # For proteins
    if ( pipeline_mode == "both" || pipeline_mode == "prot"){
        data_filtered_path <- file.path(output_dir,output_dir_filtered, output_proteome_data)
        message(data_new_sample_names_path, conditions_new_sample_names_path, thresholds_path, data_filtered_path, atLeastOne)
        df_prot_filtered <- prot_filter(data_new_sample_names_path, conditions_new_sample_names_path, thresholds_path, data_filtered_path, atLeastOne)
    }

    # For peptides
    if ( pipeline_mode == "pep"){
        pep_data_filtered_path <- file.path(output_dir,output_dir_filtered, output_peptidome_data)
        message(pep_data_new_sample_names_path, conditions_new_sample_names_path, thresholds_path, data_filtered_path, atLeastOne)
        df_pep_filtered <- pep_filter(pep_data_new_sample_names_path, conditions_new_sample_names_path, thresholds_path, pep_annotations_path, pep_data_filtered_path, peptide_occurence_filter, atLeastOne)
    }
    if ( pipeline_mode == "both"){
        pep_data_filtered_path <- file.path(output_dir,output_dir_filtered, output_peptidome_data)
        message(pep_data_new_sample_names_path, conditions_new_sample_names_path, thresholds_path, data_filtered_path, atLeastOne)
        df_pep_filtered <- pep_filter_both(pep_data_new_sample_names_path, data_filtered_path, conditions_new_sample_names_path, thresholds_path, pep_annotations_path, prot_annotations_path, pep_data_filtered_path, peptide_occurence_filter, atLeastOne)
    }


    ###### Normalization ######
    # For proteins
    if ( pipeline_mode == "both" || pipeline_mode == "prot"){
        norm_output_path <- file.path(output_dir,output_dir_normalized, output_proteome_data)
        df_prot_normalized <- prot_normalize(data_filtered_path, norm_output_path, conditions_new_sample_names_path, prot_norm_method)
    }
    # For peptides
    if ( pipeline_mode == "both" || pipeline_mode == "pep"){
        pep_norm_output_path <- file.path(output_dir,output_dir_normalized, output_peptidome_data)
        df_pep_normalized <- prot_normalize(pep_data_filtered_path, pep_norm_output_path, conditions_new_sample_names_path, pep_norm_method)
    }


    ###### Imputation of missing values ######
    # determine parameters for knn impute if needed
    k <- 0
    rowmax <- 0
    colmax <- 0

    if ( imputation_method == "knn"){
        k <- knn_k
        rowmax <- knn_rowmax
        colmax <- knn_colmax
    }

    if ( imputation_method == "rf"){
        k <- rf_knn_k
        rowmax <- rf_knn_rowmax
        colmax <- rf_knn_colmax
    }

    # For proteins
    if ( pipeline_mode == "prot" ){
        impute_output_path <- file.path(output_dir, output_dir_imputed, output_proteome_data)
        prot_impute <- prot_impute(norm_output_path, impute_output_path, imputation_method, k=k, rowmax=rowmax, colmax=colmax)
    }
    # For peptides
    if ( pipeline_mode == "pep"){
        pep_impute_output_path <- file.path(output_dir, output_dir_imputed, output_peptidome_data)
        pep_impute <- prot_impute(pep_norm_output_path, pep_impute_output_path, imputation_method, k=k, rowmax=rowmax, colmax=colmax)
    }
    # For proteins + peptides
    if ( pipeline_mode == "both" ){
        impute_both <- concate_normalize_and_impute(data_filtered_path, pep_data_filtered_path, prot_annotations_path, pep_annotations_path, output_dir, output_both_data, conditions_new_sample_names_path, both_norm_method, imputation_method, k=k, rowmax=rowmax, colmax=colmax)
    }


    ###### Plots ######
    # For proteins
    if ( pipeline_mode == "both" || pipeline_mode == "prot" ){
        before_filter_path <- file.path(output_dir, output_dir_new_names, output_proteome_data)
        after_filter_path <- file.path(output_dir, output_dir_filtered, output_proteome_data)
        after_normalization_path <- file.path(output_dir, output_dir_normalized, output_proteome_data)
        after_imputation_path <- file.path(output_dir, output_dir_imputed, output_proteome_data)
        plot_dir <- file.path(output_dir, output_dir_plots)
        prot_plots(before_filter_path, after_filter_path, after_normalization_path, after_imputation_path, conditions_new_sample_names_path, plot_dir, "prot")
    }

    # For proteins
    if ( pipeline_mode == "both" || pipeline_mode == "pep" ){
        before_filter_path <- file.path(output_dir, output_dir_new_names, output_peptidome_data)
        after_filter_path <- file.path(output_dir, output_dir_filtered, output_peptidome_data)
        after_normalization_path <- file.path(output_dir, output_dir_normalized, output_peptidome_data)
        after_imputation_path <- file.path(output_dir, output_dir_imputed, output_peptidome_data)
        plot_dir <- file.path(output_dir, output_dir_plots)
        prot_plots(before_filter_path, after_filter_path, after_normalization_path, after_imputation_path, conditions_new_sample_names_path, plot_dir, "pep")
    }

    if ( pipeline_mode == "both" ){
        after_filter_path = file.path(output_dir, "data_before_normalization.txt")
        after_normalization_path <- file.path(output_dir, "data_after_normalization.txt")
        after_imputation_path <- file.path(output_dir, output_both_data)
        concat_plots(after_filter_path, after_normalization_path, after_imputation_path, conditions_new_sample_names_path, plot_dir, "both")
    }

}