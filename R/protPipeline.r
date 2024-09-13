protPipeline <- function(output_dir, max_quant_dir, yaml_config_file) {
    # args <- commandArgs(trailingOnly = TRUE)

    # if (length(args) < 3 ){
    #     print("You need to give an output_dir, the max quant dir and a yaml config file")
    #     quit()
    # }

    # copy max_quant_dir folder to output_dir


    library(yaml)
    args = read_yaml(yaml_config_file)

    my_seed = args$seed
    conditions_path = args$conditions_path
    thresholds_path = args$thresholds_path
    group_threshold_mode = args$group_threshold_mode
    group_threshold = args$group_threshold
    global_threshold = args$global_threshold
    MNAR_filter = args$MNAR_filter
    MNAR_threshold = args$MNAR_threshold
    partial_imputation = args$partial_imputation

    proteinGroups_path = file.path(max_quant_dir, args$proteinGroups_path)
    peptides_path = file.path(max_quant_dir, args$peptides_path)
    prot_rename_db_path = args$prot_rename_db_path

    # pipeline_mode
    pipeline_mode = args$pipeline_mode

    # filter = args$filter
    atLeastOne = args$atLeastOne
    remove_reverse_identified_contaminant = args$remove_reverse_identified_contaminant
    peptide_occurence_filter = args$peptide_occurence_filter
    imputation_method = args$imputation_method

    # lib path
    # max_quant_process_path = args$max_quant_process_path
    # prot_filter_path = args$prot_filter_path
    # prot_normalize_path = args$prot_normalize_path
    # prot_impute_path = args$prot_impute_path
    # rf_impute_lib_path = args$rf_impute_lib_path
    # prot_plots_path = args$prot_plots_path

    # output directories path
    output_dir_new_names = args$output_dir_new_names
    output_dir_filtered = args$output_dir_filtered
    output_dir_normalized = args$output_dir_normalized
    output_dir_imputed = args$output_dir_imputed
    output_dir_prot_rename = args$output_dir_prot_rename
    output_dir_plots = args$output_dir_plots

    # output files
    conditions_new_sample_names = args$conditions_new_sample_names
    prot_annotations_path = args$prot_annotations_path
    pep_annotations_path = args$pep_annotations_path
    output_proteome_data = args$output_proteome_data
    output_peptidome_data = args$output_peptidome_data
    output_both_data = args$output_both_data

    # Normalisation method
    norm_method = args$norm_method


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
    # source(max_quant_process_path)
    # source(prot_filter_path)
    # source(prot_normalize_path)
    # source(prot_impute_path)
    # source(prot_plots_path)

    ###### Create directories if needed ######

    if (!dir.exists(output_dir)){
        dir.create(output_dir, recursive = TRUE)
    }

    # Get list of folders in output_dir
    output_dir_list <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)

    # Get date and time (with seconds) for new folder
    my_prefix <- format(Sys.time(), "%Y%m%d_%H%M%S")
    new_prot_dir <- paste0("pipe_results", "_", my_prefix)
    output_dir <- file.path(output_dir, new_prot_dir)
    dir.create(output_dir, recursive = TRUE)
    dir.create(file.path(output_dir, "max_quant_data"), recursive = TRUE)
    copy_files <- function(path, dir1, dir2) {
        file.copy(file.path(dir1, path), file.path(dir2, path))
    }
    files <- list.files(max_quant_dir, recursive = TRUE, full.names = FALSE)
    lapply(files, function(x){
        copy_files(x, max_quant_dir, file.path(output_dir, "max_quant_data"))
    })
    if (file.exists(prot_rename_db_path)){
        file.copy(prot_rename_db_path, file.path(output_dir, "prot_rename_db.txt"))
    }

    # Copy files to the result folder output
    file.copy(yaml_config_file, file.path(output_dir, paste0(my_prefix, "_config.yaml")))
    file.copy(conditions_path, file.path(output_dir, paste0(my_prefix, "_conditions.txt")))

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
    if (!dir.exists(file.path(output_dir,output_dir_prot_rename))){
        dir.create(file.path(output_dir,output_dir_prot_rename), recursive = TRUE)
    }

    ###### Only keep LFQ.intensity.XXX samples needed ######
    # For proteins
    if ( pipeline_mode == "prot"){
        data_new_sample_names_path <- file.path(output_dir,output_dir_new_names,
                paste0(my_prefix, "_", output_proteome_data))
        conditions_new_sample_names_path <- file.path(output_dir,
                paste0(my_prefix, "_", conditions_new_sample_names))
        prot_annotations_path <- file.path(output_dir,
                paste0(my_prefix, "_", prot_annotations_path))
        df_prot_results <- prot_max_quant_process( proteinGroups_path, conditions_path,
                                data_new_sample_names_path, conditions_new_sample_names_path,
                                prot_annotations_path, remove_reverse_identified_contaminant,
                                my_prefix)
    }
    # For peptides
    if ( pipeline_mode == "pep"){
        pep_data_new_sample_names_path <- file.path(output_dir,output_dir_new_names,
                paste0(my_prefix, "_", output_peptidome_data))
        conditions_new_sample_names_path <- file.path(output_dir,
                paste0(my_prefix, "_", conditions_new_sample_names))
        pep_annotations_path <- file.path(output_dir,
                paste0(my_prefix, "_", pep_annotations_path))
        df_pep_results <- pep_max_quant_process( peptides_path, conditions_path, pep_data_new_sample_names_path,
                                conditions_new_sample_names_path, pep_annotations_path,
                                remove_reverse_identified_contaminant, my_prefix)
    }


    ###### Filter proteins ######
    # For proteins
    if ( pipeline_mode == "prot"){
        data_filtered_path <- file.path(output_dir,output_dir_filtered, 
                    paste0(my_prefix, "_", output_proteome_data))
        message(data_new_sample_names_path, conditions_new_sample_names_path,
                    thresholds_path, data_filtered_path, atLeastOne, my_prefix)
        df_prot_filtered <- prot_filter(data_new_sample_names_path,
                    conditions_new_sample_names_path, thresholds_path,
                    data_filtered_path, atLeastOne, my_prefix, group_threshold_mode,
                    group_threshold, global_threshold)
    }

    # For peptides
    if ( pipeline_mode == "pep"){
        pep_data_filtered_path <- file.path(output_dir,output_dir_filtered,
                    paste0(my_prefix, "_", output_peptidome_data))
        message(pep_data_new_sample_names_path, conditions_new_sample_names_path,
                    thresholds_path, data_filtered_path, atLeastOne, my_prefix)
        df_pep_filtered <- pep_filter(pep_data_new_sample_names_path,
                    conditions_new_sample_names_path, thresholds_path,
                    pep_annotations_path, pep_data_filtered_path,
                    peptide_occurence_filter, atLeastOne, my_prefix,
                    group_threshold_mode, group_threshold, global_threshold)
    }

    if ( pipeline_mode == "prot" && MNAR_filter){
        prot_mnar_filtered_path <- file.path(output_dir, output_dir_filtered,
                    paste0(my_prefix, "_prot_mnar.txt"))
        message(data_new_sample_names_path,
                    conditions_new_sample_names_path,
                    data_filtered_path, MNAR_threshold)
        df_prot_mnar_filtered <- prot_mnar_filter(data_new_sample_names_path,
                    conditions_new_sample_names_path,
                    prot_mnar_filtered_path, MNAR_threshold, data_filtered_path)
    }
    # For peptides
    if ( pipeline_mode == "pep" && MNAR_filter){
        pep_data_filtered_path <- file.path(output_dir,output_dir_filtered,
                    paste0(my_prefix, "_pep_mnar.txt"))
        message(pep_data_new_sample_names_path, conditions_new_sample_names_path,
                    thresholds_path, data_filtered_path, atLeastOne, my_prefix)
        df_pep_filtered <- pep_mnar_filter(pep_data_new_sample_names_path,
                    conditions_new_sample_names_path, 
                    pep_annotations_path, pep_data_filtered_path,
                    peptide_occurence_filter, MNAR_threshold)
    }

    ###### Normalization ######
    filter_folder_path <- file.path(output_dir, output_dir_filtered)
    norm_output_path <- file.path(output_dir, output_dir_normalized,
                paste0(my_prefix, "_prot_pep.txt"))
    normalize_filter_folder(filter_folder_path, norm_output_path, norm_method,
                output_proteome_data, output_peptidome_data, my_prefix,
                output_dir, output_dir_normalized)
    
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
    impute_output_path <- file.path(output_dir, output_dir_imputed,
                paste0(my_prefix, "_prot_pep.txt"))
    partial_impute_output_path <- file.path(output_dir, output_dir_imputed,
                paste0(my_prefix, "_partialImpute_prot_pep.txt"))
    prot_impute <- prot_impute(norm_output_path, impute_output_path,
            partial_impute_output_path,
            imputation_method, k=k, rowmax=rowmax, colmax=colmax, MNAR_filter,
            thresholds_path, group_threshold_mode, group_threshold,
            my_seed, conditions_path)

    ###### Convert uniprot id to hgnc symbol ######
    rename_output_path <- file.path(output_dir, output_dir_prot_rename,
                paste0(my_prefix, "_prot_pep_renamed.txt"))
    partial_rename_output_path <- file.path(output_dir, output_dir_prot_rename,
                paste0(my_prefix, "_partialImpute_prot_pep_renamed.txt"))
    prot_rename(impute_output_path, partial_impute_output_path, rename_output_path, partial_rename_output_path, prot_rename_db_path)

    ###### Plots ######
    # For proteins
    if (pipeline_mode == "prot" ){
        print("PLOT   pipeline_mode == prot")
        before_filter_path <- file.path(output_dir, output_dir_new_names,
                    paste0(my_prefix, "_", output_proteome_data))
        after_filter_path <- file.path(output_dir, output_dir_filtered,
                    paste0(my_prefix, "_", output_proteome_data))
        after_normalization_path <- file.path(output_dir, output_dir_normalized,
                    paste0(my_prefix, "_", output_proteome_data))
        after_imputation_path <- file.path(output_dir, output_dir_imputed,
                    paste0(my_prefix, "_", output_proteome_data))
        plot_dir <- file.path(output_dir, output_dir_plots)
        prot_plots(before_filter_path, after_filter_path, after_normalization_path,
                after_imputation_path, conditions_new_sample_names_path, plot_dir, "prot", my_prefix)
    }

    # For proteins
    if (pipeline_mode == "pep" ){
        print("PLOT   pipeline_mode == pep")
        before_filter_path <- file.path(output_dir, output_dir_new_names,
                    paste0(my_prefix, "_", output_peptidome_data))
        after_filter_path <- file.path(output_dir, output_dir_filtered,
                    paste0(my_prefix, "_", output_peptidome_data))
        after_normalization_path <- file.path(output_dir, output_dir_normalized,
                    paste0(my_prefix, "_", output_peptidome_data))
        after_imputation_path <- file.path(output_dir, output_dir_imputed,
                    paste0(my_prefix, "_", output_peptidome_data))
        plot_dir <- file.path(output_dir, output_dir_plots)
        prot_plots(before_filter_path, after_filter_path, after_normalization_path, after_imputation_path,
                conditions_new_sample_names_path, plot_dir, "pep", my_prefix)
    }
}