# source("rf_impute.r")

prot_impute <- function (prot_path, output_path, output_partial_path, imputation_method, k=10, rowmax=0.4,
            colmax=0.8, MNAR_filter, partial_imputation_threshold, my_seed, conditions_path){
    message("")
    message(paste("Imputation method:", imputation_method) )
    if ( imputation_method=="knn" || imputation_method=="rf" ){
        message(" Parameters:")
        message("  - k: ",k)
        message("  - rowmax: ", rowmax)
        message("  - colmax: ", colmax)
        message("")
    }
    
    if ( partial_imputation_threshold == "1/3"){
        partial_imputation_threshold = 1/3
    }

    mat_norm = read.csv(prot_path, sep="\t")

    # Remove proteins with more than colmax missing values
    df_norm_quant_nof = as.data.frame(t(mat_norm))
    df_norm_quant <- df_norm_quant_nof[ lapply( df_norm_quant_nof, function(x) sum(is.na(x)) / length(x) ) < colmax ]

    df_result <- 0

    # Imputation with random forest
    if (imputation_method == "rf"){
        set.seed(my_seed)
        df_result <- rf_imput(df_norm_quant[,grep("mnar", colnames(df_norm_quant), invert=T)])$data
        set.seed(my_seed)
        df_result_full <- rf_imput(df_norm_quant)$data
    }

    df_conditions = read.csv(conditions_path, sep="\t")

    df_thresholds = data.frame(matrix(ncol=length(unique(df_conditions$condition)),nrow=1))
    colnames(df_thresholds) = unique(df_conditions$condition)
    for (i in 1:length(unique(df_conditions$condition))){
        df_thresholds[1,i] = sum(df_conditions$condition == unique(df_conditions$condition)[i])*partial_imputation_threshold
    }

    # Add removed proteins
    rows2add <- setdiff(colnames(df_norm_quant_nof), colnames(df_norm_quant))

    # Imputation with knn
    if (imputation_method == "knn"){
        library(impute)
        set.seed(my_seed)
        mat <- as.matrix(df_norm_quant[,grep("mnar", colnames(df_norm_quant), invert=T)])
        mat <- ifelse( mat == 'NA', NA, mat )
        mat <- ifelse( mat == '', NA, mat )
        message("Beginning of knn imputation")
        imputed_mat <- impute.knn(mat, k=k, rowmax=rowmax, colmax=colmax)
        message("End of knn imputation")
        df_knn_imput <- as.data.frame(imputed_mat$data)
        df_result <- df_knn_imput

        mat <- as.matrix(df_norm_quant)
        mat <- ifelse( mat == 'NA', NA, mat )
        mat <- ifelse( mat == '', NA, mat )
        message("Beginning of full knn imputation")
        imputed_mat <- impute.knn(mat, k=k, rowmax=rowmax, colmax=colmax)
        message("End of full knn imputation")
        df_knn_imput <- as.data.frame(imputed_mat$data)
        df_result_full <- df_knn_imput
    }

    for ( cond in unique(df_conditions$condition) ){
        cond_mask <- apply(t(df_norm_quant), 1, function(x){
            sum(!is.na(unname(unlist(x[df_conditions[df_conditions$condition == cond, "new_name"]])))) < df_thresholds[1,cond]
        })
        df_result_full[df_conditions[df_conditions$condition == cond, "new_name"], cond_mask] <- 
            df_norm_quant[df_conditions[df_conditions$condition == cond, "new_name"], cond_mask]
    }
    df_result_full[, grep("//prot//", colnames(df_result_full), value=T)] <-
        df_result[, grep("//prot//", colnames(df_result_full), value=T)]
    df_result_full[, grep("//pep//", colnames(df_result_full), value=T)] <-
        df_result[, grep("//pep//", colnames(df_result_full), value=T)]


    df_result <- cbind(df_result, df_norm_quant_nof[,rows2add])
    df_result_full <- cbind(df_result_full, df_norm_quant_nof[,rows2add])


    message(paste0("Writing ", output_path, " file"))
    write.table(t(df_result), file=output_path, sep="\t", append=F, quote=F)
    write.table(t(df_result_full), file=output_partial_path, sep="\t", append=F, quote=F)

    return(t(df_result))
}

# The concate_normalize_and_impute function takes as input the proteome filtered data,
# the peptidome filtered data, the proteome annotation data, and the peptidome annotation data.
# 1) It harmonizes and concat the proteome annotation data, and the peptidome annotation data into "annotations.txt" file.
# 2) It concatenates both data files into "data_before_normalization.txt"
# 3) It normalize it into "data_after_normalization.txt"
# 4) It imputes it into output_both_data file
concate_normalize_and_impute <- function(prot_data_path, pep_data_path, prot_annotations_path,
            pep_annotations_path, output_dir, output_both_data, conditions_new_sample_names_path, norm_method,
            imputation_method, k=10, rowmax=0.4, colmax=0.8, my_prefix=""){
    df_prot_annot = read.csv(prot_annotations_path, sep="\t")
    df_pep_annot = read.csv(pep_annotations_path, sep="\t")
    df_prot_annot$Origin = "prot"
    df_pep_annot$Origin = "pep"

    df_annots = data.frame(Protein.IDs=c(df_prot_annot$Protein.IDs, df_pep_annot$Leading.razor.protein))

    l_uniprot_ids = list()
    for (prot_id in df_prot_annot$Protein.IDs){
        main_prot = unlist(strsplit(prot_id,";"))[1]
        l_uniprot_ids = append(l_uniprot_ids, main_prot)
    }

    df_annots$Uniprot.Id = c(unlist(l_uniprot_ids), df_pep_annot$Leading.razor.protein)
    df_annots$Protein.names = c(df_prot_annot$Protein.names, df_pep_annot$Protein.names)
    df_annots$Gene.names = c(df_prot_annot$Gene.names, df_pep_annot$Gene.names)
    df_annots$Origin = c(df_prot_annot$Origin, df_pep_annot$Origin)
    rownames(df_annots) = c(rownames(df_prot_annot), rownames(df_pep_annot))

    output_annot_path = file.path(output_dir, "annotations.txt")
    message(paste0("Writing ", output_annot_path, " file"))
    write.table(df_annots, file=output_annot_path, sep="\t", append=F, quote=F)

    df_prot_data = read.csv(prot_data_path, sep="\t")
    df_pep_data = read.csv(pep_data_path, sep="\t")
    df_data = data.frame(todel=c(df_prot_data[,colnames(df_prot_data)[1]], df_pep_data[,colnames(df_prot_data)[1]]))

    for (col in colnames(df_prot_data)){
        df_data[[col]] = c(df_prot_data[,col], df_pep_data[,col])
    }

    df_data$todel <- NULL
    rownames(df_data) = c(rownames(df_prot_data), rownames(df_pep_data))
    output_data_path = file.path(output_dir, "data_before_normalization.txt")
    message(paste0("Writing ", output_data_path, " file"))
    write.table(df_data, file=output_data_path, sep="\t", append=F, quote=F)

    norm_output_path <- file.path(output_dir, "data_after_normalization.txt")
    df_normalized <- prot_normalize(output_data_path, norm_output_path, conditions_new_sample_names_path, norm_method)

    impute_output_path <- file.path(output_dir, output_both_data)
    df_imputed <- prot_impute(norm_output_path, impute_output_path, imputation_method, k=k, rowmax=rowmax, colmax=colmax)
}