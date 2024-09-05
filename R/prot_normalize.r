
prot_normalize <- function(prot_path, output_path, conditions_path, norm_method){
    message("")
    message( paste("Normalization method:", norm_method) )
    message("")

    df_temp = read.csv(prot_path, sep="\t")
    df_conditions_new_sample_names = read.csv(conditions_path, sep="\t")

    numb = length(colnames(df_temp))
    df_temp$line.id = 1:length(rownames(df_temp))
    data_unique <-  make_unique( df_temp, "line.id", "line.id", delim = ";" )
    se_after_filter <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    mat_norm <- se_after_filter@assays@data[[1]]

    if ( norm_method == "vsn" || norm_method == "vsn+quantiles" ){
        norm_vsn <- normalize_vsn(se_after_filter)
        mat_norm <- norm_vsn@assays@data[[1]]
    }

    if ( norm_method == "quantiles" || norm_method == "vsn+quantiles" ){
        dataNormVSn <- as.data.frame(mat_norm)
        mat_norm <- normalize.quantiles(as.matrix(dataNormVSn), copy=FALSE)
    }

    rownames(mat_norm) <- rownames(df_temp)

    c_sample_name_correspondance = c()
    count = 1
    for (cond in df_conditions_new_sample_names$condition){
        name = paste0(cond,"_",df_conditions_new_sample_names$replicate[count])
        c_sample_name_correspondance[[name]] = df_conditions_new_sample_names$label[count]
        count = count + 1
    }
    l_new_sample_colnames = list()
    for ( col in colnames(mat_norm) ){
        l_new_sample_colnames = append(l_new_sample_colnames, c_sample_name_correspondance[[col]])
    }
    colnames(mat_norm) <- unlist(l_new_sample_colnames)
    #colnames(mat_norm) <- colnames(df_temp)[1:length(colnames(mat_norm))]

    if (norm_method == "none"){
        message(paste0("Writing ", output_path, " file"))
        df_results = log2(df_temp[,colnames(df_temp)[1:length(colnames(mat_norm))]])
        write.table( df_results, output_path, sep="\t", append=F, quote=F )

        return(df_results)
    } else { 
    message(paste0("Writing ", output_path, " file"))
    write.table( mat_norm, output_path, sep="\t", append=F, quote=F )
    }



    return(mat_norm)
}
