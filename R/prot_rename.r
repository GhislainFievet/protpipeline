prot_rename <- function(impute_output_path, partial_impute_output_path, rename_output_path,
                    partial_rename_output_path, prot_rename_db_path){
    df_as <- read.csv(prot_rename_db_path, sep="\t")

    l_protID_hgnc <<- list()
    my_apply <- apply(df_as[nrow(df_as):1,], 1, function(x){
        a_1 <- unlist(strsplit(x['Protein.IDs'], ";"))
        a_2 <- unlist(strsplit(x['HGNC.symbol_result'], ";"))
        for (str1 in a_1){
            for (str2 in a_2[length(a_2):1]){
                l_protID_hgnc[[str1]] <<- str2
            }
        }
    })

    df_prot <- read.table(impute_output_path)
    new_index <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
        my_id <- paste0(unlist(strsplit(x, "//"))[2], "//", unlist(strsplit(x, "//"))[3])
        if (!my_base %in% names(l_protID_hgnc)){
            return(x)
        } else {
            return(paste0(l_protID_hgnc[[my_base]], "//", my_id))
        }
    }))
    rownames(df_prot) <- new_index
    write.table(df_prot, file = rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)

    df_prot <- read.table(partial_impute_output_path)
    new_index <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
        my_id <- paste0(unlist(strsplit(x, "//"))[2], "//", unlist(strsplit(x, "//"))[3])
        if (!my_base %in% names(l_protID_hgnc)){
            return(x)
        } else {
            return(paste0(l_protID_hgnc[[my_base]], "//", my_id))
        }
    }))
    rownames(df_prot) <- new_index
    write.table(df_prot, file = partial_rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)
}