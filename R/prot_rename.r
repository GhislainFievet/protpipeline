prot_rename <- function(impute_output_path, partial_impute_output_path, rename_output_path,
                    partial_rename_output_path, prot_rename_db_path){
    print("prot_rename()")
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
    print("Writing rename_output_path")
    print(rename_output_path)
    write.table(df_prot, rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)

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
    print("Writing partial_rename_output_path")
    print(partial_rename_output_path)
    write.table(df_prot, partial_rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)
}

prot_rename_biomart <- function(impute_output_path, partial_impute_output_path, rename_output_path,
                    partial_rename_output_path){

    df_prot <- read.table(impute_output_path)

    uniprot_ids <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
    }))
    uniprot_ids <- sub("-\\d+$", "", uniprot_ids)

    # Use biomaRt to get the HGNC symbols
    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    results <- biomaRt::getBM(
        attributes = c("uniprotswissprot", "hgnc_symbol"),
        filters = "uniprotswissprot",
        values = uniprot_ids,
        mart = ensembl
    )

    print(head(results))
    print(dim(results))
    new_index <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
        my_base <- unlist(strsplit(my_base, "-"))[1]
        my_id <- paste0(unlist(strsplit(x, "//"))[2], "//", unlist(strsplit(x, "//"))[3])
        if (!my_base %in% results$uniprotswissprot || results$hgnc_symbol[results$uniprotswissprot == my_base][1] == ""){
            return(x)
        } else {
            return(paste0(results$hgnc_symbol[results$uniprotswissprot == my_base][1], "//", my_id))
        }
    }))

    rownames(df_prot) <- new_index
    write.table(df_prot, rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)

    df_prot <- read.table(partial_impute_output_path)

     uniprot_ids <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
    }))
    uniprot_ids <- sub("-\\d+$", "", uniprot_ids)

    # Use biomaRt to get the HGNC symbols
    results <- biomaRt::getBM(
        attributes = c("uniprotswissprot", "hgnc_symbol"),
        filters = "uniprotswissprot",
        values = uniprot_ids,
        mart = ensembl
    )

    print(head(results))
    print(dim(results))
    new_index <- unlist(lapply(rownames(df_prot), function(x){
        my_base <- unlist(strsplit(x, "//"))[1]
        my_base <- unlist(strsplit(my_base, "-"))[1]
        my_id <- paste0(unlist(strsplit(x, "//"))[2], "//", unlist(strsplit(x, "//"))[3])
        if (!my_base %in% results$uniprotswissprot || results$hgnc_symbol[results$uniprotswissprot == my_base][1] == ""){
            return(x)
        } else {
            return(paste0(results$hgnc_symbol[results$uniprotswissprot == my_base][1], "//", my_id))
        }
    }))

    rownames(df_prot) <- new_index
    print(paste("Writing", partial_rename_output_path))
    write.table(df_prot, partial_rename_output_path, sep = "\t", quote = FALSE, row.names = TRUE)
}