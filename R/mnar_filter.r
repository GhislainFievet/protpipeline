
prot_mnar_filter <- function (prot_path, conditions_path, my_prefix, MNAR_threshold){
    df_prot = read.csv(prot_path, sep="\t")
    
    df_conditions = read.csv(conditions_path, sep="\t")

    df_thresholds_inf = data.frame(matrix(ncol=length(unique(df_conditions$condition)), nrow=1))
    colnames(df_thresholds_inf) = unique(df_conditions$condition)
    for ( i in 1:unique(df_conditions$condition) ){
        df_thresholds_inf[1,i] = sum(df_conditions$condition == unique(df_conditions$condition)[i])*MNAR_threshold
    }
    df_thresholds_sup = data.frame(matrix(ncol=length(unique(df_conditions$condition)), nrow=1))
    colnames(df_thresholds_sup) = unique(df_conditions$condition)
    for ( i in 1:unique(df_conditions$condition) ){
        df_thresholds_sup[1,i] = sum(df_conditions$condition == unique(df_conditions$condition)[i])*(1-MNAR_threshold)
    }

    # Create samples.in.conditions, make columns list by condition
    samples.in.conditions = c()

    for ( cond in unique(df_conditions$condition) ){
        samples.in.conditions[[cond]] = df_conditions[df_conditions$condition==cond, "label"]
    }

    condition.filter_inf <- sapply(colnames(df_thresholds), function(x) {
        apply(as.matrix(df_prot[samples.in.conditions[[x]]]), 1, function(y) length(which(!is.na(y))) <= df_thresholds_inf[1,x])
    })
    condition.filter_sup <- sapply(colnames(df_thresholds), function(x) {
        apply(as.matrix(df_prot[samples.in.conditions[[x]]]), 1, function(y) length(which(!is.na(y))) >= df_thresholds_sup[1,x])
    })

    keepProt_inf <- apply( condition.filter_inf, 1, any )
    keepProt_sup <- apply( condition.filter_sup, 1, any )
    keepProt <- keepProt_inf & keepProt_sup

    prot_final <- df_prot[keepProt,]

    message(paste0("Writing ", output_path, " file"))
    write.table(prot_final, file=output_path, sep="\t", append=F, quote=F)

    return(prot_final)
}

pep_mnar_filter <- function (prot_path, conditions_path,
                thresholds_path, pep_annotations_path,
                output_path, peptide_occurence_filter=1,
                atLeastOne=TRUE, my_prefix="", group_threshold_mode,
                group_threshold, MNAR_threshold, MNAR_filter){
    message(paste(prot_path, conditions_path, thresholds_path, pep_annotations_path))
    df_prot = read.csv(prot_path, sep="\t")
    df_conditions = read.csv(conditions_path, sep="\t")
    if ( group_threshold_mode == "file"){
        df_thresholds = read.csv(thresholds_path, sep="\t")
    } else {
        df_thresholds = data.frame(matrix(ncol=length(unique(df_conditions$condition)), nrow=1))
        colnames(df_thresholds) = unique(df_conditions$condition)
        for ( i in 1:unique(df_conditions$condition) ){
            df_thresholds[1,i] = sum(df_conditions$condition == unique(df_conditions$condition)[i])*group_threshold
        }
    }
    df_annot = read.csv(pep_annotations_path, sep="\t")

    # Create samples.in.conditions, make columns list by condition
    samples.in.conditions = c()

    for ( cond in unique(df_conditions$condition) ){
        samples.in.conditions[[cond]] = df_conditions[df_conditions$condition==cond, "label"]
    }

    condition.filter <- sapply(colnames(df_thresholds), function(x) {
        # print(x)
        # print(colnames(df_prot))
        # print(samples.in.conditions[[x]])
        apply(as.matrix(df_prot[samples.in.conditions[[x]]]), 1, function(y) length(which(!is.na(y))) >= df_thresholds[1,x])
    })
    
    keepProt <- if(atLeastOne) {
        apply( condition.filter, 1, any )
    } else {
        apply( condition.filter, 1, all )
    }

    prot_final <- df_prot[keepProt,]

    prot_with_id = prot_final
    prot_with_id$id = rownames(prot_with_id)
    df_annot$id = rownames(df_annot)

    df_full_annots = merge(x = prot_with_id, y = df_annot, by = "id", all.x = TRUE)

    if ( peptide_occurence_filter != 1 ){
        c_prot_occ = c()
        for ( id in df_full_annots$Leading.razor.protein ){
            c_prot_occ[[id]] = 0
        }
        for ( id in df_full_annots$Leading.razor.protein ){
            c_prot_occ[[id]] = c_prot_occ[[id]] + 1
        }
        l_kept_pep = list()
        count = 1
        for (id in df_full_annots$Leading.razor.protein){
            if (c_prot_occ[[id]]>=peptide_occurence_filter){
                l_kept_pep = append(l_kept_pep, df_full_annots$id[count])
            }
            count = count + 1
        }
        prot_final <- prot_final[unlist(l_kept_pep),]
    }

    message(paste0("Writing ", output_path, " file"))
    write.table(prot_final, file=output_path, sep="\t", append=F, quote=F)

    return(prot_final)
}

pep_mnar_filter_both <- function (pep_path, prot_path, conditions_path,
        thresholds_path, pep_annotations_path, prot_annotations_path,
        output_path, peptide_occurence_filter=1, atLeastOne=TRUE, my_prefix="",
        group_threshold_mode, group_threshold, MNAR_threshold, MNAR_filter){
    message(paste(pep_path, prot_path, conditions_path, thresholds_path, pep_annotations_path, prot_annotations_path))
    df_pep = read.csv(pep_path, sep="\t")
    df_prot = read.csv(prot_path, sep="\t")
    df_conditions = read.csv(conditions_path, sep="\t")
    if ( group_threshold_mode == "file"){
        df_thresholds = read.csv(thresholds_path, sep="\t")
    } else {
        df_thresholds = data.frame(matrix(ncol=length(unique(df_conditions$condition)), nrow=1))
        colnames(df_thresholds) = unique(df_conditions$condition)
        for ( i in 1:unique(df_conditions$condition) ){
            df_thresholds[1,i] = sum(df_conditions$condition == unique(df_conditions$condition)[i])*group_threshold
        }
    }
    df_annot = read.csv(pep_annotations_path, sep="\t")
    df_prot_annot = read.csv(prot_annotations_path, sep="\t")
    # Create samples.in.conditions, make columns list by condition
    samples.in.conditions = c()
    for ( cond in unique(df_conditions$condition) ){
        samples.in.conditions[[cond]] = df_conditions[df_conditions$condition==cond, "label"]
    }
    condition.filter <- sapply(colnames(df_thresholds), function(x) {
        # print(x)
        # print(colnames(df_pep))
        # print(samples.in.conditions[[x]])
        apply(as.matrix(df_pep[samples.in.conditions[[x]]]), 1, function(y) length(which(!is.na(y))) >= df_thresholds[1,x])
    })
    
    keepProt <- if(atLeastOne) {
        apply( condition.filter, 1, any )
    } else {
        apply( condition.filter, 1, all )
    }
    pep_final <- df_pep[keepProt,]
    pep_with_id = pep_final
    pep_with_id$id = rownames(pep_with_id)
    df_annot$id = rownames(df_annot)
    df_full_annots = merge(x = pep_with_id, y = df_annot, by = "id", all.x = TRUE)
    prot_with_id = df_prot
    prot_with_id$id = rownames(prot_with_id)
    df_prot_annot$id = rownames(df_prot_annot)
    df_prot_full_annot = merge(x = prot_with_id, y = df_prot_annot, by = "id", all.x = TRUE)

    if ( peptide_occurence_filter != 1 ){
        l_prot_ids = list()
        for ( prot_ids in df_prot_full_annot$Protein.IDs){
            for (prot_id in unlist(strsplit(prot_ids,";"))){
                l_prot_ids = append(l_prot_ids, prot_id)
            }
        }
        c_prot_occ = c()
        
        for ( id in df_full_annots$Leading.razor.protein ){
            c_prot_occ[[id]] = 0
        }
        for ( id in df_full_annots$Leading.razor.protein ){
            c_prot_occ[[id]] = c_prot_occ[[id]] + 1
        }

        l_kept_pep = list()
        count = 1
        for (id in df_full_annots$Leading.razor.protein){
            if (c_prot_occ[[id]]>=peptide_occurence_filter){
                if (!(id %in% l_prot_ids)){
                    l_kept_pep = append(l_kept_pep, df_full_annots$id[count])
                }
            }
            count = count + 1
        }
        pep_final <- pep_final[unlist(l_kept_pep),]
    }

    message(paste0("Writing ", output_path, " file"))
    write.table(pep_final, file=output_path, sep="\t", append=F, quote=F)

    return(pep_final)
}