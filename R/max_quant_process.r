# Build first proteome file:
# - new sample names
# - keep only samples from the annotation file
#
# Build annotation file with new names

prot_max_quant_process <- function ( proteinGroups_path, conditions_path, output_prot_file, samples__old_names__new_names_file, prot_annotations_path, remove_reverse_identified_contaminant=TRUE ){
    df_prot <- read.csv(proteinGroups_path, sep="\t")
    df_annotations <- read.csv(conditions_path, sep="\t")

    if ( remove_reverse_identified_contaminant ){
        df_prot<-df_prot[!(df_prot$Reverse=="+" | df_prot$Only.identified.by.site=="+" | df_prot$Potential.contaminant=="+"),]
    }

    # Building new names, "NCY38" becomes "CR_1", "NCY54" becomes "CR_2" ...
    c_cond_occ = c()
    for ( cond in df_annotations[, 'condition'] ){
        c_cond_occ[[cond]] = 1
    }
    c_col_new_names = c()
    count = 1

    if ( length(colnames(df_annotations)) == 3 ){
        for ( label in df_annotations[,'label'] ){
            c_col_new_names[[label]] = paste0( df_annotations[count, 'condition'], "_", c_cond_occ[[df_annotations[count, 'condition']]] )
            c_cond_occ[[df_annotations[count, 'condition']]] = c_cond_occ[[df_annotations[count, 'condition']]] + 1
            count = count + 1
        }
    } else {
        count = 1
        for ( label in df_annotations[,'label'] ){
            c_col_new_names[[label]] = df_annotations[count,'new_name']
            count = count + 1
        }
    }

    # Keep interesting columns, only LFQ.intensity.XXX and only columns in the annotations file
    l_kept_cols = list()
    l_ind_kept_cols = list()
    l_finals_cols = colnames(df_prot)

    for ( ind in grep("LFQ.", colnames(df_prot)) ){
        name = colnames(df_prot)[ind]
        if ( unlist(strsplit(name,"\\."))[3] %in% df_annotations[, 'label'] ){
            l_kept_cols = append(l_kept_cols, name)
            l_ind_kept_cols = append(l_ind_kept_cols, ind)
        }
    }

    df_just_samples <- df_prot[,unlist(l_kept_cols)]
    l_new_names = list()
    for ( old_name in str_replace( colnames(df_just_samples), "LFQ.intensity.", "" ) ){
        l_new_names = append( l_new_names, c_col_new_names[old_name] )
        #print( paste( old_name, c_col_new_names[old_name] ) )
    }

    colnames(df_just_samples) <- unlist(l_new_names)

    # Create unique ID
    l_uniqueID = list()
    count = 1
    for ( prot_id in df_prot[,'Protein.IDs']){
        main_prot = unlist(strsplit(prot_id,";"))[1]
        l_uniqueID = append(l_uniqueID,paste0(main_prot, "//prot//", rownames(df_prot)[count]))
        count = count + 1
    }

    rownames(df_just_samples) = unlist(l_uniqueID)
    df_just_samples[df_just_samples == 0] <- NA

    message(paste0("Writing ", output_prot_file, " file"))
    write.table(df_just_samples, file=output_prot_file, sep="\t", append=F, quote=F)

    l_ordered_new_samples_name = list()

    for (label in df_annotations$label){
        l_ordered_new_samples_name = append(l_ordered_new_samples_name, c_col_new_names[label])
    }

    df_annotations$label = unlist(l_ordered_new_samples_name)
    if ( length(colnames(df_annotations)) == 4){
        df_annotations$new_name = NULL
    }

    message(paste0("Writing ", samples__old_names__new_names_file, " file"))
    write.table(df_annotations, file=samples__old_names__new_names_file, sep="\t", append=F, quote=F, row.names=F)

    df_prot_annotations = df_prot[, c("Protein.IDs", "Protein.names", "Gene.names")]
    rownames(df_prot_annotations) = unlist(l_uniqueID)
    message(paste0("Writing ", prot_annotations_path, " file"))
    write.table(df_prot_annotations, file=prot_annotations_path, sep="\t", append=F, quote=F)

    return(df_just_samples)
}

pep_max_quant_process <- function ( proteinGroups_path, conditions_path, output_prot_file, samples__old_names__new_names_file, prot_annotations_path, remove_reverse_identified_contaminant=TRUE){
    df_prot <- read.csv(proteinGroups_path, sep="\t")
    df_annotations <- read.csv(conditions_path, sep="\t")

    if ( remove_reverse_identified_contaminant ){
        df_prot<-df_prot[!(df_prot$Reverse == "+" | df_prot$Potential.contaminant == "+"),]
    }

    # Building new names, "NCY38" becomes "CR_1", "NCY54" becomes "CR_2" ...
    c_cond_occ = c()
    for ( cond in df_annotations[, 'condition'] ){
        c_cond_occ[[cond]] = 1
    }
    c_col_new_names = c()
    count = 1
    

    if ( length(colnames(df_annotations)) == 3 ){
        for ( label in df_annotations[,'label'] ){
            c_col_new_names[[label]] = paste0( df_annotations[count, 'condition'], "_", c_cond_occ[[df_annotations[count, 'condition']]] )
            c_cond_occ[[df_annotations[count, 'condition']]] = c_cond_occ[[df_annotations[count, 'condition']]] + 1
            count = count + 1
        }
    } else {
        count = 1
        for ( label in df_annotations[,'label'] ){
            c_col_new_names[[label]] = df_annotations[count,'new_name']
            count = count + 1
        }
    }

    # Keep interesting columns, only LFQ.intensity.XXX and only columns in the annotations file
    l_kept_cols = list()
    l_ind_kept_cols = list()
    l_finals_cols = colnames(df_prot)

    for ( ind in grep("LFQ.", colnames(df_prot)) ){
        name = colnames(df_prot)[ind]
        if ( unlist(strsplit(name,"\\."))[3] %in% df_annotations[, 'label'] ){
            l_kept_cols = append(l_kept_cols, name)
            l_ind_kept_cols = append(l_ind_kept_cols, ind)
        }
    }

    df_just_samples <- df_prot[,unlist(l_kept_cols)]
    l_new_names = list()
    for ( old_name in str_replace( colnames(df_just_samples), "LFQ.intensity.", "" ) ){
        l_new_names = append( l_new_names, c_col_new_names[old_name] )
        #print( paste( old_name, c_col_new_names[old_name] ) )
    }

    colnames(df_just_samples) <- unlist(l_new_names)

    # Create unique ID
    l_uniqueID = list()
    count = 1
    for ( prot_id in df_prot[,"Leading.razor.protein"] ){
        main_prot = unlist(strsplit(prot_id,";"))[1]
        l_uniqueID = append(l_uniqueID,paste0(main_prot, "//pep//", rownames(df_prot)[count]))
        count = count + 1
    }

    rownames(df_just_samples) = unlist(l_uniqueID)
    df_just_samples[df_just_samples == 0] <- NA

    message(paste0("Writing ", output_prot_file, " file"))
    write.table(df_just_samples, file=output_prot_file, sep = "\t", append=F, quote=F)

    l_ordered_new_samples_name = list()

    for (label in df_annotations$label){
        l_ordered_new_samples_name = append(l_ordered_new_samples_name, c_col_new_names[label])
    }

    df_annotations$label = unlist(l_ordered_new_samples_name)
    if ( length(colnames(df_annotations)) == 4){
        df_annotations$new_name = NULL
    }

    message(paste0("Writing ", samples__old_names__new_names_file, " file"))
    write.table(df_annotations, file=samples__old_names__new_names_file, sep = "\t", append=F, quote=F, row.names=F)

    df_prot_annotations = df_prot[, c("Leading.razor.protein", "Protein.names", "Gene.names")]
    rownames(df_prot_annotations) = unlist(l_uniqueID)
    message(paste0("Writing ", prot_annotations_path, " file"))
    write.table(df_prot_annotations, file=prot_annotations_path, sep = "\t", append=F, quote=F)

    return(df_just_samples)
}