prot_plots <- function(before_filter_path,
                        after_filter_path,
                        after_normalization_path,
                        after_imputation_path,
                        conditions_path,
                        plot_dir,
                        prot_or_pep,
                        my_prefix=""){
    prot_or_pep_2 = "Proteins"
    prot_or_pep_3 = "proteins"
    prot_or_pep_4 = "proteins"
    prot_or_pep_5 = "Protein"
    if ( prot_or_pep == "pep"){
        prot_or_pep_2 = "Peptides"
        prot_or_pep_3 = "peptides"
        prot_or_pep_4 = "peptides"
        prot_or_pep_5 = "Peptide"
    }
    if ( prot_or_pep == "both"){
        prot_or_pep_2 = "Proteins/Peptides"
        prot_or_pep_3 = "proteins/peptides"
        prot_or_pep_4 = "proteins_peptides"
        prot_or_pep_5 = "Protein or peptide"
    }

    output_dir <- plot_dir

    ###### load files ######
    df_before_filter = read.csv(before_filter_path, sep="\t")
    df_after_filter = read.csv(after_filter_path, sep="\t")
    df_after_normalization = read.csv(after_normalization_path, sep="\t")
    df_after_imputation = 0
    if (file.exists(after_imputation_path)){
        df_after_imputation = read.csv(after_imputation_path, sep="\t")
    }
    df_conditions_new_sample_names = read.csv(conditions_path, sep="\t")

    numb = length(colnames(df_before_filter))
    df_before_filter$line.id = 1:length(rownames(df_before_filter))
    data_unique <-  make_unique( df_before_filter, "line.id", "line.id", delim = ";" )
    se_before_filter <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    numb = length(colnames(df_after_filter))
    df_after_filter$line.id = 1:length(rownames(df_after_filter))
    data_unique <-  make_unique( df_after_filter, "line.id", "line.id", delim = ";" )
    se_after_filter <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    numb = length(colnames(df_after_normalization))
    df_after_normalization$line.id = 1:length(rownames(df_after_normalization))
    data_unique <-  make_unique( df_after_normalization, "line.id", "line.id", delim = ";" )
    se_after_normalization <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    if (file.exists(after_imputation_path)){
        numb = length(colnames(df_after_imputation))
        df_after_imputation$line.id = 1:length(rownames(df_after_imputation))
        data_unique <-  make_unique( df_after_imputation, "line.id", "line.id", delim = ";" )
        se_after_imputation <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )
    }


    ###### Export plots ######
    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_nb_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_numbers(se_after_filter) + labs(title = paste0(prot_or_pep_2, " per sample"), x = "",
         y = paste0("Number of ", prot_or_pep_3)))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_nb_before_filter.png"))
    message("Writing ", plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_numbers(se_before_filter) + labs(title = paste0(prot_or_pep_2, " per sample"), x = "",
         y = paste0("Number of ", prot_or_pep_3)))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_miss_map_after_filter.png"))
    message("Writing ", plot_path, " plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_missval(se_after_filter))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_miss_map_before_filter.png"))
    if ( length(rownames(df_before_filter)) < 10000){
        message("Writing ", plot_path, " plot")
        png(filename=plot_path, width = 1000, height = 800)
        print(plot_missval(se_before_filter))
        dev.off()
    } else {
        message("Too much data to plot ", plot_path)
    }

    temp_before_filter = se_before_filter
    temp_after_filter = se_after_filter
    temp_before_filter@assays@data[[1]] = log2(temp_before_filter@assays@data[[1]])
    temp_after_filter@assays@data[[1]] = log2(temp_after_filter@assays@data[[1]])
    if (file.exists(after_imputation_path)){
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_normalize.png"))
        message("Writing ",plot_path," plot")
        png(filename=plot_path, width = 1000, height = 800)

        print(plot_normalization(temp_before_filter, temp_after_filter, se_after_normalization, se_after_imputation))
        dev.off()
    } else {
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_normalize.png"))
        message("Writing ",plot_path," plot")
        png(filename=plot_path, width = 1000, height = 800)
        print(plot_normalization(temp_before_filter, temp_after_filter, se_after_normalization))
        dev.off()
    }

    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_normalize_before_vs_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_normalization(se_before_filter, se_after_filter))
    dev.off()

    if (file.exists(after_imputation_path)){
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_normalize_normalized_vs_imputation.png"))
        message("Writing ",plot_path," plot")
        png(filename=plot_path, width = 1000, height = 800)
        print(plot_normalization(se_after_normalization, se_after_imputation))
        dev.off()
    } else {
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_normalize_normalized_vs_imputation.png"))
        message("Writing ",plot_path," plot")
        png(filename=plot_path, width = 1000, height = 800)
        print(plot_normalization(se_after_normalization))
        dev.off()
    }

    if (file.exists(after_imputation_path)){
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_conditions_distributions.png"))
        message("Writing ", plot_path, " plot")
        png(filename= plot_path, width = 1000, height = 800)
        temp_before_filter = se_before_filter
        temp_after_filter = se_after_filter
        temp_before_filter@assays@data[[1]] = log2(temp_before_filter@assays@data[[1]])
        temp_after_filter@assays@data[[1]] = log2(temp_after_filter@assays@data[[1]])
        print(plot_imputation(temp_before_filter, temp_after_filter, se_after_normalization, se_after_imputation))
        dev.off()
    } else {
        plot_path <- file.path(output_dir, paste0(prot_or_pep, "_conditions_distributions.png"))
        message("Writing ", plot_path, " plot")
        png(filename= plot_path, width = 1000, height = 800)
        temp_before_filter = se_before_filter
        temp_after_filter = se_after_filter
        temp_before_filter@assays@data[[1]] = log2(temp_before_filter@assays@data[[1]])
        temp_after_filter@assays@data[[1]] = log2(temp_after_filter@assays@data[[1]])
        print(plot_imputation(temp_before_filter, temp_after_filter, se_after_normalization))
        dev.off()
    }

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions_before_vs_after_filter.png"))
    message("Writing ", plot_path, " plot")
    png(filename= plot_path, width = 1000, height = 800)
    print(plot_imputation(se_before_filter, se_after_filter))
    dev.off()

    if (file.exists(after_imputation_path)){
        plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions_normalized_vs_imputation.png"))
        message("Writing ", plot_path, " plot")
        png(filename= plot_path, width = 1000, height = 800)
        print(plot_imputation(se_after_normalization, se_after_imputation))
        dev.off()
    } else {
        plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions_normalized_vs_imputation.png"))
        message("Writing ", plot_path, " plot")
        png(filename= plot_path, width = 1000, height = 800)
        print(plot_imputation(se_after_normalization))
        dev.off()
    }
    
    ## Missing values plot
    tab_after <- as.data.frame(se_after_filter@assays@data[[1]])
    tab_before <- as.data.frame(se_before_filter@assays@data[[1]])

    # Proportion of missing data by sample
    missing.values <- tab_after %>%
            gather(key = "key", value = "val") %>%
            mutate(isna = is.na(val)) %>%
            group_by(key) %>%
            mutate(total = n()) %>%
            group_by(key, total, isna) %>%
            summarise(num.isna = n()) %>%
            mutate(pct = num.isna / total * 100)

    levels <- (missing.values  %>% filter(isna == T) %>% arrange(desc(pct)))$key

    percentage.plot <- missing.values %>%
        ggplot() + geom_bar(aes(x = key, y = pct, fill=isna), stat = 'identity', alpha=0.8) +
        scale_x_discrete(limits = levels) +
        scale_fill_manual(name = "",values = c('darkslateblue', 'gray53'), labels = c("Present", "Missing")) +
        coord_flip() +
        labs(title = "Ratio of missing value by sample", x =' Samples', y = "% of missing value")

    plot_path <- file.path(output_dir,paste0(prot_or_pep,"_missing_values_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(percentage.plot)
    dev.off()

    row.plot <- tab_after %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
    values = c('darkslateblue', 'gray53'),
    labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Samples", y = prot_or_pep_2, title = paste0("Missing ",prot_or_pep_3," by sample, after filter")) +
    coord_flip()

    plot_path <- file.path(output_dir,paste0(prot_or_pep,"_missing_",prot_or_pep_4,"_by_sample_map_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(row.plot)
    dev.off()

    row.plot <- tab_before %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
    values = c('darkslateblue', 'gray53'),
    labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Sample", y = prot_or_pep_2, title = paste0("Missing ", prot_or_pep_3, " by sample, before filter")) + coord_flip()

    
    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_missing_", prot_or_pep_4, "_by_sample_map_before_filter.png"))
    if (prot_or_pep == "prot"){
        message("Writing ", plot_path, " plot")
        png(filename=plot_path, width = 1000, height = 800)
        print(row.plot)
        dev.off()
    } else {
        message("Not creating ", plot_path, " because processing peptides")
    }

    if (prot_or_pep == "prot"){
        # computations for prot_correlation_plot_filter.png and prot_distribution_plot_filter.png
        if(length(tab_after[,1:length(colnames(tab_after))])!=1){
            mediane <- apply(tab_after[,1:length(colnames(tab_after))],1,median,na.rm=TRUE)
            nb.categories <- trunc(log10(dim(tab_after[,1:length(colnames(tab_after))][,1:length(colnames(tab_after))])[1]))
            stats.matrix <- matrix(nrow=1+2*(nb.categories+1),ncol=dim(tab_after[,1:length(colnames(tab_after))])[2])
            colnames(stats.matrix) <- colnames(tab_after[,1:length(colnames(tab_after))])
            intervals <- 0:nb.categories
            petites.valeurs <- paste("pv-",as.character(10^intervals),sep="")
            grandes.valeurs <- paste("gv-",as.character(10^rev(intervals)),sep="")
            rownames(stats.matrix) <- c("coeff.corr",petites.valeurs,grandes.valeurs)
            for(i in (1:dim(tab_after[,1:length(colnames(tab_after))])[2])){
                stats.matrix[1,i] <- cor(tab_after[,1:length(colnames(tab_after))][,1:length(colnames(tab_after))][,i],mediane,use="complete.obs")
                mat <- sort(tab_after[,1:length(colnames(tab_after))][,i],na.last=NA)
                for(j in intervals){
                    stats.matrix[j+2,i] <- mat[10^j]
                    stats.matrix[dim(stats.matrix)[1]-j,i] <- mat[length(mat)-10^j]
                }
            }
        }

        p_75 <- tab_after
        plot_path <- file.path(output_dir,paste0(prot_or_pep,"_correlation_plot_filter.png"))
        message("Writing ",plot_path," plot")
        png(file=plot_path, width = 1000, height = 800)
        correlation_plot_filtre <- plot( stats.matrix[1,], ylab="correlation coeff.", ylim = c(0,1), xlab="samples", lwd=1, col="blue", type="l", xaxt="n", main="Correlation to median sample expression" )
        axis(side=1, labels = colnames( p_75[,1:length(colnames(tab_after))]), at = 1:length(colnames(tab_after)), cex.axis =  0.55, las = 2 )
        print(correlation_plot_filtre)
        dev.off()

        plot_path <- file.path(output_dir,paste0(prot_or_pep,"_distribution_plot_filter.png"))
        message("Writing ",plot_path," plot")
        png(file=plot_path, width = 1000, height = 800)
        min_born <- floor(min(stats.matrix[2,])) -5
        max_born <- ceiling(max(stats.matrix[9,]))+20
        distrib_plot_filtre <- plot(stats.matrix[2,], ylab="Expression intensity", ylim = c(min_born, max_born), xlab="Samples", lwd=1, col="blue", type="l", axes=FALSE, main=paste0(prot_or_pep_5, " expression quantiles"))
        axis(2)
        lines(stats.matrix[3,], col="Red" )
        lines(stats.matrix[4,], col="Darkgreen" )
        lines(stats.matrix[5,], col="Green" )
        lines(stats.matrix[6,], col="Orange" )
        lines(stats.matrix[7,], col="Grey" )
        lines(stats.matrix[8,], col="Black" )
        lines(stats.matrix[9,], col="Purple" )
        legend(1, max_born, legend=c("PV", "PV10", "PV100", "PV1000", "PG1000", "PG100", "PG10", "PG"), col=c("Blue", "Red", "Darkgreen", "Green", "Orange", "Grey","Black","Purple"), lty = c(1))
        axis(side=1, labels = colnames(p_75[,1:length(colnames(tab_after))]), at = 1:length(colnames(tab_after)), cex.axis = 0.75, las=2)
        print(distrib_plot_filtre)
        dev.off()
    } else {
        message("Not creating ", plot_path, " because processing peptides")
    }
}

concat_plots <- function(after_filter_path, after_normalization_path, after_imputation_path,
        conditions_path, plot_dir, prot_or_pep, my_prefix, output_dir){
    prot_or_pep_2 = "Proteins"
    prot_or_pep_3 = "proteins"
    prot_or_pep_4 = "proteins"
    prot_or_pep_5 = "Protein"
    if ( prot_or_pep == "pep"){
        prot_or_pep_2 = "Peptides"
        prot_or_pep_3 = "peptides"
        prot_or_pep_4 = "peptides"
        prot_or_pep_5 = "Peptide"
    }
    if ( prot_or_pep == "both"){
        prot_or_pep_2 = "Proteins/Peptides"
        prot_or_pep_3 = "proteins/peptides"
        prot_or_pep_4 = "proteins_peptides"
        prot_or_pep_5 = "Protein or peptide"
    }

    ###### load files ######
    df_after_filter = read.csv(after_filter_path, sep="\t")
    df_after_normalization = read.csv(after_normalization_path, sep="\t")
    df_after_imputation = read.csv(after_imputation_path, sep="\t")
    df_conditions_new_sample_names = read.csv(conditions_path, sep="\t")

    numb = length(colnames(df_after_filter))
    df_after_filter$line.id = 1:length(rownames(df_after_filter))
    data_unique <-  make_unique( df_after_filter, "line.id", "line.id", delim = ";" )
    se_after_filter <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    numb = length(colnames(df_after_normalization))
    df_after_normalization$line.id = 1:length(rownames(df_after_normalization))
    data_unique <-  make_unique( df_after_normalization, "line.id", "line.id", delim = ";" )
    se_after_normalization <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )

    numb = length(colnames(df_after_imputation))
    df_after_imputation$line.id = 1:length(rownames(df_after_imputation))
    data_unique <-  make_unique( df_after_imputation, "line.id", "line.id", delim = ";" )
    se_after_imputation <- make_se( data_unique, unlist(1:numb), df_conditions_new_sample_names )


    ###### Export plots ######
    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_nb_after_filter.png"))
    message("Writing ", plot_path, " plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_numbers(se_after_filter) + labs(title = paste0(prot_or_pep_2, " per sample"), x = "",
         y = paste0("Number of ", prot_or_pep_3)))
    dev.off()


    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_miss_map_after_filter.png"))
    message("Writing ", plot_path, " plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_missval(se_after_filter))
    dev.off()


    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_normalize.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    temp_after_filter = se_after_filter
    temp_after_filter@assays@data[[1]] = log2(temp_after_filter@assays@data[[1]])
    print(plot_normalization(temp_after_filter, se_after_normalization, se_after_imputation))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_normalize_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_normalization(se_after_filter))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_normalize_normalized_vs_imputation.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(plot_normalization(se_after_normalization, se_after_imputation))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions.png"))
    message("Writing ", plot_path, " plot")
    png(filename= plot_path, width = 1000, height = 800)
    print(plot_imputation(temp_after_filter, se_after_normalization, se_after_imputation))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions_after_filter.png"))
    message("Writing ", plot_path, " plot")
    png(filename= plot_path, width = 1000, height = 800)
    print(plot_imputation(se_after_filter))
    dev.off()

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_conditions_distributions_normalized_vs_imputation.png"))
    message("Writing ", plot_path, " plot")
    png(filename= plot_path, width = 1000, height = 800)
    print(plot_imputation(se_after_normalization, se_after_imputation))
    dev.off()
    
    ## Missing values plot
    tab_after <- as.data.frame(se_after_filter@assays@data[[1]])

    # Proportion of missing data by sample
    missing.values <- tab_after %>%
            gather(key = "key", value = "val") %>%
            mutate(isna = is.na(val)) %>%
            group_by(key) %>%
            mutate(total = n()) %>%
            group_by(key, total, isna) %>%
            summarise(num.isna = n()) %>%
            mutate(pct = num.isna / total * 100)

    levels <- (missing.values  %>% filter(isna == T) %>% arrange(desc(pct)))$key

    percentage.plot <- missing.values %>%
        ggplot() + geom_bar(aes(x = key, y = pct, fill=isna), stat = 'identity', alpha=0.8) +
        scale_x_discrete(limits = levels) +
        scale_fill_manual(name = "",values = c('darkslateblue', 'gray53'), labels = c("Present", "Missing")) +
        coord_flip() +
        labs(title = "Ratio of missing value by sample", x =' Samples', y = "% of missing value")

    plot_path <- file.path(output_dir, paste0(prot_or_pep,"_missing_values_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(percentage.plot)
    dev.off()

    row.plot <- tab_after %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
    values = c('darkslateblue', 'gray53'),
    labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Samples",
    y = prot_or_pep_2, title = paste0("Missing ",prot_or_pep_3," by sample, after filter")) +
    coord_flip()

    plot_path <- file.path(output_dir, paste0(prot_or_pep, "_missing_", prot_or_pep_4, "_by_sample_map_after_filter.png"))
    message("Writing ",plot_path," plot")
    png(filename=plot_path, width = 1000, height = 800)
    print(row.plot)
    dev.off()

    
    if (prot_or_pep == "prot"){
        # computations for prot_correlation_plot_filter.png and prot_distribution_plot_filter.png
        if(length(tab_after[,1:length(colnames(tab_after))])!=1){
            mediane <- apply(tab_after[,1:length(colnames(tab_after))],1,median,na.rm=TRUE)
            nb.categories <- trunc(log10(dim(tab_after[,1:length(colnames(tab_after))][,1:length(colnames(tab_after))])[1]))
            stats.matrix <- matrix(nrow=1+2*(nb.categories+1),ncol=dim(tab_after[,1:length(colnames(tab_after))])[2])
            colnames(stats.matrix) <- colnames(tab_after[,1:length(colnames(tab_after))])
            intervals <- 0:nb.categories
            petites.valeurs <- paste("pv-",as.character(10^intervals),sep="")
            grandes.valeurs <- paste("gv-",as.character(10^rev(intervals)),sep="")
            rownames(stats.matrix) <- c("coeff.corr",petites.valeurs,grandes.valeurs)
            for(i in (1:dim(tab_after[,1:length(colnames(tab_after))])[2])){
                stats.matrix[1,i] <- cor(tab_after[,1:length(colnames(tab_after))][,1:length(colnames(tab_after))][,i],mediane,use="complete.obs")
                mat <- sort(tab_after[,1:length(colnames(tab_after))][,i],na.last=NA)
                for(j in intervals){
                    stats.matrix[j+2,i] <- mat[10^j]
                    stats.matrix[dim(stats.matrix)[1]-j,i] <- mat[length(mat)-10^j]
                }
            }
        }

        p_75 <- tab_after
        plot_path <- file.path(output_dir, paste0(prot_or_pep,"_correlation_plot_filter.png"))
        message("Writing ",plot_path," plot")
        png(file=plot_path, width = 1000, height = 800)
        correlation_plot_filtre <- plot( stats.matrix[1,], ylab="correlation coeff.", ylim = c(0,1), xlab="samples", lwd=1, col="blue", type="l", xaxt="n", main="Correlation to median sample expression" )
        axis(side=1, labels = colnames( p_75[,1:length(colnames(tab_after))]), at = 1:length(colnames(tab_after)), cex.axis =  0.55, las = 2 )
        print(correlation_plot_filtre)
        dev.off()

        plot_path <- file.path(output_dir, paste0(prot_or_pep,"_distribution_plot_filter.png"))
        message("Writing ",plot_path," plot")
        png(file=plot_path, width = 1000, height = 800)
        min_born <- floor(min(stats.matrix[2,])) -5
        max_born <- ceiling(max(stats.matrix[9,]))+20
        distrib_plot_filtre <- plot(stats.matrix[2,], ylab="Expression intensity", ylim = c(min_born, max_born), xlab="Samples", lwd=1, col="blue", type="l", axes=FALSE, main=paste0(prot_or_pep_5, " expression quantiles"))
        axis(2)
        lines(stats.matrix[3,], col="Red" )
        lines(stats.matrix[4,], col="Darkgreen" )
        lines(stats.matrix[5,], col="Green" )
        lines(stats.matrix[6,], col="Orange" )
        lines(stats.matrix[7,], col="Grey" )
        lines(stats.matrix[8,], col="Black" )
        lines(stats.matrix[9,], col="Purple" )
        legend(1, max_born, legend=c("PV", "PV10", "PV100", "PV1000", "PG1000", "PG100", "PG10", "PG"), col=c("Blue", "Red", "Darkgreen", "Green", "Orange", "Grey","Black","Purple"), lty = c(1))
        axis(side=1, labels = colnames(p_75[,1:length(colnames(tab_after))]), at = 1:length(colnames(tab_after)), cex.axis = 0.75, las=2)
        print(distrib_plot_filtre)
        dev.off()
    }
    if (prot_or_pep == "pep"){
        message("Not creating ", plot_path, " because processing peptides")
    }
    if (prot_or_pep == "both"){
        message("Not creating ", plot_path, " because processing proteins + peptides")
    }
}