rf_imput <- function(df, ntree=100, k = 10, rowmax = 0.4, colmax = 0.8, my_prefix="") {
    message("This function imputates missing data with knn impute then Random Forest, and compute some accuracy metrics.")
    message("")
    message("    Input:")
    message("        - df, must be a dataframe. In case of proteome array: row should be samples and columns should be proteins.")
    message("        - ntree, optional, number of tree for the random forest. Default is 100")
    message("        - k, optional, number of closest neighboors for the impute.knn. Default is 10")
    message("        - rowmax, optional, rowmax argument of impute.knn. Default is 0.4")
    message("        - colmax, optional, colmax argument of impute.knn. Default is 0.8")
    message("")
    message("    Output:")
    message("        - $data, the knn impute + random forest imputed dataframe")
    message("        - $full_rf_imputation, the dataframe with all values are form the knn+rf imputation (not only the missing values)")
    message("        - $knn_imput, the knn imputed dataframe")
    message("        - $quantile_mape, the 'full' quantile_mape")
    message("        - $no_overfitting_quantile_mape, the 'no overfitting' quantile_mape")
    message("        - $mae, the 'full' mae")
    message("        - $no_overfitting_mae, the 'no overfitting' mae")
    message("        - $rmse, the 'full' rmse")
    message("        - $no_overfitting_rmse, the 'no overfitting' rmse")
    message("")
    message("")
    
    message("Installing necessary packages and libraries")

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    if (!requireNamespace("randomForest", quietly = TRUE))
        install.packages("randomForest")
    
    if (!requireNamespace("itertools", quietly = TRUE))
        install.packages("itertools")
    
    if (!requireNamespace("Metrics", quietly = TRUE))
        install.packages("Metrics")
    
#     BiocManager::install("impute", force=TRUE)
    library(impute)
    library(randomForest)
    library(itertools)
    library(Metrics)
    
    message("End of packages and libraries installation")
    
    message("")
    message("")

    df_columns = colnames(df)
    df_string_int_cols = c()
    for (i in 1:length(df)){
        df_string_int_cols = append(df_string_int_cols,paste("X",toString(i),sep=""))
    }
    colnames(df) = df_string_int_cols
    
    df_knn_imput <- data.frame(df)
    df_rf_imput <- data.frame(df)
    df_rf_full_imput <- data.frame(df)

    mat <- as.matrix(df_knn_imput)
    mat <- ifelse( mat == 'NA', NA, mat )
    mat <- ifelse( mat == '', NA, mat )

    message("Starting the knn imputation (Random Forest needs a first naive imputation method)")
    imputed_mat <- impute.knn(mat,k = k, rowmax = rowmax, colmax = colmax)
    message("End of knn imputation")
    df_knn_imput <- as.data.frame(imputed_mat$data)
        
    message("")
    message("")
        
    message("Starting the random forest imputation")
    count <- 0
    col_number = length(colnames(df_knn_imput))

	message(colnames(df)[1:5])

    for (col in colnames(df)){
        count <- count + 1
        if (count == 1 || count == 2 || count == 10 || count %% 100 == 0){
            message(count, " on ",col_number," columns")
        }

        #message(col)
        #message(df[,col])
        #message(as.matrix(df[col])[1:10])

        if (NA %in% as.matrix(df[,col])){
            new_formula <- paste(col," ~.")

            rf <- randomForest(as.formula(paste(col," ~.")), data = df_knn_imput, ntree = ntree)
            a_predictions <- predict(rf, df_knn_imput)
            df_rf_full_imput[,col] <- a_predictions

            ind_count <- 1
            for (ind in as.matrix(df[,col])){
                if (is.na(ind)){
                    df_rf_imput[ind_count,col] <- a_predictions[ind_count]
                }
                ind_count <- ind_count + 1
            }
        } else {
            rf <- randomForest(as.formula(paste(col," ~.")), data = df_knn_imput, ntree = ntree)
            a_predictions <- predict(rf, df_knn_imput)
            df_rf_full_imput[,col] <- a_predictions
        }
    }
    message("End of the random forest imputation")
    message(col_number, " on ",col_number)
    message("")
    message("")
    message("Computing the 95% range interval of values")
    it <- ihasNext(product(a=rownames(df_rf_imput), b=colnames(df_rf_imput)))
    a_expr <- c()
    points_number = length(rownames(df_rf_imput)) * length(colnames(df_rf_imput))
    count = 0
    while (hasNext(it)) {
        x <- nextElem(it)
        count = count + 1
        if (count == 1 || count == 2 || count == 10 || count ==100  || count == 1000 || count %% 10000 == 0){
            message(count, " on ",points_number)
        }
        if ( !is.na(df[x$a, x$b]) ){
            a_expr = append(a_expr, df[x$a, x$b])
        }
    }

    message(points_number, " on ",points_number)

    quant_025 <- unname(quantile(a_expr,0.025))
    quant_975 <- unname(quantile(a_expr,0.975))
    interval_size <- quant_975 - quant_025

    message("")
    message("2.5% quantile is ", quant_025)
    message("97.5% quantile is ", quant_975)
    message("95% range interval value is ", interval_size)
    message("")
    message("")
    message("Let's make some accuracy metric computation")
    message("")
    message("")
    it <- ihasNext(product(a=rownames(df_rf_imput), b=colnames(df_rf_imput)))
    a_quantile_mapes = c()
    a_non_zero_quantile_mapes = c()
    a_trues = c()
    a_preds = c()
    a_non_zero_trues = c()
    a_non_zero_preds = c()
    count = 0

    points_number = length(rownames(df_rf_imput)) * length(colnames(df_rf_imput))

    while (hasNext(it)) {
        x <- nextElem(it)
        count = count + 1
        if (count == 1 || count == 2 || count == 10 || count == 100  || count == 1000 || count %% 10000 == 0){
            message(count, " on ",points_number)
        }
        if ( !is.na(df[x$a, x$b]) ){
            a_quantile_mapes <- append(a_quantile_mapes,abs(df[x$a, x$b] -  df_rf_full_imput[x$a, x$b])/interval_size)
            if (df[x$a, x$b]!= df_rf_full_imput[x$a, x$b]){
                a_non_zero_quantile_mapes <- append(a_non_zero_quantile_mapes, abs(df[x$a, x$b] -  df_rf_full_imput[x$a, x$b])/interval_size)
                a_non_zero_trues <- append(a_non_zero_trues, df[x$a, x$b])
                a_non_zero_preds <- append(a_non_zero_preds, df_rf_full_imput[x$a, x$b])
            }
            a_trues <- append(a_trues, df[x$a, x$b]) 
            a_preds <- append(a_preds, df_rf_full_imput[x$a, x$b])
        }
    }
    message(points_number, " on ",points_number)

    message("")
    message("To compensate overfitting we consider for each metric:")
    message("1. 'full' metric is computed on all available data ")
    message("2. 'no overfitting' metric is only computed when actual value != predicted value")
    message("")
    message("                 If 'full' is very different from 'no overfitting', it means the model overfitted.")
    message("                 In case of overfitting, accuracy metrics are probably misleading, but the imputation method still stays valid.")
    message("")

    mean_quantile_mape = mean(a_quantile_mapes)
    message("full quantile_mape: ", mean_quantile_mape)
    no_zero_mean_quantile_mape = mean(a_non_zero_quantile_mapes)
    message("no overfitting quantile_mape: ", no_zero_mean_quantile_mape)
    classic_mae = mae(a_trues, a_preds)
    message("full mae: ", classic_mae)
    non_zero_mae = mae(a_non_zero_trues, a_non_zero_preds)
    message("no overfitting mae: ", non_zero_mae)
    classic_rmse = mae(a_trues, a_preds)
    message("full rmse: ", classic_rmse)
    non_zero_rmse = mae(a_non_zero_trues, a_non_zero_preds)
    message("no overfitting rmse: ", non_zero_rmse)

    colnames(df_rf_imput) = df_columns
    colnames(df_rf_full_imput) = df_columns
    colnames(df_knn_imput) = df_columns
    
  return(list("data"=df_rf_imput,
              "full_rf_imputation"=df_rf_full_imput,
              "knn_imput"=df_knn_imput,
              "quantile_mape"=mean_quantile_mape,
              "no_overfitting_quantile_mape"=no_zero_mean_quantile_mape,
              "mae"=classic_mae,
              "no_overfitting_mae"=non_zero_mae,
              "rmse"=classic_rmse,
              "no_overfitting_rmse"=non_zero_rmse
             ))
}
