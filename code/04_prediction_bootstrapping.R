pROC <- function(continuous_mod,test_data,
                 n_iter=1000,E_percent=5,
                 boost_percent=50,
                 parallel=FALSE,ncores=4,rseed=FALSE,
                 sub_sample=FALSE,sub_sample_size=10000){
  
  if (methods::is(continuous_mod,"RasterLayer")) {
    if (continuous_mod@data@min == continuous_mod@data@max) {
      stop("\nModel with no variability.\n")
    }
    if (is.data.frame(test_data) || is.matrix(test_data)) {
      test_data <- stats::na.omit(raster::extract(continuous_mod,
                                                  test_data))
      
    }
    vals <- continuous_mod[!is.na(continuous_mod[])]
  }
  if(is.numeric(continuous_mod)){
    vals <- continuous_mod
    if (!is.numeric(test_data))
      stop("If continuous_mod is of class numeric,
           test_data must be numeric...")
  }
  if(sub_sample){
    nvals <- length(vals)
    if(sub_sample_size> nvals) sub_sample_size <- nvals
    vals <- base::sample(vals,size = sub_sample_size)
  }
  ndigits <- proc_precision(mod_vals = vals,
                            test_data = test_data)
  
  tomult <- as.numeric(paste0("1e+",ndigits))
  test_value <- test_data*tomult
  test_value <- round(as.vector(test_value))
  
  vals2 <- round(vals*tomult)
  classpixels <- as.data.frame(base::table(vals2),
                               stringsAsFactors = F)
  names(classpixels) <- c("value",
                          "count")
  classpixels$value <- as.numeric(classpixels$value)
  classpixels <- data.frame(stats::na.omit(classpixels))
  value <- count <- totpixperclass <- NULL
  classpixels <- classpixels %>%
    dplyr::mutate(value  = rev(value),
                  count = rev(count),
                  totpixperclass = cumsum(count),
                  percentpixels = totpixperclass/sum(count)) %>%
    dplyr::arrange(value)
  
  #if(nrow(classpixels)>1500){
  #  classpixels <- classpixels %>%
  #    dplyr::sample_n(1500) %>% dplyr::arrange(value)
  #}
  
  error_sens <- 1 - (E_percent/100)
  models_thresholds <- classpixels[, "value"]
  fractional_area <- classpixels[, "percentpixels"]
  n_data <- length(test_value)
  n_samp <- ceiling((boost_percent/100) * (n_data))
  
  big_classpixels <- matrix(rep(models_thresholds,
                                each = n_samp),
                            ncol = length(models_thresholds))
  
  
  calc_aucDF <- function(big_classpixels,
                         fractional_area,
                         test_value,
                         n_data, n_samp,
                         error_sens,rseed=NULL) {
    if(is.numeric(rseed)) set.seed(rseed)
    trapz_roc <- function(x,y){
      size_x <- length(x)
      xp <- c(x, x[size_x:1])
      yp <- c(numeric(size_x), y[size_x:1])
      nda <- 2 * size_x
      p1 <- sum(xp[1:(nda - 1)] * yp[2:nda]) + xp[nda] * yp[1]
      p2 <- sum(xp[2:nda] * yp[1:(nda - 1)]) + xp[1] * yp[nda]
      return(0.5 * (p1 - p2))
    }
    
    rowsID <- sample(x = n_data,
                     size = n_samp,
                     replace = TRUE)
    
    test_value1 <- test_value[rowsID]
    omssion_matrix <- big_classpixels > test_value1
    sensibility <- 1 - colSums(omssion_matrix)/n_samp
    xyTable <- data.frame(fractional_area, sensibility)
    xyTable <- rbind(xyTable,c(0,0))
    xyTable <- xyTable[order(xyTable$fractional_area,
                             decreasing = F),]
    auc_model <- try(trapz_roc(xyTable$fractional_area,
                               xyTable$sensibility),silent = TRUE)
    if(methods::is(try(auc_model), "try-error")) auc_model <- NA
    
    
    if(error_sens>0){
      less_ID <- which(xyTable$sensibility <= error_sens)
      xyTable <- xyTable[-less_ID, ]
      auc_pmodel <- try(trapz_roc(xyTable$fractional_area,
                                  xyTable$sensibility),silent = TRUE)
      
      if(methods::is(try(auc_pmodel),"try-error")) auc_pmodel <- NA
      
      auc_prand <- try(trapz_roc(xyTable$fractional_area,
                                 xyTable$fractional_area),silent = TRUE)
      if(methods::is(try(auc_prand), "try-error")) auc_prand <- NA
      
    }
    else{
      auc_pmodel <- auc_model
      auc_prand <- 0.5
      
    }
    
    
    
    auc_ratio <- auc_pmodel/auc_prand
    auc_table <- data.frame(auc_model,
                            auc_pmodel,
                            auc_prand,
                            auc_ratio)
    return(auc_table)
  }
  
  if (parallel) {
    n_cores <- ntbox::nc(ncores)
    
    #furrr::furrr_options(packages = c("Rcpp","ntbox"))
    plan(multisession,workers=n_cores)
    options(future.globals.maxSize= 8500*1024^2)
    
    roc_env <- new.env()
    niter_big <- floor(n_iter/n_cores)
    n_runs <- rep(niter_big, n_cores)
    sum_n_runs <- sum(n_runs)
    n_runs[1] <- n_runs[1] + (n_iter - sum_n_runs)
    
    for (i in 1:length(n_runs)) {
      x <- as.character(i)
      roc_env[[x]] %<-% {
        x1 <- 1:n_runs[i]
        auc_matrix1 <- x1 %>%
          purrr::map_df(~calc_aucDF(big_classpixels,
                                    fractional_area,
                                    test_value,
                                    n_data, n_samp,
                                    error_sens,rseed=NULL))
      }
    }
    partial_AUC <- as.list(roc_env)
    rm(roc_env)
    partial_AUC <- do.call(rbind.data.frame,
                           partial_AUC)
    rownames(partial_AUC) <- NULL
    future::plan(future::sequential)
  }
  else {
    partial_AUC <- 1:n_iter %>%
      purrr::map_df(function(i){
        proc <- calc_aucDF(big_classpixels,
                           fractional_area,
                           test_value,
                           n_data,
                           n_samp,
                           error_sens,rseed = i)
      })
  }
  mauc <-  mean(partial_AUC$auc_model, na.rm = TRUE)
  maucp <- mean(partial_AUC$auc_ratio, na.rm = TRUE)
  proc <- sum(partial_AUC$auc_ratio <= 1, na.rm = TRUE)/
    length(partial_AUC$auc_ratio[!is.na(partial_AUC$auc_ratio)])
  
  p_roc <- c(mauc,maucp, proc)
  names(p_roc) <- c("Mean_AUC",
                    paste("Mean_pAUC_ratio_at_",
                          E_percent,
                          "%", sep = ""),
                    "P_value")
  p_roc_res <- list(pROC_summary = p_roc,
                    pROC_results = partial_AUC)
  return(p_roc_res)
}


proc_precision <- function(mod_vals,test_data){
  
  min_vals <- min(mod_vals,na.rm = TRUE)
  
  percentil_test <- unique(sort(stats::na.omit(test_data)))[2]
  
  
  #percentil_test <- stats::quantile(test_data,
  #                                  probs=0.1)
  partition_flag <- mean(c(min_vals,
                           percentil_test))
  fflag <- stringr::str_detect(partition_flag, "e")
  if (length(fflag)>0L && fflag) {
    ndigits <- stringr::str_split(partition_flag, "e-")[[1]]
    ndigits <- as.numeric(ndigits)[2] #- 1
  }
  else {
    med <- stringr::str_extract_all(partition_flag, pattern = "[0-9]|[.]")
    med <- unlist(med)
    med <- med[-(1:which(med == "."))]
    med1 <- which(med != 0)
    ndigits <- ifelse(med1[1] <= 2, 3, 4)
  }
  return(ndigits)
}

random_forest_function <- function(data, out_dir, no_iter, split_ratio){
  library(enmSdmX)
  library(randomForest)
  library(caret)
  require(caTools)
  library(pdp)
  library(metrica)
  library(blockCV)
  library(vip)
  library(fastshap)
  library(spatialsample)
  library(sf)
  library(pROC)
  
  #create storage for all performance measures
  auc <- c()
  auc_list <- c() 
  pauc_mean_list <- c()
  pauc_ratio_list <- c()
  pauc_pval_list <- c()
  metrics_subregion_df_iter <- c() 
  fscore_list <- c()
  sensitivity_list <- c()
  specificity_list <- c()
  tss_list <- c()
  oob_list <- c()
  imp_df_iter <- c() 
  pd_df_iter <- c()
  pred_storage_iter <- c()
  pred_storage_subregion_iter <- c()
  stor.data <- data
  
  data <- subset(data, select=-c(geometry))#remove geom col
  
  for(i in 1:no_iter){
    set.seed(i*989) #set diff seed number for each model iteration, since this will split your data in different test-train sets / set up how the tree is split each time
    #data <- stor.data #reset data variable for new subset choice
    
    #randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
    rows <- sample(nrow(data))
    data <- data[rows,]
    
    #could add a bootstrapping step here by running code below with a subset of data
    sample = sample.split(data$presence, SplitRatio = split_ratio)
    train = subset(data, sample == TRUE)
    #test  = subset(data, sample == FALSE)
    
    #remove missing data
    train_complete <- train[complete.cases(train),]; train_complete$presence <- as.factor(train_complete$presence)
    #test_complete <- test[complete.cases(test),]; test_complete$presence <- as.factor(test_complete$presence)
    
    # sdm_rf <- trainRF(data = subset(train_complete, select=-c(subregion)),
    #                   resp = 'presence',
    #                   preds = 2:ncol(subset(train_complete, select=-c(subregion))))
    models_dir <- paste(out_dir,'/models',sep="")
    models <- list.files(path=models_dir, pattern="rf", all.files=FALSE, full.names=TRUE,recursive=TRUE)
    sdm_rf <- readRDS(file=models[i])
    pred_train <- as.data.frame(predict(sdm_rf, newdata=subset(train_complete, select=-c(presence,subregion)), 'prob')) #take out response and subregion column from test data
    
    # Store predictions for ensemble model evaluation
    train_complete[,1] <- as.numeric(as.character(train_complete[,1]))
    pred_storage_train <- as.data.frame(cbind(rownames(train_complete), as.numeric(pred_train[,2]), train_complete[,1]))
    colnames(pred_storage_train) <- c("row_id", "predictions", "true_values")
    
    # Calculate model performance measures for *main model*
    auc <- pROC::roc(response=as.numeric(train_complete[,1]), predictor=as.numeric(pred_train[,2]), auc = TRUE)
    print(auc$auc)
    #pROC::plot.roc(roc(response=as.numeric(test_complete[,1]), predictor=as.numeric(pred[,2]), levels=c(1,2)))
    ## Add pROC
    # Calculate prediction vectors needed for pROC/AUC
    pred_storage_prescence_train <- pred_storage_train[which(as.numeric(pred_storage_train$true_values)==1),]
    # Calculate partial roc
    partial_roc <- pROC(continuous_mod=as.numeric(pred_storage_train$predictions),
                        test_data = as.numeric(pred_storage_prescence_train$predictions),
                        n_iter=200,E_percent=100,
                        boost_percent=50,
                        parallel=FALSE)
    # Store all pauc information
    pauc_mean <- as.numeric(partial_roc$pROC_summary[1])
    pauc_ratio <- as.numeric(partial_roc$pROC_summary[2])
    pauc_pval <- as.numeric(partial_roc$pROC_summary[3])
    print(pauc_mean)
    pauc_bootstrap_results <- partial_roc$pROC_results
    # Calculate best threshold for typical auc calculation to use with metrica package
    best.threshold <- pROC::coords(auc, "best", ret = "threshold")
    # Use best threshold to build 1/0 prediction dataset (instead of continuous)
    metrica.format <- data.frame(cbind(ifelse(train_complete[,1]==1,1,0)),ifelse(as.numeric(pred_train[,2])>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
    fscore <- metrica::fscore(data = metrica.format, obs = labels, pred = predictions)$fscore
    sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
    specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
    tss <- sensitivity + specificity - 1
    oob_error <- sdm_rf$err.rate[500,1] #first value is # of trees
    
    # Store models to build full prediction rasters later on, remember need to build "models" folder in output folder!
    #saveRDS(sdm_rf, paste0(out_dir,"/models/rf_model", i, ".rda"))
    print(i)
    
    auc_list <- c(auc_list, auc$auc)
    pauc_mean_list <- c(pauc_mean_list, pauc_mean)
    pauc_ratio_list <- c(pauc_ratio_list, pauc_ratio)
    pauc_pval_list <- c(pauc_pval_list, pauc_pval)
    fscore_list <- c(fscore_list, fscore)
    sensitivity_list <- c(sensitivity_list, sensitivity)
    specificity_list <- c(specificity_list, specificity)
    tss_list <- c(tss_list, tss)
    oob_list <-c(oob_list, oob_error) #store oob error rate
    
    #Start writing / storing output!
    metrics.output.dataframe <- cbind(auc_list, pauc_mean_list, pauc_ratio_list, pauc_pval_list, fscore_list, sensitivity_list, specificity_list, tss_list); colnames(metrics.output.dataframe) <- c("auc", "pauc_mean", "pauc_ratio", "pauc_pval", "fscore", "sensitivity", "specificity", "tss")
    
    # Subregion measures
    if("subregion" %in% colnames(data)){
      # Create storage dataframe for seven performance measures to fill in for each subregion.
      metrics_subregion_df = data.frame(matrix(vector(), 0, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))),
                                        row.names = NULL, stringsAsFactors=F)
      pred_storage_subregion <- c()
      
      for (k in 1:length(unique(data$subregion))) { #loop through each subregion
        # Build subregion test set (part of larger test set that is in subregion k)
        train_subregion <- train_complete[which(train_complete$subregion==as.character(k)),]
        train_subregion <- subset(train_subregion, select=-c(subregion))
         
        # Create storage dataframes to get overwritten during each iteration of the loop
        subregion_loop_df <- data.frame(matrix(vector(), 1, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))), 
                                        stringsAsFactors=F, row.names=NULL)
        pred_storage_subregion_loop <- c()
        
        # Build prediction dataframe for subregion k and store if has values, otherwise feed NAs
        if (dim(train_subregion)[1] != 0) {
          pred.subregion <- as.data.frame(predict(sdm_rf, newdata=subset(train_subregion, select=-c(presence)), 'prob')) #take out response from test data
          pred.subregion <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,2]), train_subregion[,1]))
          colnames(pred.subregion) <- c("row_id", "predictions", "true_values")
          # Store predictions for ensemble model evaluation
          pred_storage_subregion_loop <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,2]), train_subregion[,1]))
          pred_storage_subregion_loop$subregion <- k
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }else{
          pred_storage_subregion_loop <- as.data.frame(cbind(c(NA),c(NA),c(NA),k))
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }
        
        # Check to make sure subregion has both an occurrence point and background point, otherwise cannot calculate AUC
        if(1 %in% train_subregion$presence && 0 %in% train_subregion$presence) {
          # Calculate metrics
          auc.subregion <- pROC::roc(response=as.numeric(train_subregion[,1]), predictor=as.numeric(pred.subregion[,2]), auc = TRUE) 
          pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
          #print(pred.subregion)
          #print(pred_subregion_prescence_train)
          if(dim(pred_subregion_prescence_train)[1]>2){#length(pred.subregion$predictions)>2 && 
            ## Add pROC
            # Calculate prediction vectors needed for pROC/AUC
            pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
            partial_roc_subregion <- pROC(continuous_mod=as.numeric(pred.subregion$predictions),
                                test_data = as.numeric(pred_subregion_prescence_train$predictions),
                                n_iter=200,E_percent=100,
                                boost_percent=50,
                                parallel=FALSE)
            # Store all pauc information
            pauc_mean_subregion <- as.numeric(partial_roc_subregion$pROC_summary[1])
            pauc_ratio_subregion <- as.numeric(partial_roc_subregion$pROC_summary[2])
            pauc_pval_subregion <- as.numeric(partial_roc_subregion$pROC_summary[3])
            pauc_bootstrap_results_subregion <- partial_roc_subregion$pROC_results
          }else{
            pauc_mean_subregion <- NA
            pauc_ratio_subregion <- NA
            pauc_pval_subregion <- NA
            pauc_bootstrap_results_subregion <- NA
          }
          best.threshold.subregion <- pROC::coords(auc.subregion, "best", ret = "threshold")
          metrica.format.subregion <- data.frame(cbind(ifelse(train_subregion[,1]==1,1,0)),ifelse(as.numeric(pred.subregion[,2])>=best.threshold.subregion[1,1],1,0)); colnames(metrica.format.subregion) <- c("labels","predictions"); rownames(metrica.format.subregion) <- 1:dim(metrica.format.subregion)[1]
          # Additional check to make sure subregion has both an occurrence point PREDICTED and a background point PREDICTED, otherwise cannot calculate metrica measures
          if(1 %in% metrica.format.subregion$predictions && 0 %in% metrica.format.subregion$predictions){
            fscore.subregion <- metrica::fscore(data = metrica.format.subregion, obs = labels, pred = predictions)$fscore
            sensitivity.subregion <- metrica::recall(data = metrica.format.subregion, obs = labels, pred = predictions)$recall 
            specificity.subregion <- metrica::specificity(data = metrica.format.subregion, obs = labels, pred = predictions)$spec
            tss.subregion <- sensitivity.subregion + specificity.subregion - 1
          }else{ #If second condition not met fill everything with NAs besides AUC
            fscore.subregion <- NA
            sensitivity.subregion <- NA
            specificity.subregion <- NA
            tss.subregion <- NA
          }
          
          #could probably compile into above once working well
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1] #print number of test points predicted on to get a sense of data quantity for each subregion (should be on average ~20% of total amount)
          subregion_loop_df$auc <- auc.subregion$auc
          subregion_loop_df$pauc_mean <- pauc_mean_subregion
          subregion_loop_df$pauc_ratio <- pauc_ratio_subregion
          subregion_loop_df$pauc_pval <- pauc_pval_subregion
          subregion_loop_df$fscore <- fscore.subregion
          subregion_loop_df$sensitivity <- sensitivity.subregion
          subregion_loop_df$specificity <- specificity.subregion
          subregion_loop_df$tss <- tss.subregion
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        } else { #If first condition not met fill everything with NAs (including AUC)
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1]
          subregion_loop_df$auc <- NA
          subregion_loop_df$pauc_mean <- NA
          subregion_loop_df$pauc_ratio <- NA
          subregion_loop_df$pauc_pval <- NA
          subregion_loop_df$fscore <- NA
          subregion_loop_df$sensitivity <- NA
          subregion_loop_df$specificity <- NA
          subregion_loop_df$tss <- NA
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        }
      }
      metrics_subregion_df$iter <- i
      metrics_subregion_df_iter <- rbind(metrics_subregion_df_iter, metrics_subregion_df)
      pred_storage_subregion$iter <- i
      pred_storage_subregion_iter <- rbind(pred_storage_subregion_iter, pred_storage_subregion)
    }
    
  }
  write.csv(metrics.output.dataframe, paste0(out_dir, "/rf_baseline_metrics_cropped.csv"))
  write.csv(oob_list, paste0(out_dir, "/rf_oob_cropped.csv"))
  write.csv(metrics_subregion_df_iter, paste0(out_dir, "/rf_subregion_metrics_cropped.csv"))
}
brt_function <- function(data, out_dir, no_iter, split_ratio){
  library(enmSdmX)
  library(randomForest)
  library(caret)
  require(caTools)
  library(pdp)
  library(metrica)
  library(blockCV)
  library(vip)
  library(fastshap)
  library(spatialsample)
  library(sf)
  library(pROC)
  
  #create storage for all performance measures
  auc <- c()
  auc_list <- c() 
  metrics_subregion_df_iter <- c() 
  pauc_mean_list <- c()
  pauc_ratio_list <- c()
  pauc_pval_list <- c()
  fscore_list <- c()
  sensitivity_list <- c()
  specificity_list <- c()
  tss_list <- c()
  oob_list <- c()
  imp_df_iter <- c() 
  pd_df_iter <- c()
  pred_storage_iter <- c()
  pred_storage_subregion_iter <- c()
  stor.data <- data
  
  data <- subset(data, select=-c(geometry))#remove geom col
  
  for(i in 1:no_iter){
    set.seed(i*989) #set diff seed number for each model iteration, since this will split your data in different test-train sets / set up how the tree is split each time
    #data <- stor.data
    
    #randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
    rows <- sample(nrow(data))
    data <- data[rows,]
    
    #add a bootstrapping step here by running code below with a subset of data
    sample = sample.split(data$presence, SplitRatio = split_ratio)
    train = subset(data, sample == TRUE)
    train <- as.data.frame(train)
    #remove missing data
    train_complete <- train[complete.cases(train),]#; train_complete$presence <- as.factor(train_complete$presence)
    #(1) train model OF (below)
    # sdm_brt <- trainBRT(data = subset(train_complete, select=-c(subregion)),
    #                     resp = 'presence',
    #                     preds = names(subset(train_complete, select=-c(presence,subregion))),
    #                     w = FALSE,
    #                     out = c("model")#,
    #                     #minTrees = 800,
    #                     #maxTrees = 14000,
    #                     #verbose = TRUE,
    #                     #treeComplexity = c(1:16)
    # )
    #(2) load trained models
    models_dir <- paste(out_dir,'/models',sep="")
    models <- list.files(path=models_dir, pattern="brt", all.files=FALSE, full.names=TRUE,recursive=TRUE)
    sdm_brt <- readRDS(file=models[i])
    if (is.null(sdm_brt)) { #test to check model converged, pass to next iteration of for loop if not
      next
    }
    pred_train <- as.data.frame(predict(sdm_brt, newdata=subset(train_complete, select=-c(presence,subregion)), type = 'response')) #take out response and subregion column from test data
    pred_storage_train <- as.data.frame(cbind(rownames(train_complete), as.numeric(pred_train[,1]), train_complete[,1]))
    colnames(pred_storage_train) <- c("row_id", "predictions", "true_values")
    
    auc <- pROC::roc(response=as.numeric(train_complete[,1]), predictor=as.numeric(pred_train[,1]), auc = TRUE) 
    #print(auc$auc)
    ## Add pROC
    # Calculate prediction vectors needed for pROC/AUC
    pred_storage_prescence_train <- pred_storage_train[which(as.numeric(pred_storage_train$true_values)==1),]
    # Calculate partial roc
    partial_roc <- pROC(continuous_mod=as.numeric(pred_storage_train$predictions),
                        test_data = as.numeric(pred_storage_prescence_train$predictions),
                        n_iter=200,E_percent=100,#-(100*as.numeric(minimum)),
                        boost_percent=50,
                        parallel=FALSE)
    #print(partial_roc$pROC_summary)
    #print(head(partial_roc$pROC_results))
    # Store all pauc information
    pauc_mean <- as.numeric(partial_roc$pROC_summary[1])
    pauc_ratio <- as.numeric(partial_roc$pROC_summary[2])
    pauc_pval <- as.numeric(partial_roc$pROC_summary[3])
    pauc_bootstrap_results <- partial_roc$pROC_results
    # Calculate best threshold for typical auc calculation to use with metrica package
    best.threshold <- pROC::coords(auc, "best", ret = "threshold") # Find best threshold to use with metrica dataset
    # Use best threshold to build 1/0 prediction dataset (rather than just probabilities)
    metrica.format <- data.frame(cbind(ifelse(train_complete[,1]==1,1,0)),ifelse(as.numeric(pred_train[,1])>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
    if(1 %in% metrica.format$predictions && 0 %in% metrica.format$predictions){
      fscore <- metrica::fscore(data = metrica.format, obs = labels, pred = predictions)$fscore
      sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
      specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
      tss <- sensitivity + specificity - 1
      oob_error <- sdm_brt$err.rate[500,1] #first value is # of trees
    }else{
      #auc$auc <- NA
      fscore <- NA
      sensitivity <- NA
      specificity <- NA
      tss <- NA
      oob_error <- sdm_brt$err.rate[500,1] #first value is # of trees
    }
    
    # Store models to build full prediction rasters later on, remember need to build "models" folder in output folder!
    #saveRDS(sdm_brt, paste0(out_dir,"/models/brt_model", i, ".rda"))
    print(i)
    
    auc_list <- c(auc_list, auc$auc)
    pauc_mean_list <- c(pauc_mean_list, pauc_mean)
    pauc_ratio_list <- c(pauc_ratio_list, pauc_ratio)
    pauc_pval_list <- c(pauc_pval_list, pauc_pval)
    fscore_list <- c(fscore_list, fscore)
    sensitivity_list <- c(sensitivity_list, sensitivity)
    specificity_list <- c(specificity_list, specificity)
    tss_list <- c(tss_list, tss)
    oob_list <-c(oob_list, oob_error) #store oob error rate
    
    #Start writing / storing output!
    metrics.output.dataframe <- cbind(auc_list, pauc_mean_list, pauc_ratio_list, pauc_pval_list, fscore_list, sensitivity_list, specificity_list, tss_list); colnames(metrics.output.dataframe) <- c("auc", "pauc_mean", "pauc_ratio", "pauc_pval",
                                                                                                                                                                                                     "fscore", "sensitivity", "specificity", "tss")
    if("subregion" %in% colnames(data)){
      # Create storage dataframe for seven performance measures to fill in for each subregion.
      metrics_subregion_df = data.frame(matrix(vector(), 0, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))),
                                        row.names = NULL, stringsAsFactors=F)
      pred_storage_subregion <- c()
      
      for (k in 1:length(unique(data$subregion))) { #loop through each subregion
        # Build subregion test set (part of larger test set that is in subregion k)
        train_subregion <- train_complete[which(train_complete$subregion==as.character(k)),]
        train_subregion <- subset(train_subregion, select=-c(subregion))
        
        # Create storage dataframes to get overwritten during each iteration of the loop
        subregion_loop_df <- data.frame(matrix(vector(), 1, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))), 
                                        stringsAsFactors=F, row.names=NULL)
        pred_storage_subregion_loop <- c()
        
        # Build prediction dataframe for subregion k and store if has values, otherwise feed NAs
        if (dim(train_subregion)[1] != 0) {
          pred.subregion <- as.data.frame(predict(sdm_brt, newdata=subset(train_subregion, select=-c(presence)), type = 'response')) #take out response from test data
          pred.subregion <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,1]), train_subregion[,1]))
          colnames(pred.subregion) <- c("row_id", "predictions", "true_values")
          # Store predictions for ensemble model evaluation
          pred_storage_subregion_loop <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,1]), train_subregion[,1]))
          pred_storage_subregion_loop$subregion <- k
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }else{
          pred_storage_subregion_loop <- as.data.frame(cbind(c(NA),c(NA),c(NA),k))
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }
        
        # Check to make sure subregion has both an occurrence point and background point, otherwise cannot calculate AUC
        if(1 %in% train_subregion$presence && 0 %in% train_subregion$presence) {
          # Calculate metrics
          auc.subregion <- pROC::roc(response=as.numeric(train_subregion[,1]), predictor=as.numeric(pred.subregion[,1]), auc = TRUE) 
          pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
          #print(pred.subregion)
          #print(pred_subregion_prescence_train)
          if(dim(pred_subregion_prescence_train)[1]>2){#length(pred.subregion$predictions)>2 && 
            ## Add pROC
            # Calculate prediction vectors needed for pROC/AUC
            pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
            partial_roc_subregion <- pROC(continuous_mod=as.numeric(pred.subregion$predictions),
                                          test_data = as.numeric(pred_subregion_prescence_train$predictions),
                                          n_iter=200,E_percent=100,
                                          boost_percent=50,
                                          parallel=FALSE)
            # Store all pauc information
            pauc_mean_subregion <- as.numeric(partial_roc_subregion$pROC_summary[1])
            pauc_ratio_subregion <- as.numeric(partial_roc_subregion$pROC_summary[2])
            pauc_pval_subregion <- as.numeric(partial_roc_subregion$pROC_summary[3])
            pauc_bootstrap_results_subregion <- partial_roc_subregion$pROC_results
          }else{
            pauc_mean_subregion <- NA
            pauc_ratio_subregion <- NA
            pauc_pval_subregion <- NA
            pauc_bootstrap_results_subregion <- NA
          }
          best.threshold.subregion <- pROC::coords(auc.subregion, "best", ret = "threshold")
          metrica.format.subregion <- data.frame(cbind(ifelse(train_subregion[,1]==1,1,0)),ifelse(as.numeric(pred.subregion[,1])>=best.threshold.subregion[1,1],1,0)); colnames(metrica.format.subregion) <- c("labels","predictions"); rownames(metrica.format.subregion) <- 1:dim(metrica.format.subregion)[1]
          # Additional check to make sure subregion has both an occurrence point PREDICTED and a background point PREDICTED, otherwise cannot calculate metrica measures
          if(1 %in% metrica.format.subregion$predictions && 0 %in% metrica.format.subregion$predictions){
            fscore.subregion <- metrica::fscore(data = metrica.format.subregion, obs = labels, pred = predictions)$fscore
            sensitivity.subregion <- metrica::recall(data = metrica.format.subregion, obs = labels, pred = predictions)$recall 
            specificity.subregion <- metrica::specificity(data = metrica.format.subregion, obs = labels, pred = predictions)$spec
            tss.subregion <- sensitivity.subregion + specificity.subregion - 1
          }else{ #If second condition not met fill everything with NAs besides AUC
            fscore.subregion <- NA
            sensitivity.subregion <- NA
            specificity.subregion <- NA
            tss.subregion <- NA
          }
          
          #could probably compile into above once working well
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1] #print number of test points predicted on to get a sense of data quantity for each subregion (should be on average ~20% of total amount)
          subregion_loop_df$auc <- auc.subregion$auc
          subregion_loop_df$pauc_mean <- pauc_mean_subregion
          subregion_loop_df$pauc_ratio <- pauc_ratio_subregion
          subregion_loop_df$pauc_pval <- pauc_pval_subregion
          subregion_loop_df$fscore <- fscore.subregion
          subregion_loop_df$sensitivity <- sensitivity.subregion
          subregion_loop_df$specificity <- specificity.subregion
          subregion_loop_df$tss <- tss.subregion
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        } else { #If first condition not met fill everything with NAs (including AUC)
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1]
          subregion_loop_df$auc <- NA
          subregion_loop_df$pauc_mean <- NA
          subregion_loop_df$pauc_ratio <- NA
          subregion_loop_df$pauc_pval <- NA
          subregion_loop_df$fscore <- NA
          subregion_loop_df$sensitivity <- NA
          subregion_loop_df$specificity <- NA
          subregion_loop_df$tss <- NA
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        }
      }
      metrics_subregion_df$iter <- i
      metrics_subregion_df_iter <- rbind(metrics_subregion_df_iter, metrics_subregion_df)
      pred_storage_subregion$iter <- i
      pred_storage_subregion_iter <- rbind(pred_storage_subregion_iter, pred_storage_subregion)
    }
    
  }
  write.csv(metrics.output.dataframe, paste0(out_dir, "/brt_baseline_metrics_cropped.csv"))
  write.csv(oob_list, paste0(out_dir, "/brt_oob_cropped.csv"))
  write.csv(metrics_subregion_df_iter, paste0(out_dir, "/brt_subregion_metrics_cropped.csv"))
}
maxent_function <- function(data, out_dir, no_iter, split_ratio){
  library(enmSdmX)
  library(randomForest)
  library(caret)
  require(caTools)
  library(pdp)
  library(metrica)
  library(blockCV)
  library(vip)
  library(fastshap)
  library(spatialsample)
  library(ecospat) #for maxent variable importance function
  library(sf)
  library(pROC)
  
  #create storage for all performance measures
  auc <- c()
  auc_list <- c() 
  pauc_mean_list <- c()
  pauc_ratio_list <- c()
  pauc_pval_list <- c()
  metrics_subregion_df_iter <- c() 
  fscore_list <- c()
  sensitivity_list <- c()
  specificity_list <- c()
  tss_list <- c()
  oob_list <- c()
  imp_df_iter <- c() 
  pd_df_iter <- c()
  pred_storage_iter <- c()
  pred_storage_subregion_iter <- c()
  stor.data <- data
  
  data <- subset(data, select=-c(geometry))#remove geom col
  
  for(i in 1:no_iter){
    set.seed(i*989) #set diff seed number for each model iteration, since this will split your data in different test-train sets / set up how the tree is split each time
    #data <- stor.data #reset data variable for new subset choice
    
    #randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
    rows <- sample(nrow(data))
    data <- data[rows,]
    
    #could add a bootstrapping step here by running code below with a subset of data
    sample = sample.split(data$presence, SplitRatio = split_ratio)
    train = subset(data, sample == TRUE)
    #test  = subset(data, sample == FALSE)
    
    #try stratified sampling
    # sample <- createDataPartition(data$presence, p = split_ratio, list = FALSE)
    # train <- data[ sample,]
    # test  <- data[-sample,]
    
    #remove missing data
    train_complete <- train[complete.cases(train),]; train_complete$presence <- as.numeric(train_complete$presence)
    
    #(1) train model OR (below)
    ####tune hyperparameters - note that this weights presence to equal background (parameter w = TRUE) 
    # snail.ent <- trainMaxNet(
    #   data=subset(train_complete, select=-c(subregion)),
    #   resp='presence',
    #   preds=colnames(subset(train_complete, select=-c(presence,subregion))), #take out response and subregion column from test data ##%%%%remove full zero columns (ag)
    #   #factors='state.numbers',
    #   regMult=c(seq(0.5, 10, by = 0.5)),
    #   classes='lqpht',
    #   verbose=TRUE,
    #   out = c('model')
    # )
    #(2) load trained models
    models_dir <- paste(out_dir,'/models',sep="")
    models <- list.files(path=models_dir, pattern="maxent", all.files=FALSE, full.names=TRUE,recursive=TRUE)
    snail.ent <- readRDS(file=models[i])
    pred_train <- as.data.frame(raster::predict(snail.ent, subset(train_complete, select=-c(presence,subregion)), type='logistic')) #take out response and subregion column from test data
    
    # Store predictions for ensemble model evaluation
    pred_storage_train <- as.data.frame(cbind(rownames(train_complete), as.numeric(pred_train[,1]), train_complete[,1]))
    colnames(pred_storage_train) <- c("row_id", "predictions", "true_values")
    
    if(1 %in% train_complete$presence && 0 %in% train_complete$presence && is.na(pred_train)[1]==F) {
      # Calculate model performance measures for *main model*
      auc <- pROC::roc(response=as.numeric(train_complete[,1]), predictor=as.numeric(pred_train[,1]), auc = TRUE)
      ## Add pROC
      # Calculate prediction vectors needed for pROC/AUC
      pred_storage_prescence_train <- pred_storage_train[which(as.numeric(pred_storage_train$true_values)==1),]
      #minimum <- min(as.numeric(pred_storage$predictions))
      if(length(pred_storage_prescence_train$predictions)>2){
        # Calculate partial roc
        partial_roc <- pROC(continuous_mod=as.numeric(pred_storage_train$predictions),
                            test_data = as.numeric(pred_storage_prescence_train$predictions),
                            n_iter=200,E_percent=100,#-(100*as.numeric(minimum)),
                            boost_percent=50,
                            parallel=FALSE)
        # Store all pauc information
        pauc_mean <- as.numeric(partial_roc$pROC_summary[1])
        pauc_ratio <- as.numeric(partial_roc$pROC_summary[2])
        pauc_pval <- as.numeric(partial_roc$pROC_summary[3])
        pauc_bootstrap_results <- partial_roc$pROC_results
      }else{
        pauc_mean <- NA
        pauc_ratio <- NA
        pauc_pval <- NA
        pauc_bootstrap_results <- NA
      }
      # Calculate best threshold for typical auc calculation to use with metrica package
      best.threshold <- pROC::coords(auc, "best", ret = "threshold") # Find best threshold to use with metrica dataset
      # Use best threshold to build 1/0 prediction dataset (rather than just probabilities)
      metrica.format <- data.frame(cbind(ifelse(train_complete[,1]==1,1,0)),ifelse(as.numeric(pred_train[,1])>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
      #print(metrica.format)
      #print(table(metrica.format))
      if(1 %in% metrica.format$predictions && 0 %in% metrica.format$predictions){
        fscore <- metrica::fscore(data = metrica.format, obs = labels, pred = predictions)$fscore
        sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
        specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
        tss <- sensitivity + specificity - 1
      }else{
        #auc$auc <- NA
        fscore <- NA
        sensitivity <- NA
        specificity <- NA
        tss <- NA
      }
      oob_error <- snail.ent$err.rate[500,1] #first value is # of trees
      #auc$auc <- NA
    }  
    # Store models to build full prediction rasters later on, remember need to build "models" folder in output folder!
    #saveRDS(snail.ent, paste0(out_dir,"/models/maxent_model", i, ".rda"))
    
    #there is probably a nicer way to do this but now I tag each df with an iteration # and bind it to the rest of the iterations
    print(i)
    auc_list <- c(auc_list, auc$auc)
    pauc_mean_list <- c(pauc_mean_list, pauc_mean)
    pauc_ratio_list <- c(pauc_ratio_list, pauc_ratio)
    pauc_pval_list <- c(pauc_pval_list, pauc_pval)
    fscore_list <- c(fscore_list, fscore)
    sensitivity_list <- c(sensitivity_list, sensitivity)
    specificity_list <- c(specificity_list, specificity)
    tss_list <- c(tss_list, tss)
    oob_list <-c(oob_list, oob_error) #store oob error rate
    
    #Start writing / storing output!
    metrics.output.dataframe <- cbind(auc_list, pauc_mean_list, pauc_ratio_list, pauc_pval_list, fscore_list, sensitivity_list, specificity_list, tss_list)
    colnames(metrics.output.dataframe) <- c("auc", "pauc_mean", "pauc_ratio", "pauc_pval","fscore", "sensitivity", "specificity", "tss")
    
    if("subregion" %in% colnames(data)){
      # Create storage dataframe for seven performance measures to fill in for each subregion.
      metrics_subregion_df = data.frame(matrix(vector(), 0, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))),
                                        row.names = NULL, stringsAsFactors=F)
      pred_storage_subregion <- c()
      
      for (k in 1:length(unique(data$subregion))) { #loop through each subregion
        # Build subregion test set (part of larger test set that is in subregion k)
        train_subregion <- train_complete[which(train_complete$subregion==as.character(k)),]
        train_subregion <- subset(train_subregion, select=-c(subregion))
        
        # Create storage dataframes to get overwritten during each iteration of the loop
        subregion_loop_df <- data.frame(matrix(vector(), 1, 10, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'pauc_mean', 'pauc_ratio', 'pauc_pval', 'fscore', 'sensitivity', 'specificity', 'tss'))), 
                                        stringsAsFactors=F, row.names=NULL)
        pred_storage_subregion_loop <- c()
        
        # Build prediction dataframe for subregion k and store if has values, otherwise feed NAs
        if (dim(train_subregion)[1] != 0) {
          pred.subregion <- as.data.frame(raster::predict(snail.ent, subset(train_subregion, select=-c(presence)), type='logistic')) #take out response from test data
          pred.subregion <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,1]), train_subregion[,1]))
          colnames(pred.subregion) <- c("row_id", "predictions", "true_values")
          # Store predictions for ensemble model evaluation
          pred_storage_subregion_loop <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,1]), train_subregion[,1]))
          pred_storage_subregion_loop$subregion <- k
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }else{
          pred_storage_subregion_loop <- as.data.frame(cbind(c(NA),c(NA),c(NA),k))
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }
        
        # Check to make sure subregion has both an occurrence point and background point, otherwise cannot calculate AUC
        if(1 %in% train_subregion$presence && 0 %in% train_subregion$presence) {
          # Calculate metrics
          auc.subregion <- pROC::roc(response=as.numeric(train_subregion[,1]), predictor=as.numeric(pred.subregion[,1]), auc = TRUE) 
          pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
          #print(pred.subregion)
          #print(pred_subregion_prescence_train)
          if(dim(pred_subregion_prescence_train)[1]>2){#length(pred.subregion$predictions)>2 && 
            ## Add pROC
            # Calculate prediction vectors needed for pROC/AUC
            pred_subregion_prescence_train <- pred.subregion[which(as.numeric(pred.subregion$true_values)==1),]
            partial_roc_subregion <- pROC(continuous_mod=as.numeric(pred.subregion$predictions),
                                          test_data = as.numeric(pred_subregion_prescence_train$predictions),
                                          n_iter=200,E_percent=100,
                                          boost_percent=50,
                                          parallel=FALSE)
            # Store all pauc information
            pauc_mean_subregion <- as.numeric(partial_roc_subregion$pROC_summary[1])
            pauc_ratio_subregion <- as.numeric(partial_roc_subregion$pROC_summary[2])
            pauc_pval_subregion <- as.numeric(partial_roc_subregion$pROC_summary[3])
            pauc_bootstrap_results_subregion <- partial_roc_subregion$pROC_results
          }else{
            pauc_mean_subregion <- NA
            pauc_ratio_subregion <- NA
            pauc_pval_subregion <- NA
            pauc_bootstrap_results_subregion <- NA
          }
          best.threshold.subregion <- pROC::coords(auc.subregion, "best", ret = "threshold")
          metrica.format.subregion <- data.frame(cbind(ifelse(train_subregion[,1]==1,1,0)),ifelse(as.numeric(pred.subregion[,1])>=best.threshold.subregion[1,1],1,0)); colnames(metrica.format.subregion) <- c("labels","predictions"); rownames(metrica.format.subregion) <- 1:dim(metrica.format.subregion)[1]
          # Additional check to make sure subregion has both an occurrence point PREDICTED and a background point PREDICTED, otherwise cannot calculate metrica measures
          if(1 %in% metrica.format.subregion$predictions && 0 %in% metrica.format.subregion$predictions){
            fscore.subregion <- metrica::fscore(data = metrica.format.subregion, obs = labels, pred = predictions)$fscore
            sensitivity.subregion <- metrica::recall(data = metrica.format.subregion, obs = labels, pred = predictions)$recall 
            specificity.subregion <- metrica::specificity(data = metrica.format.subregion, obs = labels, pred = predictions)$spec
            tss.subregion <- sensitivity.subregion + specificity.subregion - 1
          }else{ #If second condition not met fill everything with NAs besides AUC
            fscore.subregion <- NA
            sensitivity.subregion <- NA
            specificity.subregion <- NA
            tss.subregion <- NA
          }
          
          #could probably compile into above once working well
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1] #print number of test points predicted on to get a sense of data quantity for each subregion (should be on average ~20% of total amount)
          subregion_loop_df$auc <- auc.subregion$auc
          subregion_loop_df$pauc_mean <- pauc_mean_subregion
          subregion_loop_df$pauc_ratio <- pauc_ratio_subregion
          subregion_loop_df$pauc_pval <- pauc_pval_subregion
          subregion_loop_df$fscore <- fscore.subregion
          subregion_loop_df$sensitivity <- sensitivity.subregion
          subregion_loop_df$specificity <- specificity.subregion
          subregion_loop_df$tss <- tss.subregion
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        } else { #If first condition not met fill everything with NAs (including AUC)
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(train_subregion)[1]
          subregion_loop_df$auc <- NA
          subregion_loop_df$pauc_mean <- NA
          subregion_loop_df$pauc_ratio <- NA
          subregion_loop_df$pauc_pval <- NA
          subregion_loop_df$fscore <- NA
          subregion_loop_df$sensitivity <- NA
          subregion_loop_df$specificity <- NA
          subregion_loop_df$tss <- NA
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        }
      }
      metrics_subregion_df$iter <- i
      metrics_subregion_df_iter <- rbind(metrics_subregion_df_iter, metrics_subregion_df)
      pred_storage_subregion$iter <- i
      pred_storage_subregion_iter <- rbind(pred_storage_subregion_iter, pred_storage_subregion)
    }
    
  }
                                                                                                                                                                                                 
  write.csv(metrics.output.dataframe, paste0(out_dir, "/maxent_baseline_metrics_cropped.csv"))
  write.csv(oob_list, paste0(out_dir, "/maxent_oob_cropped.csv"))
  write.csv(metrics_subregion_df_iter, paste0(out_dir, "/maxent_subregion_metrics_cropped.csv"))
}

# use function!
files <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=TRUE, recursive=FALSE)
output_folder_names <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=FALSE, recursive=FALSE)
files
file_incides <- c(7)#c(1,2,4,5,7:12)
files[file_incides]

for (i in file_incides) {
  print(i)
  #all_data <- readRDS(files[i])
  all_data <- readRDS(files[12]) #change for nationally fit models tested on state data!
  output_folder_name <- paste(tools::file_path_sans_ext(output_folder_names[i]), "_models", sep = "", collapse = NULL)
  random_forest_function(all_data, paste("~/Desktop/r&r_predictions/", output_folder_name, sep=""), 10, 0.8)
  brt_function(all_data, paste("~/Desktop/r&r_predictions/", output_folder_name, sep=""), 10, 0.8)
  maxent_function(all_data, paste("~/Desktop/r&r_predictions/", output_folder_name, sep=""), 10, 0.8)
  print(paste("File", match(i,file_incides), "of", length(files[file_incides]), "done", sep =" "))
}


#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#prediction_df_creation <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/prediction_datasets/national.rds")
#predictions_storage <- data.frame(rep(0,dim(prediction_df_creation)[1]))#; colnames(predictions_storage) <- c("fill")

prediction_function <- function(models_dir, prediction_df_name, model_type, out_dir){
  prediction_df <- readRDS(paste("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/prediction_datasets/",prediction_df_name,sep=""))
  #print(head(prediction_df))
  models_dir_models <- paste(models_dir,"/models",sep="")
  models <- list.files(path=models_dir_models, pattern=model_type, all.files=FALSE, full.names=TRUE,recursive=TRUE)
  #print(models)
  for(i in 1:length(models)){
    if(model_type=="rf"){
      trained_model <- readRDS(file=models[i])
      predictions <- predict(trained_model, newdata=prediction_df[,c(3:dim(prediction_df)[2])], 'prob') #rf
      predictions_col <- predictions[,2]
      #predictions_storage[,paste0("predictions_", i)] <- predictions_col
      prediction_summary_table <- cbind(prediction_df[,1:2], predictions_col)
      colnames(prediction_summary_table) <- c("x", "y", "pred")
      pred_raster <- rasterFromXYZ(prediction_summary_table)
      out_name <- paste(models_dir,'/prediction_rasters/',model_type,"_raster",i,sep="")
      writeRaster(pred_raster, filename=out_name, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
      print(i)
    }
    if(model_type=="maxent"){
      trained_model <- readRDS(file=models[i])
      predictions <- as.data.frame(raster::predict(trained_model, prediction_df[,c(3:dim(prediction_df)[2])], type='logistic')) #maxent
      predictions_col <- predictions[,1]
      #predictions_storage[,paste0("predictions_", i)] <- predictions_col
      prediction_summary_table <- cbind(prediction_df[,1:2], predictions_col)
      colnames(prediction_summary_table) <- c("x", "y", "pred")
      pred_raster <- rasterFromXYZ(prediction_summary_table)
      out_name <- paste(models_dir,'/prediction_rasters/',model_type,"_raster",i,sep="")
      writeRaster(pred_raster, filename=out_name, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
      print(i)
    }
    if(model_type=="brt"){
      trained_model <- readRDS(file=models[i])
      predictions <- as.data.frame(predict(trained_model, newdata=prediction_df[,c(3:dim(prediction_df)[2])], type = 'response')) #brt
      predictions_col <- predictions[,1]
      #predictions_storage[,paste0("predictions_", i)] <- predictions_col
      prediction_summary_table <- cbind(prediction_df[,1:2], predictions_col)
      colnames(prediction_summary_table) <- c("x", "y", "pred")
      pred_raster <- rasterFromXYZ(prediction_summary_table)
      out_name <- paste(models_dir,'/prediction_rasters/',model_type,"_raster",i,sep="")
      writeRaster(pred_raster, filename=out_name, format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
      print(i)
    }
  }
  #calculate mean and sd
  all_tiffs <- list.files(path=paste(models_dir,"/prediction_rasters/",sep=""), pattern=model_type, all.files=FALSE, full.names=TRUE,recursive=TRUE)
  all_tiffs_stack <- raster::stack(all_tiffs)
  mean_raster <- calc(all_tiffs_stack, fun = mean)
  #sd_raster <- calc(all_tiffs_stack, fun = sd)
  writeRaster(mean_raster, filename=paste(models_dir,'/prediction_rasters/',model_types[k],"_mean",sep=""), format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
  #writeRaster(sd_raster, filename=paste(out_dir,model_type,"_",tools::file_path_sans_ext(prediction_df),"_",basename(models_dir),"sd",sep=""), format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
}

#prediction_function(models_dir, prediction_df_name, model_type, out_dir)
folders <- list.dirs("~/Desktop/r&r_predictions")[-1]
folders <- folders[c(1,4,8,11,15,19,22,25,28,31)]
folders <- folders[4]
folders
model_types <- c("rf")
#model_types <- c("brt", "maxent")
for(i in 1:length(folders)){
  for(k in 1:length(model_types)){
    extent <- word(basename(folders[i]), 2, sep = "_")
    prediction_function(folders[i], paste(extent,".rds",sep=""), model_types[k],"~/Desktop/r&r_predictions/")
    print(paste(model_types[k],"done",sep=" "))
  }
  print(paste("Folder",i,"done"))
}

#testing
models_dir <- folders[1]
prediction_df_name <- paste(word(basename(folders[1]), 2, sep = "_"),".rds",sep="")
model_type <- "rf"
out_dir <- "~/Desktop/r&r_figures/"
