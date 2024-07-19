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


# MODEL FUNCTION USE NOTES:
#     make sure response variable (presence/bkg) is column 1, and lat/long columns have been excluded
#     build two folders for output: (1) main output folder, (2) within that folder build a subfolder named "models"
#     more instructions for function use at the bottom!

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
  
  #_______________________________________________________________________________________________
  #_______________________________________________________________________________________________
  # Prep for spatial cross validation
  analysis_data <- st_as_sf(data,
                            agr = "constant", crs = 4326)
  
  # Spatial Clustering w Kmeans Cross-validation workflow
  set.seed(1)
  clusters <- spatial_clustering_cv(analysis_data, v = 10, cluster_function = "kmeans") #k-means clustering to identify 10 cross-validation folds
  #for loop to create a dataframe that assigns a fold number to each municipality
  splits_df <- c()
  for(i in 1:10){
    new_df <- assessment(clusters$splits[[i]]) #extract municipalities in fold number i
    new_df$spatialCV_fold <- i
    splits_df <- rbind(splits_df, new_df) #bind all municipalities x folds together
  }
  data <- st_drop_geometry(splits_df) #drop shapefiles

  # CHECK: Only do the analysis on fold with presence and bg
  final.folds <- c()
  for (ii in 1:length(unique(data$spatialCV_fold))) {
    if (1 %in% data$presence[which(data$spatialCV_fold==ii)] && 
        0 %in% data$presence[which(data$spatialCV_fold==ii)]){
      final.folds <- c(final.folds, ii)
    }
  }
  print(final.folds)
  #_______________________________________________________________________________________________
  #_______________________________________________________________________________________________
  
  for(i in final.folds){
    set.seed(i*989) #set diff seed number for each model iteration, since this will split your data in different test-train sets / set up how the tree is split each time

    #randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
    rows <- sample(nrow(data))
    data <- data[rows,]
    
    #with spatial_cv
    train <- data[data$spatialCV_fold != i,]
    train <- subset(train, select=-c(spatialCV_fold))
    test  <- data[data$spatialCV_fold == i,]
    test <- subset(test, select=-c(spatialCV_fold))

    #remove missing data
    train_complete <- train[complete.cases(train),]; train_complete$presence <- as.factor(train_complete$presence)
    test_complete <- test[complete.cases(test),]; test_complete$presence <- as.factor(test_complete$presence)
    
    #### Train model
    # sdm_rf <- trainRF(data = subset(train_complete, select=-c(subregion)),
    #                   resp = 'presence',
    #                   preds = 2:ncol(subset(train_complete, select=-c(subregion))))
    # models_dir <- paste(out_dir,'/models',sep="")
    # models <- list.files(path=models_dir, pattern="rf", all.files=FALSE, full.names=TRUE,recursive=TRUE)
    # ii <- match(i,final.folds)
    # sdm_rf <- readRDS(file=models[ii])
    pred <- as.data.frame(predict(sdm_rf, newdata=subset(test_complete, select=-c(presence,subregion)), 'prob')) #take out response and subregion column from test data
    pred_train <- as.data.frame(predict(sdm_rf, newdata=subset(train_complete, select=-c(presence,subregion)), 'prob')) #take out response and subregion column from train data
    
    # Store predictions for ensemble model evaluation
    test_complete[,1] <- as.numeric(as.character(test_complete[,1]))
    train_complete[,1] <- as.numeric(as.character(train_complete[,1]))
    pred_storage <- as.data.frame(cbind(rownames(test_complete), as.numeric(pred[,2]), test_complete[,1]))
    colnames(pred_storage) <- c("row_id", "predictions", "true_values")
    pred_storage_train <- as.data.frame(cbind(rownames(train_complete), as.numeric(pred_train[,2]), train_complete[,1]))
    colnames(pred_storage_train) <- c("row_id", "predictions", "true_values")

    # Check fold has 
    if(1 %in% test_complete$presence && 0 %in% test_complete$presence) {
      # Calculate model performance measures for *main model*
      auc <- pROC::roc(response=as.numeric(test_complete[,1]), predictor=as.numeric(pred[,2]), auc = TRUE)
      print(auc$auc)
      #pROC::plot.roc(roc(response=as.numeric(test_complete[,1]), predictor=as.numeric(pred[,2]), levels=c(1,2)))
      ## Add pROC
      # Calculate prediction vectors needed for pROC/AUC
      pred_storage_prescence <- pred_storage[which(as.numeric(pred_storage$true_values)==1),]
      pred_storage_prescence_train <- pred_storage_train[which(as.numeric(pred_storage_train$true_values)==1),]
      #minimum <- min(as.numeric(pred_storage$predictions))
      if(length(pred_storage_prescence$predictions)>2){
        # Calculate partial roc
        partial_roc <- pROC(continuous_mod=as.numeric(pred_storage_train$predictions),
                            test_data = as.numeric(pred_storage_prescence$predictions),
                            n_iter=200,E_percent=100,#-(100*as.numeric(minimum)),
                            boost_percent=50,
                            parallel=FALSE)
        # Store all pauc information
        pauc_mean <- as.numeric(partial_roc$pROC_summary[1])
        pauc_ratio <- as.numeric(partial_roc$pROC_summary[2])
        pauc_pval <- as.numeric(partial_roc$pROC_summary[3])
        print(pauc_mean)
        pauc_bootstrap_results <- partial_roc$pROC_results
      }else{
        pauc_mean <- NA
        pauc_ratio <- NA
        pauc_pval <- NA
        pauc_bootstrap_results <- NA
      }
      # Calculate best threshold for typical auc calculation to use with metrica package
      best.threshold <- pROC::coords(auc, "best", ret = "threshold")
      # Use best threshold to build 1/0 prediction dataset (instead of continuous)
      metrica.format <- data.frame(cbind(ifelse(test_complete[,1]==1,1,0)),ifelse(as.numeric(pred[,2])>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
      if(1 %in% metrica.format$predictions && 0 %in% metrica.format$predictions){
        fscore <- metrica::fscore(data = metrica.format, obs = labels, pred = predictions)$fscore
        sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
        specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
        tss <- sensitivity + specificity - 1
        oob_error <- sdm_rf$err.rate[500,1] #first value is # of trees
      }else{
        #auc$auc <- NA
        #pauc_mean <- NA
        #pauc_ratio <- NA
        #pauc_pval <- NA
        #pauc_bootstrap_results <- NA
        fscore <- NA
        sensitivity <- NA
        specificity <- NA
        tss <- NA
        oob_error <- NA
      }
    }  
    # Store models to build full prediction rasters later on, remember need to build "models" folder in output folder!
    saveRDS(sdm_rf, paste0(out_dir,"/models/rf_model", i, ".rda"))

    if("subregion" %in% colnames(data)){
      # Create storage dataframe for seven performance measures to fill in for each subregion.
      metrics_subregion_df = data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'fscore', 'sensitivity', 'specificity', 'tss'))),
                                    row.names = NULL, stringsAsFactors=F)
      pred_storage_subregion <- c()
      
      for (k in 1:length(unique(data$subregion))) { #loop through each subregion
        test_subregion <- test_complete[which(test_complete$subregion==as.character(k)),]
        test_subregion <- subset(test_subregion, select=-c(subregion))
        
        # Create storage dataframes to get overwritten during each iteration of the loop
        subregion_loop_df <- data.frame(matrix(vector(), 1, 7, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'fscore', 'sensitivity', 'specificity', 'tss'))), 
                                        stringsAsFactors=F, row.names=NULL)
        pred_storage_subregion_loop <- c()
        
        # Build prediction dataframe for subregion k and store if has values, otherwise feed NAs
        if (dim(test_subregion)[1] != 0) {
          pred.subregion <- as.data.frame(predict(sdm_rf, newdata=subset(test_subregion, select=-c(presence)), 'prob')) #take out response from test data
          # Store predictions for ensemble model evaluation
          pred_storage_subregion_loop <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,2]), test_subregion[,1]))
          pred_storage_subregion_loop$subregion <- k
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }else{
          pred_storage_subregion_loop <- as.data.frame(cbind(c(NA),c(NA),c(NA),k))
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }
        
        # Check to make sure subregion has both an occurrence point and background point, otherwise cannot calculate AUC
        if(1 %in% test_subregion$presence && 0 %in% test_subregion$presence) {
          # Calculate metrics
          auc.subregion <- pROC::roc(response=as.numeric(test_subregion[,1]), predictor=as.numeric(pred.subregion[,2]), auc = TRUE) 
          best.threshold.subregion <- pROC::coords(auc.subregion, "best", ret = "threshold")
          metrica.format.subregion <- data.frame(cbind(ifelse(test_subregion[,1]==1,1,0)),ifelse(as.numeric(pred.subregion[,2])>=best.threshold.subregion[1,1],1,0)); colnames(metrica.format.subregion) <- c("labels","predictions"); rownames(metrica.format.subregion) <- 1:dim(metrica.format.subregion)[1]
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
          subregion_loop_df$total_no <- dim(test_subregion)[1] #print number of test points predicted on to get a sense of data quantity for each subregion (should be on average ~20% of total amount)
          subregion_loop_df$auc <- auc.subregion$auc
          subregion_loop_df$fscore <- fscore.subregion
          subregion_loop_df$sensitivity <- sensitivity.subregion
          subregion_loop_df$specificity <- specificity.subregion
          subregion_loop_df$tss <- tss.subregion
          #Bind to larger data set to print
          metrics_subregion_df <- rbind(metrics_subregion_df, subregion_loop_df)
          pred_storage_subregion <- rbind(pred_storage_subregion, pred_storage_subregion_loop)
        } else { #If first condition not met fill everything with NAs (including AUC)
          subregion_loop_df$subregion <- k
          subregion_loop_df$total_no <- dim(test_subregion)[1]
          subregion_loop_df$auc <- NA
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
    
    pred_fun = function(object, newdata) {
      as.data.frame(predict(object, newdata, 'prob'))[,2]
    }
  
    imp_df <- as.data.frame(vi_shap(sdm_rf, feature_names = names(subset(train_complete, select=-c(presence,subregion))), train = train_complete, pred_wrapper = pred_fun))

    #data to make partial dependence plots
    pd_df = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),
                       row.names = NULL, stringsAsFactors=F)
    
    for (j in 1:nrow(imp_df)) { #loop through each variable
      
      output <- as.data.frame(pdp::partial(sdm_rf, pred.var = imp_df[j,1], prob = TRUE, train = train_complete))
      
      loop_df <- data.frame(matrix(vector(), nrow(output), 3,
                                   dimnames=list(c(), c('variable', 'value','yhat'))), stringsAsFactors=F,
                            row.names=NULL)
      
      loop_df$value <- output[[1]]
      loop_df$yhat <- output[[2]]
      loop_df$variable <- imp_df[j,1]
      
      pd_df <- rbind(pd_df, loop_df)
      
    }
    
    #there is probably a nicer way to do this but now I tag each df with an iteration # and bind it to the rest of the iterations
    print(i)
    imp_df$iter <- i; pd_df$iter <- i; pred_storage$iter <- i
    auc_list <- c(auc_list, auc$auc)
    pauc_mean_list <- c(pauc_mean_list, pauc_mean)
    pauc_ratio_list <- c(pauc_ratio_list, pauc_ratio)
    pauc_pval_list <- c(pauc_pval_list, pauc_pval)
    fscore_list <- c(fscore_list, fscore)
    sensitivity_list <- c(sensitivity_list, sensitivity)
    specificity_list <- c(specificity_list, specificity)
    tss_list <- c(tss_list, tss)
    oob_list <-c(oob_list, oob_error)
  
    imp_df_iter <- rbind(imp_df_iter, imp_df)
    pd_df_iter <- rbind(pd_df_iter, pd_df)
    pred_storage_iter <- rbind(pred_storage_iter, pred_storage)
  }
  
  # Store metrics
  metrics.output.dataframe <- cbind(auc_list, pauc_mean_list, pauc_ratio_list, pauc_pval_list, fscore_list, sensitivity_list, specificity_list, tss_list); colnames(metrics.output.dataframe) <- c("auc", "pauc_mean", "pauc_ratio", "pauc_pval",
                                                                                                                                                                                                   "fscore", "sensitivity", "specificity", "tss")
  # Write output files
  write.csv(metrics.output.dataframe, paste0(out_dir, "/rf_baseline_metrics.csv"))
  write.csv(oob_list, paste0(out_dir, "/rf_oob.csv"))
  write.csv(metrics_subregion_df_iter, paste0(out_dir, "/rf_subregion_metrics.csv"))
  write.csv(pred_storage_subregion_iter, paste0(out_dir, "/rf_subregion_preds.csv"))
  write.csv(imp_df_iter, paste0(out_dir, "/rf_imp_df.csv")) 
  write.csv(pd_df_iter, paste0(out_dir,"/rf_pdp.csv")) 
  write.csv(pred_storage_iter, paste0(out_dir,"/rf_preds.csv")) 
}

# use function!
#reminder of function call & inputs
#random_forest_function(data, out_dir, no_iter, split_ratio)
all_data <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sudeste.rds")
random_forest_function(all_data, "~/Desktop/r&r_model_runs/tenagophila_sudeste_models/", 10, 0.8)

#automate changing directories
files <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=TRUE, recursive=FALSE)
output_folder_names <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=FALSE, recursive=FALSE)
files
file_incides <- c(4:13)
files[file_incides]
file_incides

for (i in file_incides) {
  print(i)
  all_data <- readRDS(files[i])
  output_folder_name <- paste(tools::file_path_sans_ext(output_folder_names[i]), "_models", sep = "", collapse = NULL)
  random_forest_function(all_data, paste("~/Desktop/r&r_model_runs/", output_folder_name, sep=""), 10, 0.8)
  print(paste("File", i, "of", length(files[file_incides]), "done", sep =" "))
}

# Function use instructions and examples
# (1) Simple
#     Ex. random_forest_function(all_data, NA, "~/Desktop/output_data/straminea_national_models/", 25, 0.8, FALSE, FALSE)
# (2) Add extrapolation data tests
#     To ask the model to also predict and calculate performance measures for an extrapolation data set (or really any
#     data set that is in the same format as all_data) simply give the model a second dataset after "all_data." If
#     you do not want to give the model a second dataset, write NA.
#     Ex. random_forest_function(all_data, extrapolation_data, "~/Desktop/output_data/straminea_national_models_w_extrapolation/", 25, 0.8, FALSE, FALSE)
# (3) Add subregion column as predictor or not
#     To ask the model to include subregion column *as a predictor* say TRUE. To exclude, say FALSE. This does not impact
#     whether or not the function will calculate performance measures for each subregion. That happens automatically if a 
#     subregion column can be found in the provided dataset.
#     Ex. random_forest_function(all_data, NA, "~/Desktop/output_data/straminea_national_models_w_subregion/", 25, 0.8, TRUE, FALSE)
# (4) Limited data test
#     To ask the model to run the "limited data test" say TRUE. To run on full all_data, say FALSE. This is when we want to 
#     test out whether limiting data quantity impacts performance (i.e. data from the entire national area, but with the same 
#     amount of occ/bg/rows as minas gerais). This is not yet coded to run in addition to the normal predictions. If you say true,
#     then it will reduce the occ/bg data at the start of the function and output the results for this new, limited data set in 
#     "baseline metrics." Also note, as of now you yourself need to manually change *how much* the dataset should be limited by
#     inside the function. Check the Brazil Schisto sheets on the "occ/bg dist" tab for quantities to input.
#     Ex. random_forest_function(all_data, NA, "~/Desktop/output_data/straminea_national_models_limited_data/", 25, 0.8, FALSE, TRUE)

#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________
#_______________________________________________________________________________________________________________________

rf_variable_importance <- function(data){
 library(ggplot2)
  library(tidyr)
  library(dplyr)
  
   sum_import <- data %>%
    group_by(var) %>%
    summarise(mean.importance = mean(MeanDecreaseGini, na.rm = TRUE),
              sd.importance = sd(MeanDecreaseGini, na.rm = TRUE),
              n.importance = n()) %>%
    mutate(se.importance = sd.importance / sqrt(n.importance),
           lower.ci.importance = mean.importance - qt(1 - (0.05 / 2), n.importance - 1) * se.importance,
           upper.ci.importance = mean.importance + qt(1 - (0.05 / 2), n.importance - 1) * se.importance)
  
  importance_plot <- ggplot(sum_import, aes(x = reorder(var, mean.importance), y = mean.importance)) +
    geom_point(size = 3) + xlab('feature') + ylab('importance') +
    geom_errorbar(aes(ymin = lower.ci.importance, ymax = upper.ci.importance), position = "dodge", width = 0.4, size = 1.5) +
    #scale_color_manual(values=c('grey',colors[c(3,5,8)])) + 
    coord_flip() +
    theme_bw(base_size = 14)
  
  print(importance_plot)
}
