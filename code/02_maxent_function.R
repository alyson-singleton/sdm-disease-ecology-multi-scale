# load pROC functions 
source("./code/00_functions.R")

# MODEL FUNCTION USE NOTES:
#     make sure response variable (presence/bkg) is column 1, and lat/long columns have been excluded
#     build two folders for output: (1) main output folder, (2) within that folder build a subfolder named "models"
#     more instructions for function use at the bottom!

maxent_function <- function(data, extrapolation_data, out_dir, no_iter, split_ratio, add_subregion, limited_data_test){
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
  auc_extrapolate_list <- c() 
  fscore_list <- c()
  fscore_extrapolate_list <- c()
  sensitivity_list <- c()
  sensitivity_extrapolate_list <- c()
  specificity_list <- c()
  specificity_extrapolate_list <- c()
  tss_list <- c()
  tss_extrapolate_list <- c()
  oob_list <- c()
  imp_df_iter <- c() 
  pd_df_iter <- c()
  pred_storage_iter <- c()
  pred_storage_subregion_iter <- c()
  stor.data <- data
  
  if(limited_data_test==TRUE){
    #option here to subset to smaller size for data quantity v heterogeneity tests
    presence.rows.sample <- data[sample(nrow(data[which(data$presence==1),]),162),]#currently manual, could update to be an input
    background.rows.sample <- data[sample(nrow(data[which(data$presence==0),]),324),]#^
    data <- rbind(presence.rows.sample, background.rows.sample)
    print(dim(data)) #check subset worked correctly
  }
  
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
  
  #i <- 2
  #data <- all_data
  for(i in final.folds){
    set.seed(i*989) #set diff seed number for each model iteration, since this will split your data in different test-train sets / set up how the tree is split each time
    #data <- stor.data #reset data variable for new subset choice
    
    #randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
    rows <- sample(nrow(data))
    data <- data[rows,]
    
    #could add a bootstrapping step here by running code below with a subset of data
    #sample = sample.split(data$presence, SplitRatio = split_ratio)
    #train = subset(data, sample == TRUE)
    #test  = subset(data, sample == FALSE)
    
    #try stratified sampling
    # sample <- createDataPartition(data$presence, p = split_ratio, list = FALSE)
    # train <- data[ sample,]
    # test  <- data[-sample,]
    
    #with spatial_cv
    train <- data[data$spatialCV_fold != i,]
    train <- subset(train, select=-c(spatialCV_fold))
    #; train <- train[, -14] #take out "spatialCV_fold" column
    test  <- data[data$spatialCV_fold == i,]
    test <- subset(test, select=-c(spatialCV_fold))
    #; test <- test[, -14] #take out "spatialCV_fold" column
  
    #remove missing data
    train_complete <- train[complete.cases(train),]; train_complete$presence <- as.numeric(train_complete$presence)
    test_complete <- test[complete.cases(test),]; test_complete$presence <- as.numeric(test_complete$presence)
    
    #check if running extrapolation tests (i.e. given an extrapolation dataset) and removes missing data if exists (rather than NA)
    if(class(extrapolation_data)=="data.frame"){
      extrapolation_data <- extrapolation_data[complete.cases(extrapolation_data),]; extrapolation_data$presence <- as.numeric(extrapolation_data$presence)
    }
    
    ####tune hyperparameters - note that this weights presence to equal background (parameter w = TRUE) 
    ## ask if you want to include subregion as a categorical variable
    if(add_subregion==TRUE){ #if TRUE/YES, leaves subregion column in as a predictor
      snail.ent <- trainMaxNet(
        data=train_complete,
        resp='presence',
        preds=colnames(train_complete)[-1], #keep subregion categorical variable
        factors='subregion',
        regMult=c(seq(0.5, 5, by = 0.5), 7.5, 10),
        classes='lqph',
        verbose=TRUE,
        out = c('model')
      )
      pred <- as.data.frame(raster::predict(snail.ent, test_complete[,-c(1)], factors='subregion', type='logistic'))  #take out response from test data
    }else{ #otherwise, removes subregion column in as a predictor (column 13)
      snail.ent <- trainMaxNet(
        data=subset(train_complete, select=-c(subregion)),
        resp='presence',
        preds=colnames(subset(train_complete, select=-c(presence,subregion))), #take out response and subregion column from test data ##%%%%remove full zero columns (ag)
        #factors='state.numbers',
        regMult=c(seq(0.5, 10, by = 0.5)),
        classes='lqpht',
        verbose=TRUE,
        out = c('model')
      )
      pred <- as.data.frame(raster::predict(snail.ent, subset(test_complete, select=-c(presence,subregion)), type='logistic')) ##%%%%remove full zero columns (ag)
      pred_train <- as.data.frame(raster::predict(snail.ent, subset(train_complete, select=-c(presence,subregion)), type='logistic')) #take out response and subregion column from test data
      #check if running extrapolation tests (i.e. given an extrapolation dataset), if so as model to predict on extrapolated dataset and 
      #       store seperately for own set of performance measures
      if(class(extrapolation_data)=="data.frame"){
        pred.extrapolate <- as.data.frame(raster::predict(snail.ent, newdata=subset(extrapolation_data, select=-c(presence)), type='logistic')) #take out response and ask model to predict outside of subregion
      }
    }
    
    # Store predictions for ensemble model evaluation
    pred_storage <- as.data.frame(cbind(rownames(test_complete), as.numeric(pred[,1]), test_complete[,1]))
    colnames(pred_storage) <- c("row_id", "predictions", "true_values")
    pred_storage_train <- as.data.frame(cbind(rownames(train_complete), as.numeric(pred_train[,1]), train_complete[,1]))
    colnames(pred_storage_train) <- c("row_id", "predictions", "true_values")
    
    if(1 %in% test_complete$presence && 0 %in% test_complete$presence && is.na(pred)[1]==F) {
      # Calculate model performance measures for *main model*
      auc <- pROC::roc(response=as.numeric(test_complete[,1]), predictor=as.numeric(pred[,1]), auc = TRUE)
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
        print(partial_roc$pROC_summary)
        print(head(partial_roc$pROC_results))
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
      metrica.format <- data.frame(cbind(ifelse(test_complete[,1]==1,1,0)),ifelse(as.numeric(pred[,1])>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
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
    saveRDS(snail.ent, paste0(out_dir,"/models/maxent_model", i, ".rda"))
    
    # Check if running extrapolation tests (i.e. given an extrapolation dataset), and calculate performance metric if so
    if(class(extrapolation_data)=="data.frame"){
      auc.extrapolate <- pROC::roc(response=as.numeric(extrapolation_data[,1]), predictor=as.numeric(pred.extrapolate[,1]), auc = TRUE)
      best.threshold.extrapolate <- pROC::coords(auc.extrapolate, "best", ret = "threshold")
      metrica.format.extrapolate <- data.frame(cbind(ifelse(extrapolation_data[,1]==1,1,0)),ifelse(as.numeric(pred.extrapolate[,1])>=best.threshold.extrapolate[1,1],1,0)); colnames(metrica.format.extrapolate) <- c("labels","predictions"); rownames(metrica.format.extrapolate) <- 1:dim(metrica.format.extrapolate)[1]
      if(1 %in% metrica.format.extrapolate$predictions && 0 %in% metrica.format.extrapolate$predictions){
        fscore.extrapolate <- metrica::fscore(data = metrica.format.extrapolate, obs = labels, pred = predictions)$fscore
        sensitivity.extrapolate <- metrica::recall(data = metrica.format.extrapolate, obs = labels, pred = predictions)$recall 
        specificity.extrapolate <- metrica::specificity(data = metrica.format.extrapolate, obs = labels, pred = predictions)$spec
        tss.extrapolate <- sensitivity.extrapolate + specificity.extrapolate - 1
      }else{
        fscore.extrapolate <- NA
        sensitivity.extrapolate <- NA
        specificity.extrapolate <- NA
        tss.extrapolate <- NA
      }
    }
    
    # Calculate metrics for all subregions if given a subregion column (pretty much always give a subregion column at this point)
    # NOTE: This is a bit of a mess and could probably be cleaned up if you care to take a crack at it! Otherwise it is functional. :)
    # To clarify: This is calculated no matter if subregion is used as a categorical variable, which is the reason we need to include
    #             it in the given datasets even without using subregion as a predictor variable.
    if("subregion" %in% colnames(data)){
      # Create storage dataframe for seven performance measures to fill in for each subregion.
      metrics_subregion_df = data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'fscore', 'sensitivity', 'specificity', 'tss'))),
                                        row.names = NULL, stringsAsFactors=F)
      pred_storage_subregion <- c()
      
      for (k in 1:length(unique(data$subregion))) { #loop through each subregion
        # Build subregion test set (part of larger test set that is in subregion k)
        if(add_subregion==TRUE){ #leave subregion column in if want as a predictor
          test_subregion <- test_complete[which(test_complete$subregion==as.character(k)),]
        }else{ #take subregion column out after subsetting if do *not* want as a predictor 
          test_subregion <- test_complete[which(test_complete$subregion==as.character(k)),]
          test_subregion <- subset(test_subregion, select=-c(subregion)) ##%%%%remove full zero columns (ag)
        }
        
        # Create storage dataframes to get overwritten during each iteration of the loop
        subregion_loop_df <- data.frame(matrix(vector(), 1, 7, dimnames=list(c(), c('subregion', 'total_no', 'auc', 'fscore', 'sensitivity', 'specificity', 'tss'))), 
                                        stringsAsFactors=F, row.names=NULL)
        pred_storage_subregion_loop <- c()
        
        # Build prediction dataframe for subregion k and store if has values, otherwise feed NAs
        if (dim(test_subregion)[1] != 0) {
          pred.subregion <- as.data.frame(raster::predict(snail.ent, subset(test_subregion, select=-c(presence)), type='logistic')) #take out response from test data
          # Store predictions for ensemble model evaluation
          pred_storage_subregion_loop <- as.data.frame(cbind(rownames(pred.subregion), as.numeric(pred.subregion[,1]), test_subregion[,1]))
          pred_storage_subregion_loop$subregion <- k
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }else{
          pred_storage_subregion_loop <- as.data.frame(cbind(c(NA),c(NA),c(NA),k))
          colnames(pred_storage_subregion_loop) <- c("row_id", "predictions", "true_values", "subregion")
        }
        
        # Check to make sure subregion has both an occurrence point and background point, otherwise cannot calculate AUC
        if(1 %in% test_subregion$presence && 0 %in% test_subregion$presence && is.na(pred)[1]==F) {
          # Calculate metrics
          auc.subregion <- pROC::roc(response=as.numeric(test_subregion[,1]), predictor=as.numeric(pred.subregion[,1]), auc = TRUE) 
          best.threshold.subregion <- pROC::coords(auc.subregion, "best", ret = "threshold")
          metrica.format.subregion <- data.frame(cbind(ifelse(test_subregion[,1]==1,1,0)),ifelse(as.numeric(pred.subregion[,1])>=best.threshold.subregion[1,1],1,0)); colnames(metrica.format.subregion) <- c("labels","predictions"); rownames(metrica.format.subregion) <- 1:dim(metrica.format.subregion)[1]
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
          
          #could compile into above once working well
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
    
    #variable importance
    #imp_df <- as.data.frame(ecospat.maxentvarimport(snail.ent, test_complete[,-c(1)], nperm=5)); imp_df$var <- rownames(imp_df)
    pred_fun = function(object, newdata) {
      as.data.frame(raster::predict(object, newdata, type='logistic'))[,1]
    }
    imp_df <- as.data.frame(vi_shap(snail.ent, feature_names = names(subset(train_complete, select=-c(presence))), train = train_complete, pred_wrapper = pred_fun))
    #imp$iteration <- i
    #mxt_shap_importance <- rbind(mxt_shap_importance, imp)
    
    #data to make partial dependence plots
    pd_df = data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),
                       row.names = NULL, stringsAsFactors=F)
    
    for (j in 1:nrow(imp_df)) { #loop through each variable
      
      output <- as.data.frame(pdp::partial(snail.ent, pred.var = imp_df[j,1], prob = TRUE, train = train_complete, type = 'regression', center=FALSE))

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
    oob_list <-c(oob_list, oob_error) #store oob error rate
    
    #extrapolated data metrics
    if(class(extrapolation_data)=="data.frame"){
      auc_extrapolate_list <- c(auc_extrapolate_list, auc.extrapolate$auc)
      fscore_extrapolate_list <- c(fscore_extrapolate_list, fscore.extrapolate)
      sensitivity_extrapolate_list <- c(sensitivity_extrapolate_list, sensitivity.extrapolate)
      specificity_extrapolate_list <- c(specificity_extrapolate_list, specificity.extrapolate)
      tss_extrapolate_list <- c(tss_extrapolate_list, tss.extrapolate)
    }  
    
    imp_df_iter <- rbind(imp_df_iter, imp_df)
    pd_df_iter <- rbind(pd_df_iter, pd_df)
    pred_storage_iter <- rbind(pred_storage_iter, pred_storage)
    
  }
  #Start writing / storing output!
  metrics.output.dataframe <- cbind(auc_list, pauc_mean_list, pauc_ratio_list, pauc_pval_list, fscore_list, sensitivity_list, specificity_list, tss_list); colnames(metrics.output.dataframe) <- c("auc", "pauc_mean", "pauc_ratio", "pauc_pval",
                                                                                                                                                                                                   "fscore", "sensitivity", "specificity", "tss")
  write.csv(metrics.output.dataframe, paste0(out_dir, "/maxent_baseline_metrics.csv"))
  write.csv(oob_list, paste0(out_dir, "/maxent_oob.csv"))
  
  # Extrapolated data output files
  if(class(extrapolation_data)=="data.frame"){
    extrapolate.metrics.output.dataframe <- cbind(auc_extrapolate_list, fscore_extrapolate_list, sensitivity_extrapolate_list, specificity_extrapolate_list, tss_extrapolate_list); colnames(extrapolate.metrics.output.dataframe) <- c("auc", "fscore", "sensitivity", "specificity", "tss")
    write.csv(extrapolate.metrics.output.dataframe, paste0(out_dir, "/maxent_extrapolate_metrics.csv"))
  }
  
  # Subregion data output files
  if("subregion" %in% colnames(data)) {
    write.csv(metrics_subregion_df_iter, paste0(out_dir, "/maxent_subregion_metrics.csv"))
    write.csv(pred_storage_subregion_iter, paste0(out_dir, "/maxent_subregion_preds.csv"))
  }
  
  write.csv(imp_df_iter, paste0(out_dir, "/maxent_imp_df.csv")) 
  write.csv(pd_df_iter, paste0(out_dir,"/maxent_pdp.csv"))
  write.csv(pred_storage_iter, paste0(out_dir,"/maxent_preds.csv")) 
  
}

# use function!
#reminder of function call & inputs
#maxent_function(data, extrapolation_data, out_dir, no_iter, split_ratio, add_subregion, limited_data_test)
all_data <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/glabrata_sudeste.rds")
maxent_function(all_data, NA, "~/Desktop/", 10, 0.8, FALSE, FALSE)

#automate changing directories
files <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=TRUE, recursive=FALSE)
output_folder_names <- list.files(path="~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets", pattern="*.rds", full.names=FALSE, recursive=FALSE)

file_incides <- c(2:13)
#files <- files[file_incides]
files
file_incides

for (i in file_incides) {
  print(i)
  all_data <- readRDS(files[i])
  output_folder_name <- paste(tools::file_path_sans_ext(output_folder_names[i]), "_models", sep = "", collapse = NULL)
  maxent_function(all_data, NA, paste("~/Desktop/r&r_model_runs/", output_folder_name, sep=""), 10, 0.8, FALSE, FALSE)
  print(paste("File", i, "of", length(files), "done", sep =" "))
}

# Function use instructions and examples
# (1) Simple
#     Ex. maxent_function(all_data, NA, "~/Desktop/output_data/straminea_national_models/", 25, 0.8, FALSE, FALSE)
# (2) Add extrapolation data tests
#     To ask the model to also predict and calculate performance measures for an extrapolation data set (or really any
#     data set that is in the same format as all_data) simply give the model a second dataset after "all_data." If
#     you do not want to give the model a second dataset, write NA.
#     Ex. maxent_function(all_data, extrapolation_data, "~/Desktop/output_data/straminea_national_models_w_extrapolation/", 25, 0.8, FALSE, FALSE)
# (3) Add subregion column as predictor or not
#     To ask the model to include subregion column *as a predictor* say TRUE. To exclude, say FALSE. This does not impact
#     whether or not the function will calculate performance measures for each subregion. That happens automatically if a 
#     subregion column can be found in the all_data set.
#     Ex. maxent_function(all_data, NA, "~/Desktop/output_data/straminea_national_models_w_subregion/", 25, 0.8, TRUE, FALSE)
# (4) Limited data test
#     To ask the model to run the "limited data test" say TRUE. To run on full all_data, say FALSE. This is when we want to 
#     test out whether limiting data quantity impacts performance (i.e. data from the entire national area, but with the same 
#     amount of occ/bg/rows as minas gerais). This is not yet coded to run in addition to the normal predictions. If you say true,
#     then it will reduce the occ/bg data at the start of the function and output the results for this new, limited data set in 
#     "baseline metrics." Also note, as of now you yourself need to manually change *how much* the dataset should be limited by
#     inside the function. Check the Brazil Schisto sheets on the "occ/bg dist" tab for quantities to imput.
#     Ex. maxent_function(all_data, NA, "~/Desktop/output_data/straminea_national_models_limited_data/", 25, 0.8, FALSE, TRUE)
