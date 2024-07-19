## Build predictions maps and store for figure building

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
# Run baby models on all available training data for best prediction maps
all_data <- readRDS("~/Desktop/late_feb_runs/datasets/glabrata_mg.rds")
# prepare data
data <- all_data
#data <- dataset
#randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
rows <- sample(nrow(data))
data <- data[rows,]
train <- subset(data, select=-c(geometry))#remove geom col

# RF
train_complete <- train[complete.cases(train),]; train_complete$presence <- as.factor(train_complete$presence)
sdm_rf <- trainRF(data = subset(train_complete, select=-c(subregion)), 
                  resp = 'presence', 
                  preds = 2:ncol(subset(train_complete, select=-c(subregion))))

# MaxEnt
train_complete <- train[complete.cases(train),]; train_complete$presence <- as.numeric(train_complete$presence)
snail.ent <- trainMaxNet(
  data=train_complete,
  resp='presence',
  preds=colnames(subset(train_complete, select=-c(presence,subregion))), #take out response and subregion column from test data ##%%%%remove full zero columns (ag)
  #factors='state.numbers',
  regMult=c(seq(0.5, 10, by = 0.5)),
  classes='lqpht',
  verbose=TRUE,
  out = c('model')
)

# BRT
train_complete <- train[complete.cases(train),]#; train_complete$presence <- as.factor(train_complete$presence)
sdm_brt <- trainBRT(data = subset(train_complete, select=-c(subregion)), 
                    resp = 'presence', 
                    preds = names(train_complete)[2:ncol(subset(train_complete, select=-c(subregion)))],
                    w = FALSE,
                    out = c("model")#,
                    #minTrees = 800,
                    #maxTrees = 14000,
                    #verbose = TRUE,
                    #treeComplexity = c(1:16)
)

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
## Load environmental data for all of Brazil (or correct subregion)

#a file containing all the rasters of covariates used in the model.
env_data <- list.files(path="~/Desktop/brazil_schisto_raster_folders/national_reduced_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
env_data <- list.files(path="~/Desktop/brazil_schisto_raster_folders/sp_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
env_data <- list.files(path="~/Desktop/brazil_schisto_raster_folders/mg_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
e <- raster::stack(env_data)
# sac1 <- blockCV::cv_spatial_autocor(r = e, 
#                            num_sample = 5000)
#create a data-frame that has the covariate for each grid cell in Brazil (or your area of interest)
prediction_df <- as.data.frame(rasterToPoints(e)) ##this gets raster value for every grid cell of Brazil
# colnames(prediction_df) <- c('x', 'y','agriculture','bio13_prec','bio2_diurn', 'bio4_temp_', 'bio5_temp_', 'forest', 'hnd',
#                              'pasture', 'pH', 'soilWater','urban_mining')
# colnames(prediction_df) <- c('x', 'y','agriculture', 'bio11_temp', 'bio16_prec', 'bio17_prec', 'bio4_temp_', 'forest', 'hnd',
#                              'pasture', 'pH', 'soilWater','urban_mining')
# colnames(prediction_df) <- c('x', 'y','agriculture','bio11_temp', 'bio17_prec', 'bio2_diurn', 'bio4_temp_', 'clay', 'forest', 'hnd',
#                              'pasture', 'pH', 'soilWater', 'total_population', 'urban_mining')
# colnames(prediction_df) <- c('x', 'y','agriculture','bio11_temp', 'bio17_prec', 'bio2_diurn', 'bio4_temp_', 'clay', 'forest', 'hnd',
#                              'pasture', 'pH', 'soilWater', 'total_population', 'urban_mining')
# colnames(prediction_df) <- c('x', 'y', 'ag_mosiac', 'bio11_temp', 'bio17_prec', 'bio2_diurn', 'bio4_temp_', 'clay',
#                              'distance_high_pop', 'hnd', 'pH', 'soilWater')
colnames(prediction_df) <- c('x', 'y', 'ag_mosiac', 'bio11_temp', 'bio16_prec', 'bio17_prec', 'bio4_temp_', 'clay',
                             'distance_high_pop', 'hnd', 'pH', 'soilWater')
## Missingness functions
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

is.na.data.frame <- function(x){
  do.call(cbind, lapply(x, is.na))
}

#Remove NA and NaN values from prediction_df
prediction_df[is.nan(prediction_df)] <- 0
prediction_df[is.na(prediction_df)] <- 0

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
# In sample AUC, sensitivity, and specificity for cropping analysis

model <- readRDS(file="~/Desktop/subregion_tests/tenagophila_national_models/models/brt_model1.rda")
dataset <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sp.rds")
dataset <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_national.rds")

#randomizing data b/c otherwise data goes from all 1s to all 0s, helps a little with model performance
rows <- sample(nrow(dataset))
dataset <- dataset[rows,]
train <- subset(dataset, select=-c(geometry))#remove geom col

#RANDOM FOREST
train_set <- train[complete.cases(train),]
train_set$presence <- as.factor(train_set$presence) #RF
pred <- as.data.frame(predict(model, newdata=subset(train_set, select=-c(presence,subregion)), 'prob'))[,2] #RF
auc <- pROC::roc(response=as.numeric(train_set[,1]), predictor=as.numeric(pred), levels=c(1,2), auc = TRUE, plot=T) #RF
auc$auc

#MAXENT
train_set <- train[complete.cases(train),]; train_set$presence <- as.numeric(train_set$presence) #MaxEnt
pred <- as.data.frame(raster::predict(model, subset(train_set, select=-c(presence,subregion)), type='logistic')) #MaxEnt
auc <- pROC::roc(response=as.numeric(train_set[,1]), predictor=as.numeric(pred[,1]), auc = TRUE, plot=T) #MaxEnt
auc$auc

#BRT
train_set <- train[complete.cases(train),]#; train_complete$presence <- as.factor(train_complete$presence)
pred <- as.data.frame(predictEnmSdm(model, newdata=subset(train_set, select=-c(presence,subregion))))[,1]#, type = 'response')) #BRT
auc <- pROC::roc(response=as.numeric(train_set[,1]), predictor=as.numeric(pred), auc = TRUE, plot=T) #BRT
auc$auc

#ALL
best.threshold <- pROC::coords(auc, "best", ret = "threshold") # Find best threshold to use with metrica dataset
metrica.format <- data.frame(cbind(ifelse(train_set[,1]==1,1,0)),ifelse(as.numeric(pred)>=best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
sensitivity
specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
specificity
tss <- sensitivity + specificity - 1
tss
# add to storage dataframe
col <- c(auc$auc,sensitivity,specificity,"State", "Tenagophila", "BRT")
#newfigvals <- c()
newfigvals <- rbind(newfigvals, col)
newfigvals

newfigvals <- as.data.frame(newfigvals)
colnames(newfigvals) <- c("auc", "sensitivity", "specificity", "model_geo_ex", "species", "model_type")
rownames(newfigvals) <- c(1:12)
newfigvals[,1] <- as.numeric(newfigvals[,1]); newfigvals[,2] <- as.numeric(newfigvals[,2]); newfigvals[,3] <- as.numeric(newfigvals[,3])
newfigvals[,c(1:3)] <- round(newfigvals[,c(1:3)], 4)
newfigvals
saveRDS(newfigvals, "~/Desktop/late_feb_runs/datasets/cropping_performance_measures.rds")
test<-readRDS("~/Desktop/late_feb_runs/datasets/cropping_performance_measures.rds")
#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
# Loop for predictions over iterations

#initialize storage objects 
predictions_sum <- data.frame(rep(0,dim(prediction_df)[1]))
predictions_storage <- data.frame(rep(0,dim(prediction_df)[1])); colnames(predictions_storage) <- c("fill")

#Grab set of models for averaging (subregion x model type)
models <- list.files(path="~/Desktop/week_four_runs/tenagophila_national_models/models", pattern="brt", all.files=FALSE, full.names=TRUE,recursive=TRUE)

#Evaluate and store predictions
for(i in 1:length(models)){
  trained_model <- readRDS(file=models[i]) 
  predictions <- predict(sdm_rf, newdata=prediction_df[,c(3:dim(prediction_df)[2])], 'prob') #rf
  predictions <- as.data.frame(raster::predict(snail.ent, prediction_df[,c(3:dim(prediction_df)[2])], type='logistic')) #maxent
  predictions <- as.data.frame(predict(sdm_brt, newdata=prediction_df[,c(3:dim(prediction_df)[2])], type = 'response')) #brt
  predictions_sum <- predictions_sum + predictions[,1] #2 for rf, 1 for maxent, 1 for brt
  predictions_storage[,paste0("predictions_", i)] <- predictions[,1]
}

#predictions_test <- as.data.frame(raster::predict(trained_model, prediction_df[,c(3:13)])) #maxent

#remove fill column
predictions_storage <- predictions_storage[,2:dim(predictions_storage)[2]]

#calculate mean and sd, could add more if interested in future
prediction_summary_table <- cbind(prediction_df[,1:2], rowMeans(predictions_storage), apply(predictions_storage, 1, sd))
colnames(prediction_summary_table) <- c("x", "y", "mean", "sd")
prediction_summary_storage <- cbind(prediction_summary_table, predictions_storage)

#subset data from to only longitude, latitude, and probability
mean_tiff_df <- prediction_summary_table[, c("x", "y", "mean")]
sd_tiff_df <- prediction_summary_table[, c("x", "y", "sd")]

#convert probabilities to rasters
mean_raster <- rasterFromXYZ(mean_tiff_df)
sd_raster <- rasterFromXYZ(sd_tiff_df)

# just predictions no mean no uncertainty
#2 for rf, 1 for maxent, 1 for brt
prediction_summary_table <- cbind(prediction_df[,1:2], predictions[,1]); colnames(prediction_summary_table) <- c("x","y","predictions")
pred_raster <- rasterFromXYZ(prediction_summary_table)

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
#Save rasters / predictions / models

writeRaster(mean_raster, filename='~/Desktop/week_four_runs/prediction_maps/brt_tenagophila_mean_nov8.tif', format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
writeRaster(sd_raster, filename='~/Desktop/week_four_runs/prediction_maps/brt_tenagophila_sd_nov8.tif', format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
writeRaster(pred_raster, filename='~/Desktop/late_feb_runs/prediction_maps/brt_glabrata_mg_fulltrain.tif', format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
out_dir <- "~/Desktop/late_feb_runs/"
saveRDS(sdm_brt, paste0(out_dir,"model_storage/brt_glabrata_mg_fulltrain.rda"))

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
## Quick check / test figs

zlim <- range(c(0,1))
str_name <-'~/Desktop/late_feb_runs/prediction_maps/maxent_glabrata_national_fulltrain.tif' 
pred_raster <- brick(str_name)
par(mar = c(1, 1, 1, 1))
plot(pred_raster)

#set standard scale for this code block
breakpoints2 <- c(seq(0, 0.9730816, (0.9730816-0)/8))
colors2 <- c(RColorBrewer::brewer.pal(9, "BuPu"))[2:9]
#zlim <- range(c(0.104,upperbound))

plot(pred_raster,  
     main="",
     axes = FALSE, box = FALSE,
     breaks = breakpoints2, col = colors2,
     cex.main = 1.7, cex.lab=1.7#, zlim=zlim
)

#___________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________
# crop national rasters to state shapefiles

### mg shapefile
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),] #mesoregion as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),])
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)
mg_shapefile <- geo_extent
### sp shapefile
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),] #mesoregion as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),])
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)
sp_shapefile <- geo_extent
### national tenagopihla rasters
load.raster <- raster("~/Desktop/r&r_predictions/tenagophila_national_models/prediction_rasters/maxent_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/tenagophila_national_models/prediction_rasters/rf_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/tenagophila_national_models/prediction_rasters/brt_mean.tif")

### national straminea rasters
load.raster <- raster("~/Desktop/r&r_predictions/straminea_national_models/prediction_rasters/maxent_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/straminea_national_models/prediction_rasters/rf_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/straminea_national_models/prediction_rasters/brt_mean.tif")

### national glabrata rasters
load.raster <- raster("~/Desktop/r&r_predictions/glabrata_national_models/prediction_rasters/maxent_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/glabrata_national_models/prediction_rasters/rf_mean.tif")
load.raster <- raster("~/Desktop/r&r_predictions/glabrata_national_models/prediction_rasters/brt_mean.tif")

plot(load.raster)
masked <- mask(x = load.raster, mask = geo_extent)
raster.crop <- crop(masked, geo_extent)

colors2 <- c(RColorBrewer::brewer.pal(9, "BuPu"))[2:9]
plot(raster.crop,
     #main="Boosted Regression Tree",
     axes = FALSE, box = FALSE, legend = FALSE,
     breaks = c(seq(raster.crop@data@min, raster.crop@data@max, (raster.crop@data@max-raster.crop@data@min)/8)), col = colors2,
     #xlab="Longitude", ylab="Latitude",
     cex.main = 2, cex.lab=2, zlim=range(c(raster.crop@data@min,raster.crop@data@max)))
#plot(geo_extent, add = TRUE)
writeRaster(raster.crop, filename='~/Desktop/r&r_predictions/glabrata_mg_models/prediction_rasters/brt_national_cropped.tif', format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)
