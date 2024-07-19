#------------------------------------------------------
#------------------------------------------------------
### Prepare rds data for input into ML functions. Could be changed into functions later on. 
#       Might need to be project specific.
#------------------------------------------------------
#------------------------------------------------------

#load arial if being weird
library(showtext)
font_add("Arial", "/Library/Fonts/Arial.ttf")  # Use the actual file path
showtext_auto()

# install packages, load libraries
pacman::p_load(tidyverse, dplyr, kableExtra, knitr, sf, raster, terra, spData, 
               spDataLarge, tmap, leaflet, ggplot2, spThin, mapdata, rnaturalearth, 
               rnaturalearthdata, devtools, brazilmaps, mapview, grid, 
               gridExtra, RColorBrewer, randomForest, caret, caTools, pdp,
               tidyr, dismo, dismotools, sp, enmSdm, grpregOverlap, stats,readr)

#------------------------------------------------------
# CHOICE 1: Which occurrence data to use? Either run Expert section *OR* GBIF section.
#           This is also where you could easily create a combo Expert / GBIF dataset.
#------------------------------------------------------

# EXPERT SNAIL DATA
# Load RDS file
#environmental_covariate_df <- readRDS(file = "~/Desktop/doctorate/ch1 brazil schisto/multi-scale-sdm-schisto/final_data/clean_data/bra_fio_all_snails_all_covariates.rds")
environmental_covariate_df <- readRDS(file = "~/Desktop/doctorate/ch1 brazil schisto/multi-scale-sdm-schisto/final_data/clean_data/all_expert_all_gbif_data_dec12.rds")

# Quick thinning test
gbif.test.glab <- read.table("~/Downloads/0244816-230224095556074.csv", sep="\t", header=TRUE)
gbif.test.stram <- read.table("~/Downloads/0244824-230224095556074.csv", sep="\t", header=TRUE)

gbif.test <- gbif.test.stram
gbif.test <- as.data.frame(gbif.test)
gbif.test[,22] <- as.numeric(gbif.test[,22])
gbif.test[,23] <- as.numeric(gbif.test[,23])
s <- gridSample(gbif.test[c(23,22)], rast, n=1)
s.thin <- gbif.test[row.names(s),]
thin.occ <- occ.target[row.names(s),]; thin.occ <- thin.occ[complete.cases(thin.occ), ]


# Only keep variables of interest
# caroline's fav list: 'bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec'
# best convergence for small models: 'bio2_diurn', 'bio4_temp_', 'bio5_temp_', 'bio13_prec'
# mixture: 'bio2_diurn', 'bio4_temp_', 'bio11_temp', 'bio17_prec'

variable_list <- c('bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                   'ag_mosiac', 'distance_high_pop',
                   'hnd', 'soilWater', 'pH', 'clay',
                   'species', 'year', 'latitude', 'longitude', 'source')
environmental_covariate_df <- environmental_covariate_df[,variable_list]

# Catch all to remove any incomplete rows
environmental_covariate_df <- environmental_covariate_df[complete.cases(environmental_covariate_df), ]
environmental_covariate_df_stor <- environmental_covariate_df
saveRDS(environmental_covariate_df, "~/Desktop/full_dataset.rds")
#------------------------------------------------------

# GBIF v Expert SNAIL DATA

# # Pick GBIF or expert for occurrence data counting tests
environmental_covariate_df <- environmental_covariate_df %>%
  filter(source == "expert") #expert, gbif

#------------------------------------------------------
# CHOICE 2: What is our geographic extent of interest? Only run the six lines that correspond to
#           the region you want the model to run on. Store as "geo_extent" which is what will be
#           referenced in the rest of the code. Also store appropriate subregion geometry set to 
#           use later. Brazil shapefiles pulled from Brazilian govern website. Perhaps worth 
#           putting on GitHub.
#------------------------------------------------------

# Set multiplier
multiplier <- 2

### National
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Pais_2021/BR_Pais_2021.shp")
geo_extent_shp <- st_union(geo_extent_shp)
subregions <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_UF_2021/BR_UF_2021.shp") #state as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)

### Regional (Sudeste)
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_UF_2021/BR_UF_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$NM_REGIAO=="Sudeste"),] #state as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$NM_REGIAO=="Sudeste"),])
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)

### State (Minas Gerais)
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),] #mesoregion as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),])
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)

### State (Sao Paulo)
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),] #mesoregion as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),])
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)

### State (Pernambuco)
geo_extent_shp <- read_sf("~/Desktop/doctorate/mordecai_lab/brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="PE"),] #mesoregion as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="PE"),])
geo_extent <- st_as_sf(geo_extent_shp)
geo_extent$x <- st_transform(geo_extent$x, 4326)

#------------------------------------------------------
# CHOICE 3: Choose target species of interest.
#           Expert options: 'glabrata', 'tenagophila', 'straminea'
#           GBIF options: 'Biomphalaria glabrata', 'Biomphalaria tenagophila', 'Biomphalaria straminea'
#------------------------------------------------------

# Build target species occurrence rows
# Set target species
#target_species_vector <- c('glabrata', 'Biomphalaria glabrata', 'B. glabrata')
target_species_vector <- c('tenagophila', 'Biomphalaria tenagophila', 'B. tenagophila')
#target_species_vector <- c('straminea', 'Biomphalaria straminea', 'B. straminea')

# Grab target species rows
occ.target <- environmental_covariate_df %>%
  filter(species %in% target_species_vector)

#Expert v GBIF
# occ.target <- occ.target %>%
#   filter(source == "expert") #expert, gbif

# Clean lat long columns for thinning procedure
occ_target_lat_lon <- as.data.frame(occ.target[,c("species", "longitude", "latitude")])
occ_target_lat_lon$latitude <- as.numeric(occ_target_lat_lon$latitude)
occ_target_lat_lon$longitude <- as.numeric(occ_target_lat_lon$longitude)

# Load raster to use for grid cell references
rast <- raster("~/Desktop/doctorate/ch1 brazil schisto/schisto_archives/regional_rasters/national_mean_EVI.tif")
# One point per grid cell
s <- gridSample(occ_target_lat_lon[2:3], rast, n=1)
thin.occ <- occ.target[row.names(s),]; thin.occ <- thin.occ[complete.cases(thin.occ), ]

# Subset to region of interest (national, regional, state, etc.)
occ_sf_points <- st_as_sf(x = thin.occ, 
                          coords = c("longitude", "latitude"), crs = 4326)
occ_sf_points$geometry <- st_set_crs(occ_sf_points$geometry, 4326)
occ_by_region <- st_intersects(occ_sf_points$geometry, geo_extent$x, sparse = FALSE)
occ_sf_points_region <- occ_sf_points[occ_by_region[,1],]
thin.occ <- occ_sf_points_region

# Add subregion column
thin.occ.subregion <- st_intersects(thin.occ$geometry, subregions$geometry, sparse = FALSE)
thin.occ$subregion <- apply(thin.occ.subregion, 1, which)
#thin.occ <- st_set_geometry(thin.occ, NULL) #keep geomegry for spatial cv function
dim(thin.occ) #check number of occurrence points in region of interest

#------------------------------------------------------
#------------------------------------------------------
# INTERLUDE INVESTIGATION
# plot species data on the map
e <- raster("~/Desktop/brazil_schisto_raster_folders/sp/sp_hnd.tif")
plot(e)#, xlim=c(-53,-38), ylim=c(-25,-12)) # plot raster data
plot(thin.occ, 
     pch = 16, col=alpha("blue",0.2), 
     add=TRUE) # add absence points
#plot(analysis_data[which(analysis_data$presence==1), ], pch = 16, col= alpha("red",0.2), add=TRUE) # add presence points

# spatial blocking by specified range with random assignment
set.seed(5)
thin.occ.blocking <- thin.occ
sb2 <- blockCV::spatialBlock(speciesData = thin.occ.blocking,
                             species = "source",
                             rasterLayer = e,
                             rows = 5,
                             cols = 5,
                             k = 9)
                             # selection = "systematic",
                             # biomod2Format = TRUE)
thin.occ.blocking$spatialCV_fold <- sb2$foldID
thin.occ.blocking.large.fold <- thin.occ.blocking[which(thin.occ.blocking$spatialCV_fold==5),]
set.seed(106)
rows <- sample(nrow(thin.occ.blocking.large.fold),78)
thin.occ.blocking.sampled.fold <- thin.occ.blocking.large.fold[rows,]
thin.occ <- rbind(thin.occ.blocking[which(thin.occ.blocking$spatialCV_fold!=5),],thin.occ.blocking.sampled.fold)
thin.occ <- subset(thin.occ, select=-c(spatialCV_fold))
#------------------------------------------------------
#------------------------------------------------------

# # IF NEED TO LIMIT OCC DATA FOR PROPER OCC/BG RATIO DO HERE
# # keep "best?" function of predictor variables?
set.seed(115) #100, 110, 11, 5 (4), 50, 4, 10, 12, 13, 20, 200, 27, 33, 55, 66, 77, 88, 99
rows <- sample(nrow(thin.occ), 115) #115 #134 #218
thin.occ <- thin.occ[rows,]
dim(thin.occ)

plot(e)#, xlim=c(-53,-38), ylim=c(-25,-12)) # plot raster data
plot(thin.occ, 
     pch = 16, col=alpha("blue",0.2), 
     add=TRUE) #

# Build extrapolation occ dataset (anything not in region of interest post thinning, could consider not thinning for this set?)
#     NOTE: Depending on the comparison we want to make, we could have it be those outside the region *and* inside,
#     to compare national predicting national to minas gerais predicting national (not minas gerais predicting national-minas gerais?)
extrapolation.occ <- occ_sf_points[apply(occ_by_region,1,isFALSE),]
#for now do not build subregion column for extrapolation points
#extrapolation.occ.subregion <- st_intersects(extrapolation.occ$geometry, subregions$geometry, sparse = FALSE)
#extrapolation.occ$subregion <- apply(extrapolation.occ.subregion, 1, which)
#extrapolation.occ <- st_set_geometry(extrapolation.occ, NULL) #keep geomegry for spatial cv function
dim(extrapolation.occ) #check number of occurrence points outside region of interest

#------------------------------------------------------
# Build background species / background mask rows
#------------------------------------------------------

# Limit to background species (non target species, non wo species specification)
bg_df <- environmental_covariate_df_stor %>%
  filter(!species %in% c(target_species_vector, '', 'sp.'))

# Read in template raster and list setup
rast <- raster("~/Desktop/doctorate/ch1 brazil schisto/schisto_archives/regional_rasters/national_mean_EVI.tif")
bg_species_list <- unique(bg_df$species)

# Extract number of Insecta (+ supplemental) background points per grid cell (i.e., weighted bias mask)
bg_points <- bg_df %>% dplyr::select(c(longitude, latitude)) %>%
  as.matrix()

bg_df$index <- c(1:dim(bg_df)[1])

# Build data set that counts background points by grid cell to build background mask and outputs summary stats
#       about each grid cell grouping
bg_longlat <- cellFromXY(rast, bg_points) %>% as.data.frame() %>%
  cbind(bg_df$year) %>%
  cbind(bg_df$index) %>%
  mutate(count = 1) %>% setNames(c("cell","year","index","count")) %>%
  group_by(cell) %>% dplyr::summarise(count = sum(count),
                                      max_year = max(as.numeric(year)), #to use when choosing land use variable values
                                      avg_year = mean(as.numeric(year)), #to use when choosing land use variable values
                                      max_index = max(index)) %>%
  arrange(desc(count)) %>%
  mutate(longitude = xFromCell(rast, cell),  # Acquire longitude (x) and latitude (y) from cell centroids
         latitude = yFromCell(rast, cell)) %>%
  dplyr::select(-cell) %>% # Cell number is now obsolete, since will be working from (x,y) as an sf object
  filter(!is.na(longitude) & !is.na(latitude)) # Remove the NA locations

# Build geometries
bg_mask_sf_full <- st_as_sf(bg_longlat, coords = c("longitude","latitude"),
                            agr = "constant", crs = 4326) 

# Subset to region of interest (national, regional, state, etc.)
bg_by_region <- st_intersects(bg_mask_sf_full$geometry, geo_extent$x, sparse = FALSE)
bg_sf_points_region <- bg_mask_sf_full[bg_by_region[,1],]

# Build subregion column
bg_mask_sf_subregion <- st_intersects(bg_sf_points_region$geometry, subregions$geometry, sparse = FALSE)
bg_sf_points_region$subregion <- apply(bg_mask_sf_subregion, 1, which)
bg_mask_sf <- bg_sf_points_region
dim(bg_mask_sf) #check number of background points in region of interest

# Random sample bg without replacement from weighted bias mask at (1.5/2x occ) multiplier
set.seed(9)

bg_mask_weights <- bg_mask_sf %>%
  mutate(weight = count/sum(count))

bg_mask_df <- bg_mask_sf[sample(nrow(bg_mask_weights),
                                size = multiplier * nrow(thin.occ),
                                replace = FALSE,
                                prob = bg_mask_weights$weight),]

# Link to environmental covariate data
bg_df_wcovariates <- cbind(bg_df[c(bg_mask_df$max_index),], bg_mask_df$geometry, bg_mask_df$subregion)
colnames(bg_df_wcovariates)[dim(bg_df_wcovariates)[2]] <- 'subregion'

# Build extrapolation bg dataset (anything not in region of interest post thinning, could consider not thinning for this set?)
#       Same discussion points as above. Would also need to change here.
extrapolation.bg <- bg_mask_sf_full[apply(bg_by_region,1,isFALSE),]
dim(extrapolation.bg) #check number of background points outside region of interest
# Don't make subregion for extrapolation datasets yet
#extrapolation.bg.subregion <- st_intersects(extrapolation.bg$geometry, subregions$geometry, sparse = FALSE)
#extrapolation.bg$subregion <- apply(extrapolation.bg.subregion, 1, which)

# Random sample bg without replacement from weighted bias mask at (1.5/2x occ) multiplier
extrapolation.bg.weights <- extrapolation.bg %>%
  mutate(weight = count/sum(count))

extrapolation.bg <- extrapolation.bg[sample(nrow(extrapolation.bg.weights),
                                size = multiplier * nrow(extrapolation.occ),
                                replace = FALSE,
                                prob = extrapolation.bg.weights$weight),]

#extrapolation.bg <- st_set_geometry(extrapolation.bg, NULL)#keep geometry for spatial cv function
extrapolation.bg <- bg_df[c(extrapolation.bg$max_index),]
#question.. do we want to ask it to predict all potential background points or just the multiplier * 
#       the extrapolated occurrences? for now have it set as the latter to be as comparable as possible

#------------------------------------------------------
# Build final input "all_data" set
#------------------------------------------------------

#Build outcome row and only retain predictor columns
thin.occ$presence <- rep(1, nrow(thin.occ))
thin.occ <- data.frame(thin.occ)
#thin.occ <- subset(thin.occ, select=-c(species, year, source))
thin.occ <- thin.occ[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                        'ag_mosiac', 'distance_high_pop', 'hnd', 'soilWater', 'pH', 'clay',
                        'geometry', 'subregion')]#, 'source', 'species')]

bg_df_wcovariates$presence <- rep(0, nrow(bg_df_wcovariates))
bg_df_wcovariates <- data.frame(bg_df_wcovariates)
#bg_df_wcovariates <- subset(bg_df_wcovariates, select=-c(species, year, source, index))
bg_df_wcovariates <- bg_df_wcovariates[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                                          'ag_mosiac', 'distance_high_pop', 'hnd', 'soilWater', 'pH', 'clay',
                                          'geometry', 'subregion')]#, 'source', 'species')]
all_data <- rbind(thin.occ, bg_df_wcovariates)
# Add subregion effects as categorical variable of numbers
all_data$subregion <- as.factor(all_data$subregion)

#sanity check
dim(all_data)
table(all_data$presence)
names(all_data)

saveRDS(all_data, "~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sp_limited.rds")
saveRDS(test, "~/Desktop/schisto_manuscript_figs/jan_mlk_runs/sp_occurences_thinned_mapping_stram_tenag.rds")

test <- readRDS("~/Desktop/schisto_manuscript_figs/jan_mlk_runs/sp_occurences_thinned_mapping.rds")
#test2 <- readRDS("~/Desktop/schisto_manuscript_figs/jan_mlk_runs/mg_straminea_occurences_thinned_mapping.rds")
test <- test[,c(1,12:15)]
thin.occ <- thin.occ[,c(1,12:15)]
test <- rbind(test, thin.occ)
test <- test[which(test$species!="B. straminea"),]
thin.occ[which(thin.occ$species=="glabrata"),5] <- c("B. glabrata")

#------------------------------------------------------
# Build extrapolation dataset for subregional model tests
#------------------------------------------------------

# extrapolation.occ$presence <- rep(1, nrow(extrapolation.occ))
# extrapolation.occ <- extrapolation.occ[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec', 'hnd',
#                         'total_population', 'clay', 'forest', 'fresh_water', 'pasture', 'agriculture')]
# extrapolation.bg$presence <- rep(0, nrow(extrapolation.bg))
# extrapolation.bg <- extrapolation.bg[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec', 'hnd',
#                   'total_population', 'clay', 'forest', 'fresh_water', 'pasture', 'agriculture')]
# extrapolation_data <- rbind(extrapolation.occ, extrapolation.bg)
# # Add subregion effects as categorical variable of numbers (at the moment not in use for extrapolation dataset)
# #extrapolation_data$subregion <- as.factor(extrapolation_data$subregion)

test <- readRDS("~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sp_expert_limited.rds")
head(test)
test <- test[,c(1:13)]
head(test)
saveRDS(test, "~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sp_gbif.rds")
