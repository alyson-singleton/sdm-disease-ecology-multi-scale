#------------------------------------------------------
#------------------------------------------------------
### Script used to create preprocessed_datasets above. 
# Posted for reproducibility. Recommend to just use the 
# preprocessed_datasets in the data folder, rather than 
# recreate from scratch.
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

# Load RDS file
environmental_covariate_df <- readRDS(file = "../data/full_dataset.rds")

# Only keep variables of interest
variable_list <- c('bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                   'ag_mosiac', 'distance_high_pop',
                   'hnd', 'soilWater', 'pH', 'clay',
                   'species', 'year', 'latitude', 'longitude', 'source')
environmental_covariate_df <- environmental_covariate_df[,variable_list]

# Catch all to remove any incomplete rows
environmental_covariate_df <- environmental_covariate_df[complete.cases(environmental_covariate_df), ]
environmental_covariate_df_stor <- environmental_covariate_df

#------------------------------------------------------
# CHOICE 1: Which occurrence data to use? All? Expert only? GBIF only?
#------------------------------------------------------

# Pick GBIF or expert if want to seperate, otherwise leave commented out
# environmental_covariate_df <- environmental_covariate_df %>%
#   filter(source == "expert") #options: expert, gbif

#------------------------------------------------------
# CHOICE 2: Choose the geographic extent of interest.  Store as "geo_extent" which is what will be
#           referenced throughout the rest of the file. Also store appropriate subregion geometry set to 
#           use later. Brazil shapefiles pulled from Brazilian govern website.
#------------------------------------------------------

### National
geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Pais_2021/BR_Pais_2021.shp")
geo_extent_shp <- st_union(geo_extent_shp)
subregions <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_UF_2021/BR_UF_2021.shp") #state as subregion
subregions$geometry <- st_transform(subregions$geometry, 4326)
geo_extent <- st_as_sf(geo_extent_shp) 
geo_extent$x <- st_transform(geo_extent$x, 4326)

# ### Regional (Sudeste)
# geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_UF_2021/BR_UF_2021.shp")
# subregions <- geo_extent_shp[which(geo_extent_shp$NM_REGIAO=="Sudeste"),] #state as subregion
# subregions$geometry <- st_transform(subregions$geometry, 4326)
# geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$NM_REGIAO=="Sudeste"),])
# geo_extent <- st_as_sf(geo_extent_shp) 
# geo_extent$x <- st_transform(geo_extent$x, 4326)
# 
# ### State (Minas Gerais)
# geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
# subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),] #mesoregion as subregion
# subregions$geometry <- st_transform(subregions$geometry, 4326)
# geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="MG"),])
# geo_extent <- st_as_sf(geo_extent_shp) 
# geo_extent$x <- st_transform(geo_extent$x, 4326)
# 
# ### State (Sao Paulo)
# geo_extent_shp <- read_sf("~/Desktop/doctorate/ch1 brazil schisto/brazil_shapefiles/BR_Mesorregioes_2021/BR_Mesorregioes_2021.shp")
# subregions <- geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),] #mesoregion as subregion
# subregions$geometry <- st_transform(subregions$geometry, 4326)
# geo_extent_shp <- st_union(geo_extent_shp[which(geo_extent_shp$SIGLA=="SP"),])
# geo_extent <- st_as_sf(geo_extent_shp) 
# geo_extent$x <- st_transform(geo_extent$x, 4326)

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
dim(thin.occ) #check number of occurrence points in region of interest

#------------------------------------------------------
# If need to limit no. of occurrences for Figure E in S1 text uncomment next four lines:
# set.seed(115) 
# rows <- sample(nrow(thin.occ), 115)
# thin.occ <- thin.occ[rows,]
# dim(thin.occ)

#------------------------------------------------------
# Build background species / background mask rows
#------------------------------------------------------
# Limit to background species (non target species, non wo species specification)
bg_df <- environmental_covariate_df_stor %>%
  filter(!species %in% c(target_species_vector, '', 'sp.'))

# Read in template raster and list setup
rast <- raster("~/Desktop/doctorate/ch1 brazil schisto/schisto_archives/regional_rasters/national_mean_EVI.tif")
bg_species_list <- unique(bg_df$species)

# Extract number of possible background points per grid cell (i.e., weighted bias mask)
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

# Random sample bg without replacement from weighted bias mask at 2x occ multiplier
multiplier <- 2

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

#------------------------------------------------------
# Build final input "all_data" set
#------------------------------------------------------

#Build outcome row and only retain predictor columns
thin.occ$presence <- rep(1, nrow(thin.occ))
thin.occ <- data.frame(thin.occ)
thin.occ <- thin.occ[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                        'ag_mosiac', 'distance_high_pop', 'hnd', 'soilWater', 'pH', 'clay',
                        'geometry', 'subregion')]

bg_df_wcovariates$presence <- rep(0, nrow(bg_df_wcovariates))
bg_df_wcovariates <- data.frame(bg_df_wcovariates)
bg_df_wcovariates <- bg_df_wcovariates[,c('presence','bio4_temp_', 'bio11_temp', 'bio16_prec', 'bio17_prec',
                                          'ag_mosiac', 'distance_high_pop', 'hnd', 'soilWater', 'pH', 'clay',
                                          'geometry', 'subregion')]
all_data <- rbind(thin.occ, bg_df_wcovariates)

# Add subregion effects as categorical variable of numbers
all_data$subregion <- as.factor(all_data$subregion)

# Sanity check
dim(all_data)
table(all_data$presence)
names(all_data)

# Save file
saveRDS(all_data, "~/Desktop/doctorate/ch1 brazil schisto/late_feb_runs/datasets/tenagophila_sp_limited.rds")