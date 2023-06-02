## Survey and commercial catch data
#----------------------------------

# Standard deviation of survey data
sigma_survey = 0.5

# Catchability survey
q_survey = 1

# I load these datasets from Maximes spatial IPM
# this way, survey data and commercial data are 
# in the good format to feed spatial IPM

## Survey data
#-------------
## Load survey data for sampling locations
load('4_full_MSE/data/Data_Geostat_4class.RData')

DF <- as_tibble(Data_Geostat)

Data_Geostat = cbind(
  "size_class" = DF[, "size_bin"],
  "year" = DF[, "Year"],
  "Catch_N" = DF[, "Catch"],
  "AreaSwept_km2" = DF[, "AreaSwept_km2"],
  "Vessel" = 0,
  "Lat" = DF[, "Lat"],
  "Lon" = DF[, "Lon"]
)

colnames(Data_Geostat) <-
  c("size_class",
    "year",
    "Catch_N",
    "AreaSwept_km2",
    "Vessel",
    "Lat",
    "Lon")

Data_Geostat = Data_Geostat %>%
  filter(year %in% 2015:2016)

## Commercial catch data
#-----------------------
load(paste0(project_spatialIPM,"02_transformed_data/movement/COG/COG_smoother_ADFG/catch_fishery_mov_intersect.RData"))
catch_N <- catch_fishery_mov_intersect

# Choose if we implement fisheries catches (Fisheries_catches <-TRUE) or not in the model 
Fisheries_catches <-TRUE
mov <- FALSE

# Choose if we implement simulated fisheries catches(Catches_sim_random <-TRUE) ro not in the model 
Catches_sim_random <- FALSE
Catches_sim_fractionAb <- FALSE

heterMov <- FALSE
homoMov <- FALSE

# Smoother
#smoother <- TRUE # smoother is with knots
#ADFG <- TRUE # smoother is with ADFG cells
#KNOT_i <- FALSE
