## Stock assessment configuration
#--------------------------------
SA <- "none"
# "spatialIPM": spatially-explicit model IPM
# "nonspatialIPM": non spatial model IPM
# "GMACS": standard stock assessment model
# "none": no feedback loop

## For spatial MSE
#-----------------
## Spatial settings
Data_Set <- 'Snow_crab'

## Number of knot/stations
n_x = c(30, 50, 75, 100, 150, 200, 300)[3]

# Output from Calc_Kmeans with knots for a triangulated mesh
# Calc_Kmeans determines the location for a set of knots for
# approximating spatial variation
# n_x: the number of knots to select
# nstart the number of times that the k-means algorithm is run while
# searching for the best solution (default=100)
# ter.max the number of iterations used per k-means algorithm
Kmeans_Config = list("randomseed" = 1,
                     "nstart" = 100,
                     "iter.max" = 1e3)

# Define studied region
# builds an object used to determine areas to extrapolation
# densities to when calculating indices
strata.limits <- data.frame('STRATA' = "All_areas")
Region = "Eastern_Bering_Sea"
Extrapolation_List = FishStatsUtils::make_extrapolation_info(Region = Region,
                                                             strata.limits = strata.limits)


## For GMACS
#-----------
# Name of the R script used to incorporate GMACS specifications
Gmacs_Input_file <- "Input_GMACS_MSE"

source(file = file.path(here::here(), "5_GMACS", "functions", 
                        "write_GMACS_files.R"))

