# Assessment model set up for GMACS

# Local declarations ----
fsep <- .Platform$file.sep

#  I. Folder/ files set up ----
# ---------------------------------------------------------------------------- #
# Path to the base model input files (the ones to be modified)
dir_GMACS <- file.path(here::here(), "5_GMACS", fsep = fsep)
pathFrom <- file.path(dir_GMACS, "model", "build",
                      "assessment_files", fsep = fsep)

# Name of the model/MSE.
# This is used to create a MSE specific folder where GMACS will be run
# and the outputs saved - Specifically, results will be saved in
# ~/5_GMACS/model/build/model_name
model_name <- "MSE_1"

# Overwrite input files?
overwrite <- TRUE

# Print statement
verbose <- TRUE
# ----------------------------------------------------------------------------

# II. Data file related stuf ----
# ---------------------------------------------------------------------------- #
# File names
# ======================================= #
# Name of the input data file to be written
datfilenameTo <- "snow.dat"

# Main set up
# ======================================= #
StartY_Ass <- 1989             # Start year for the assessment period
nsex <- 1                      # Number of sex to consider in the assessment
nMaturity <- 2                 # Number of maturity state
nShell <- 1                    # Number of shell condition state
F_Fleet_names <- "Pot_Fishery" # Name(s) for the fishing fleet(s)
Survey_names <- "Survey_1"     # Name(s) for the survey fleet(s)
N_seas <- 1                    # Season to comput N output

# Natural mortality
# ======================================= #
M_type <- 1

# Catch
# ======================================= #
# How to split the catches per season?
# Currently one option (equally split between seasons)
# but we can imagine to set up a distribution law on this ??
Split_CatchDF  <- 0

# Catch description for each data frame
CatchSex <-      c(0, 0, 0)          # Set as unsexed
CatchCV <-       c(0.01, 0.01, 0.01) # Catch CVs
CatchType <-     c(0, 0, 0)          # Total
CatchUnit <-     c(2, 2, 2)          # Numbers
Catchmult <-     c(1, 1, 1)          # Use them as data
Catcheffort <-   c(0, 0, 0)          # Effort
Discard_Morta <- c(1, 1, 1)          # Discard mortality

# Survey
# ======================================= #
# Number of survey data frame
N_SurveyDF <- 1

# Data type for each abundance index
# (1: total selectivity; 2:retention*selectivity)
Sv_type <- rep(1, N_SurveyDF)

# Survey description
SurvSex <- c(1)      # Set as unsexed
SurvMaturity <- c(0) # Surv maturity (0: combined)
SurvCV <- c(0.1)     # Survey CV
SurvUnit <- c(2)     # Numbers
SurvTiming <- c(0)   # Survey timing

# Size composition data
# ======================================= #
# Number of size frequency data frame
# c(fisheries, surveys)
N_SizeFreq_df <- c(3,1)

# Nsamp for each data frame
# Here this will be 4 data frame: 
# 3 for each season of the fishery and one for the survey
# This allow to play a bit
Nsam_SizeC <- c(100,100,100,100)

# Growth increment data
# ======================================= #
# The input data depends upon the GrowthObsType (0,1,2,3)
# 0: No growth data observation
# 1: Growth increment data observation
# 2: Data about growth increment based on CMR data (size-at-release; size-at-recaptures, and time-at-liberty)
# 3: Values of growth increment for each size class
GrowthObsType <- 0
  # If GrowthObsType > 0
# NGrowthObs <- NA # Number of growth observation

# ENVIRONMENTAL DATA
# ======================================= #
# Number of environmental indices
NenvIndics <- 0
  # If NenvIndics > 0
# EnvYrs <- NA # Year ranges for each index (One line per index - From 1 to N env indices)

# ---------------------------------------------------------------------------- #

#  III. Control file related stuff ----
# ---------------------------------------------------------------------------- #
# File names
# ======================================= #
# Name of the input data file to be written
ctlfilenameTo <- "snow.ctl"

# Get the estimated values for the parameters from the last assessment
# and set those values in the MSE
# ======================================= #
# This goes for theta, growth increment, molt proba, selectivity, retention
useEstParamVals <- c(TRUE, TRUE, TRUE, TRUE, TRUE)

# Turn off the estimation for all or some parameters
# ======================================= #
TurnOffest <- c(TRUE, TRUE, TRUE, TRUE, TRUE)

# How to initialize at unfished conditions
# ======================================= #
# 1 => unfished
# 2 => Free parameters
# 3 => Free parameters scaled
InitializeUnfished <- 3

# Length-weight relationship type
# ======================================= #
# 1 => LW relationship
# 2 => vector of data
# 3 => matrix of data
  # For options 2 & 3 set the vector/matrix as
  # wt_at_len_U if 1 sex
  # wt_at_len if 2 sex
lw_type <- 1

# Growth matrix
# ======================================= #
# Here the size transition matrix is computed outside GMACS and set as an input
# Then the probability of molting is ignore in GMACS
CustomGrowthMatrix <- 2
# Growth increment model
GrowthIncrementModel <- 1 # Linear growth model
# Number of size Inc periods (can have multiple for each sex)
nSizeIncBlock <- rep(1, nsex)
  # If multiple blocks - Indicate years of blocks for each sex (1 line per sex)
# SizeIncBlock_Yrs <- NULL

# Probability of molting model -
# ======================================= #
# Use a constant molting prob (set equal to 1) - Identical for both sexes
CustomMoltProbability <- 1
# Number of molt periods (can have multiple for each sex)
nMoltBlock <- rep(1, nsex)
# If multiple blocks - Indicate years of blocks for each sex (1 line per sex)
# MoltBlock_Yrs <- NULL
# ----------------------------------------------------------------------------
