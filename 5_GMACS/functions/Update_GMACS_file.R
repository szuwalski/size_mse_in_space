# Update Assessment model input files
# This updates both the data and the control file
# in the loop

# current env
var.to.save <- ls()

# Basic functions
.an <- function(x) {
  return(as.numeric(x))
}
.ac <- function(x) {
  out <- as.character(x)
  return(out)
}

# Directory to write the files
pathTo <- file.path(dir_GMACS, "model", "build",
                    model_name, fsep = fsep)

# Read the input files
# ---------------------------------------------------------------------------- #
info <- readLines(file.path(dir_GMACS, "scripts_mse",
                            paste0(Gmacs_Input_file,".R"), fsep = fsep))
model_name <- stringr::str_detect(string = info, pattern = "model_name <-")
model_name <- unlist(stringr::str_split(string = info[model_name], 
                                        pattern = "<-", n = 2))[2]

model_name <- base::trimws(
  base::gsub(x = model_name, pattern = "[^%[:alnum:]%s_]+", "")
)
rm(info)

gmacsFiles <- readGMACS.dat(path = file.path(pathTo, "gmacs.dat", fsep = fsep),
                            verbose = FALSE)
Dat <- readGMACSdat(FileName = file.path(pathTo, gmacsFiles$DatFileName, fsep = fsep),
                       verbose = FALSE)

Ctl <- readGMACSctl(FileName = file.path(pathTo, gmacsFiles$CtlFileName, fsep = fsep),
                       verbose = FALSE, 
                       DatFile = Dat, 
                       nyrRetro = gmacsFiles$N_Year_Retro)
# ----------------------------------------------------------------------------

# Fill in the Data file ---- 
# ---------------------------------------------------------------------------- #

# Assessment years
# ======================================= #
# round(t/year_step)
Dat[["End_Y"]] <- .an(max(catch_N$Year))
Dat[["N_year"]] <- Dat$End_Y - Dat$Start_Y + 1

# Catch
# ======================================= #
# Get the material
source(file = file.path(here::here(), 
                        "5_GMACS", "scripts_mse", 
                        paste0(Gmacs_Input_file, ".R")))

# Number of catch data frame in the assessment
# Pot fishery - We split the catch over the fishing seasons
# => 3 data frame
nfleet <- length(F_Fleet_names)
Dat[["N_CatchDF"]] <- nsex * length(which(fish_time_Yr == 1)) * nfleet

# Total catch per year split equally over the number of
# fishing season
seas_Fish <- which(fish_time_Yr == 1) # Fishing seasons

# Total catch per year
catch_Ass <- catch_N %>%
  dplyr::group_by(Year) %>%
  summarise(catch_N = sum(Catches_N))

# Split the catch over the fishing season
if (Split_CatchDF == 0) {
  catch_Ass <- catch_Ass %>%
    dplyr::group_by(Year) %>%
    tidyr::expand(catch_N, seas_Fish) %>%
    dplyr::mutate(catch_N = catch_N / length(seas_Fish))
}

# Split catch_Ass to get one data frame per season
catch_Ass <- split(catch_Ass, catch_Ass$seas_Fish)

# Number of rows per data frame
Dat[["Nrows_CatchDF"]] <- sapply(catch_Ass, nrow)

# Fill the catch
tmpCatch <- NULL

for (n in 1:Dat$N_CatchDF) {
  # n <- 1
  tmp <- NULL
  tmp <- catch_Ass[[n]]
  tmp <- tmp %>% 
    dplyr::mutate(year = Year,
                  seas = seas_Fish,
                  fleet = 1,
                  sex = CatchSex[n],
                  obs = catch_N,
                  CV = CatchCV[n],
                  Type = CatchType[n],
                  units = CatchUnit[n],
                  mult = Catchmult[n],
                  effort = Catcheffort[n],
                  discard_mortality = Discard_Morta[n])
  tmp <- tmp[,-c(1:3)]
  tmp <- tmp %>% mutate_if(is.character, as.numeric)
  Dat[["Catch"]][[n]] <- tmp
  colnames(Dat[["Catch"]][[n]]) <- c(
    "year",
    "seas",
    "fleet",
    "sex",
    "obs",
    "CV",
    "Type",
    "units",
    "mult",
    "effort",
    "discard_mortality"
  )
  names(Dat[["Catch"]])[n] <- paste0(F_Fleet_names[unique(Dat[["Catch"]][[n]]$fleet)], 
                                    "_Seas_", unique(Dat[["Catch"]][[n]]$seas))
  tmpCatch <- rbind(tmpCatch, tmp)
}

# Survey
# ======================================= #
# Season of survey
Seas_Surv <- which(survey_time_Yr==1)

# Total catch per year
Surv_Ass <- Data_Geostat %>%
  dplyr::group_by(year) %>%
  summarise(Data_Geostat = sum(Catch_N))

# Split Surv_Ass to get one data frame per fleet
Surv_Ass$fleet <- 2 # This needs to be changed later if multiple surveys
Surv_Ass <- split(Surv_Ass, Surv_Ass$fleet)

# Set the list of survey data frames 
Dat[["Surveys"]] <- list()
# Number of rows per data frame
Dat[["Nrows_SvDF"]] <- NULL

# Fill in the survey
tmpSurv <- NULL

for(n in 1:Dat$N_SurveyDF){
  tmp <- Surv_Ass[[n]]
  Dat[["Nrows_SvDF"]] <- c(Dat[["Nrows_SvDF"]], nrow(tmp))
  
  tmp <- dplyr::rename(.data = tmp, Abundance= "Data_Geostat")
  tmp <- tmp %>% 
    dplyr::mutate(Index = 1,
                  year = year,
                  seas = Seas_Surv,
                  fleet = fleet,
                  sex = SurvSex[n],
                  Mature = SurvMaturity[n],
                  Abundance = Abundance,
                  CV = SurvCV[n],
                  units = SurvUnit[n],
                  Timing = SurvTiming[n])
  tmp <- tmp %>% 
    dplyr::relocate(any_of(c(
      "Index",
      "year",
      "seas",
      "fleet",
      "sex",
      "Mature",
      "Abundance",
      "CV",
      "units",
      "Timing"
    )))
  tmp <- tmp %>% mutate_if(is.character, as.numeric)
  Dat[["Surveys"]][[n]] <- tmp
  names(Dat[["Surveys"]])[n] <- Dat$Survey_names[n]
  tmpSurv <- rbind(tmpSurv, tmp)
}

# Size Frequency data
# ======================================= #

# Catch length comps
SizeCps <- catch_N %>%
  dplyr::group_by(Year, size_bin) %>%
  summarise(tot_catch = sum(Catches_N))


SizeCps <- SizeCps %>%
  group_by(Year) %>%
  mutate(freq = tot_catch/sum(tot_catch)) %>%
  dplyr::select(size_bin, freq) %>%
  pivot_wider(names_from = size_bin,
              values_from = freq) %>%
  replace(is.na(.), 0)

SizeCps <- as.data.frame(cbind(
  SizeCps[,1],
  SizeCps[,-1][,order(as.numeric(colnames(SizeCps[,-1])))]
))

SizeFreq <- list()

for(n in 1:N_SizeFreq_df[1]){
  
  #_Year_| Season_| Fleet_| Sex_| Type_| Shell_| Maturity_| Nsamp_|
  
  
  SizeFreq[[n]] <- SizeCps %>%
    mutate(year = Year,
           seas = unique(Dat[["Catch"]][[n]]$seas),
           fleet = unique(Dat[["Catch"]][[n]]$fleet),
           sex = 1,
           Type = 1,
           Shell = 0,
           Maturity = 0,
           Nsamp = Nsam_SizeC[n])
  SizeFreq[[n]]<- SizeFreq[[n]][,colnames(SizeFreq[[n]])!="Year"]
  colnamSizeC <- c(
    "year",
    "seas",
    "fleet",
    "sex",
    "Type",
    "Shell",
    "Maturity",
    "Nsamp",
    1:Dat[["Nbins_SiseFreq"]][n])
  
  SizeFreq[[n]] <- SizeFreq[[n]] %>% 
    mutate_if(is.character, as.numeric)
  
  Dat[["Nrows_SiseFreqDF"]][n] <- dim(SizeFreq[[n]])[1]
  
  Dat[["SizeFreq"]][[n]] <- SizeFreq[[n]] %>% relocate(any_of(colnamSizeC))
  
  names(Dat[["SizeFreq"]])[n] <- paste0(F_Fleet_names[unique(Dat[["SizeFreq"]][[n]]$fleet)], 
                                        "_Seas_", unique(Dat[["Catch"]][[n]]$seas))
}

# Survey Length Comps
SizeCps <- Data_Geostat %>%
  dplyr::group_by(year, size_class) %>%
  summarise(tot_catch = sum(Catch_N))


SizeCps <- SizeCps %>%
  group_by(year) %>%
  mutate(freq = tot_catch/sum(tot_catch)) %>%
  dplyr::select(size_class, freq) %>%
  pivot_wider(names_from = size_class,
              values_from = freq) %>%
  replace(is.na(.), 0)


SizeCps <- as.data.frame(cbind(
  SizeCps[,1],
  SizeCps[,-1][,order(as.numeric(colnames(SizeCps[,-1])))]
))

SizeFreq <- list()

for(n in (N_SizeFreq_df[1]+1):(N_SizeFreq_df[1]+N_SizeFreq_df[2])){
  
  p <- n - N_SizeFreq_df[1]
  
  SizeFreq[[n]] <- SizeCps %>%
    mutate(seas = unique(Dat[["Surveys"]][[p]]$seas),
           fleet = unique(Dat[["Surveys"]][[p]]$fleet),
           sex = 1,
           Type = 1,
           Shell = 0,
           Maturity = 0,
           Nsamp = Nsam_SizeC[n])
  colnamSizeC <- c(
    "year",
    "seas",
    "fleet",
    "sex",
    "Type",
    "Shell",
    "Maturity",
    "Nsamp",
    1:Dat[["Nbins_SiseFreq"]][n])
  
  SizeFreq[[n]] <- SizeFreq[[n]] %>% 
    mutate_if(is.character, as.numeric)
  
  Dat[["Nrows_SiseFreqDF"]][n] <- dim(SizeFreq[[n]])[1]
  
  Dat[["SizeFreq"]][[n]] <- SizeFreq[[n]] %>% relocate(any_of(colnamSizeC))
  
  names(Dat[["SizeFreq"]])[n] <- paste0(Survey_names[p], 
                                        "_Seas_", unique(Dat[["Surveys"]][[p]]$seas))
}

# update N_SizeFreq_df
Dat[["N_SizeFreq_df"]] <- sum(N_SizeFreq_df)

# Write the data file
# ======================================= #
# writeGmacsdatfile(
#   Dir = pathTo,
#   FileName = gmacsFiles$DatFileName,
#   overwrite = TRUE,
#   DatFile = Dat,
#   stock = "Snow crab",
#   model_name = model_name,
#   Ass_Year = Dat[["End_Y"]]
# )
# ----------------------------------------------------------------------------


# Fill in the Control file ---- 
# ---------------------------------------------------------------------------- #

#  Get the size transition matrix
# ======================================= #
# Probably not a good idea to allow time varying in GMACS - Need to check hypotheses
if (Ctl[["bUseCustomGrowthMatrix"]] == 1 ||
    Ctl[["bUseCustomGrowthMatrix"]] == 2) {
  # GROWTH_FIXEDSIZETRANS
  if (Ctl[["bUseCustomGrowthMatrix"]] == 2)
    cat("\t!! Considering the male size transition matrix for all individuals.\n")
  if(nsex == 1){
    Ctl[["CustomGrowthMatrix"]] <- size_transition_mat_m_imm
  } else {
    Ctl[["CustomGrowthMatrix"]] <- rbind(size_transition_mat_m_imm,
                                            size_transition_mat_m_imm)
  }
} else {
  Ctl[["CustomGrowthMatrix"]] <- ""
}

# Update the time block in selectivity and retention control parameter 
# matrix
# ======================================= #
tmp <- cbind(
  c(F_Fleet_names, Survey_names), 
  1:length(c(F_Fleet_names, Survey_names))
  )
for(i in 1:dim(tmp)[1]){
  
  if(tmp[i,1] %in% F_Fleet_names){
    tmpMin_Y <- .an(min(tmpCatch[tmpCatch$fleet == .an(tmp[i,2]), "year"]))
    tmpMax_Y <- .an(max(tmpCatch[tmpCatch$fleet == .an(tmp[i,2]), "year"]))
  } else if (tmp[i,1] %in% Survey_names){
    tmpMin_Y <- .an(min(tmpSurv[tmpSurv$fleet == .an(tmp[i,2]), "year"]))
    tmpMax_Y <- .an(max(tmpSurv[tmpSurv$fleet == .an(tmp[i,2]), "year"]))
  }
    # Selectivity
  Ctl[["Selex_control"]][Ctl[["Selex_control"]]$Fleet == .an(tmp[i,2]),] <- Ctl[["Selex_control"]] %>% 
      dplyr::filter(Fleet == .an(tmp[i,2])) %>%
      dplyr::mutate(Start_Block = tmpMin_Y) %>%
      dplyr::mutate(End_Block = tmpMax_Y)
    # Retention
  Ctl[["Ret_control"]][Ctl[["Ret_control"]]$Fleet == -1*.an(tmp[i,2]),] <- 
    Ctl[["Ret_control"]] %>% 
      dplyr::filter(Fleet == -1*.an(tmp[i,2])) %>%
      dplyr::mutate(Start_Block = tmpMin_Y) %>%
      dplyr::mutate(End_Block = tmpMax_Y)
}

# Other addition controls
# ======================================= #
Ctl[["rdv_syr"]] <- StartY_Ass  # First year of recruitment estimation deviations
Ctl[["rdv_eyr"]] <- Dat$End_Y # Last year of recruitment estimation deviations
# Weights on catches for the likelihood component
Ctl[["catch_emphasis"]] <- rep(1,Dat[["N_CatchDF"]])

# Write the updated control file
# ======================================= #
# writeGmacsctlfile(Dir = pathTo,
#                   FileName = gmacsFiles$CtlFileName,
#                   DatFile = Dat,
#                   CtlFile = Ctl,
#                   stock = "Snow crab",
#                   model_name = model_name,
#                   Ass_Year = Dat$End_Y)
# ----------------------------------------------------------------------------



# Write the updated files

# Update and Write the gmacs.dat file
# ======================================= #

DatFileName <- paste0("snow_Yr",round(t/12),".dat")
CtlFileName <- paste0("snow_Yr",round(t/12),".ctl")
gmacsFiles$DatFileName <- DatFileName
gmacsFiles$CtlFileName <- CtlFileName

writeGmacs.dat(
  Dir = pathTo,
  FileName = "gmacs.dat",
  gmacsDat = gmacsFiles,
  stock = "Snow crab",
  model_name = model_name,
  Ass_Year = Dat$End_Y
)



# Write the data file
# ======================================= #
writeGmacsdatfile(
  Dir = pathTo,
  FileName = gmacsFiles$DatFileName,
  overwrite = TRUE,
  DatFile = Dat,
  stock = "Snow crab",
  model_name = model_name,
  Ass_Year = Dat[["End_Y"]]
)
# Write the updated control file
# ======================================= #
writeGmacsctlfile(Dir = pathTo,
                  FileName = gmacsFiles$CtlFileName,
                  DatFile = Dat,
                  CtlFile = Ctl,
                  stock = "Snow crab",
                  model_name = model_name,
                  Ass_Year = Dat$End_Y)

# Restore the environment as it was but saved some stuff needed in the loop
rm(list = setdiff(ls(), var.to.save))
