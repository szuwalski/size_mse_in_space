# Basic functions
.an <- function(x) {
  return(as.numeric(x))
}
.ac <- function(x) {
  out <- as.character(x)
  return(out)
}
.av <- function(x) {
  out <- as.vector(x)
  return(out)
}


#' Function to write the base input data file for a MSE
#'
#' @param pathfrom (character string)- Path to the directory of the input data
#' file to be modified.
#' @param pathTo (character string)- path to the directory where to save the new
#' input control file.
#' @param filename (character string)- Name of the input data file to be modified
#' @param model_name (character string)- Name for the model/MSE. This will be use
#' to create the folder (if it does not already exist) where GMACS outputs will be
#' housed.
#' @param overwrite (logical)- Do you want to overwrite the \code{filename} file.
#' @param filenameTo (character string)- Name of the new input data file to write.
#' @param verbose (logical)- Print statement.
#'
#' This function write the input data file for a given configuration of a MSE
set_GMACSdataFile <- function(pathFrom = NULL,
                              pathTo = NULL,
                              model_name = NULL,
                              filename = NULL,
                              overwrite = NULL,
                              filenameTo = NULL,
                              verbose = NULL) {
  # Local declarations
  fsep <- .Platform$file.sep

  # Check Input
  if (pathFrom == pathTo) {
    if (!overwrite && (is.null(filenameTo) || filenameTo == filename)) {
      message(
        "=> Something is wrong with the set up to write the data input file.\n
You can either:
\t1) provide a different name for the 'pathTo' directory
\t2) allow overwriting the input data file in the 'pathFrom' directory, or
\t3) privide a name (or a different one if already provided) for 'filenameTo'."
      )
      stop()
    }
    if (overwrite) {
      cat("=> The input data file in the 'pathFrom' directory will be overwritten !\n")
    }
  }
  if (is.null(pathTo) && !overwrite) {
    message(
      "=> Please provide the 'pathTo' directory or
allow overwriting the input data file in the 'pathFrom' directory."
    )
    stop()
  }
  if (is.null(pathTo) && overwrite) {
    cat("=> The input data file in the 'pathFrom' directory will be overwritten !\n")
  }
  
  # Create pathTo if needed
  if (!is.null(pathTo) & !dir.exists(pathTo)) {
    if (verbose)
      cat("=> Creating the folowing directory:\n", pathTo, "\n")
    dir.create(path = pathTo)
  }
  
  # Read a given data file to be modified
  Dat <- readGMACSdat(FileName = file.path(pathFrom, filename),
                      verbose = verbose)
  DatOut <- Dat
  rm(Dat)
  
  # Model dimensions
  # ======================================================================= #
  # -------------------------------------------------------------------------
  DatOut[["Start_Y"]] <- StartY_Ass # Start year
  DatOut[["End_Y"]] <- Start_Y - 1 # End year
  DatOut[["N_year"]] <- DatOut$End_Y - DatOut$Start_Y + 1
  DatOut[["N_seasons"]] <- year_step # Number of season
  DatOut[["N_fleet"]] <-
    length(c(F_Fleet_names, Survey_names)) # Number of fleets (fishing fleets and surveys)
  DatOut[["N_sexes"]] <- nsex # Number of sexes
  DatOut[["N_shell_cdt"]] <-
    nShell # Number of shell condition types (not considered)
  DatOut[["N_maturity"]] <- nMaturity # Number of maturity types
  DatOut[["N_sizeC"]] <- n_p # Number of size class in the model
  
  
  # Timing
  # Only one season is available for each biological process
  # We consider the first season within the year
  DatOut[["Recr_Season"]] <-
    which(recruit_time_Yr == 1)[1] # Season recruitment occurs
  DatOut[["Grwth_Season"]] <-
    which(Fem_molt_time_Yr == 1)[1] # Season molting & growth occur
  DatOut[["SSB_Season"]] <-
    which(mate_time_Yr == 1)[1] # Season to calculate SSB
  DatOut[["N_Season"]] <- N_seas # Season for N output
  # ======================================================================= #
  # -------------------------------------------------------------------------
  
  # Size classes definition
  # -------------------------------------------------------------------------
  if (DatOut[["N_sexes"]] == 1) {
    # Maximum size class
    DatOut[["Max_sizeC"]] <- n_p # for males only
  } else {
    DatOut[["Max_sizeC"]] <-
      rep(n_p, DatOut[["N_sexes"]]) # for males and then females
  }
  DatOut[["Size_breaks"]] <- binclass # Size breaks
  # vector giving the break points between size intervals
  # -------------------------------------------------------------------------
  
  # Natural mortality
  # -------------------------------------------------------------------------
  DatOut[["M_in_Type"]] <-
    M_type # Natural mortality per season input type
  # If M_in_Type == 1 vector by season
  # If M_in_Type == 2 matrix by year/season => matrix[N_year, N_seasons]
  DatOut[["M_Seas_prop"]] <-
    rep(1 / DatOut[["N_seasons"]], DatOut[["N_seasons"]])
  # -------------------------------------------------------------------------
  
  # Fishery and survey definition
  # -------------------------------------------------------------------------
  DatOut[["F_Fleet_names"]] <- F_Fleet_names # Fishing fleet names
  DatOut[["Survey_names"]] <- Survey_names # Survey fleet names
  
  # Type of fishing mortality per season (0: instantaneous; 1: continuous)
  DatOut[["F_Season_Type"]] <- rep(1, DatOut[["N_seasons"]])
  # -------------------------------------------------------------------------
  
  # Catch data
  # -------------------------------------------------------------------------
  nCatchDF <- 1
  DatOut[["N_CatchDF"]] <- nCatchDF # Number of catch data frame (default)
  DatOut[["Nrows_CatchDF"]] <-
    rep(1, nCatchDF) # Number of rows in the catch data frame
  
  DatOut[["Catch"]] <- list()
  namcatch <- NULL
  colnamCatch <- c(
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
  for (n in 1:DatOut$N_CatchDF) {
    # n <- 1
    DatOut[["Catch"]][[n]] <-
      matrix(-1, nrow = 1, ncol = length(colnamCatch))
    colnames(DatOut[["Catch"]][[n]]) <- colnamCatch
  }
  names(DatOut[["Catch"]]) <-
    paste0("Initial_Empty_catch_DF", 1:nCatchDF)
  # -------------------------------------------------------------------------
  
  # Relative abundance index
  # -------------------------------------------------------------------------
  # Number of Survey data frames
  DatOut[["N_SurveyDF"]] <- length(DatOut$Survey_names)
  
  if (DatOut$N_SurveyDF != length(DatOut$Survey_names)) {
    cat(
      "The number of relative abundance indices does not correspond to the
number of survey names. GMACS will consider CPUE index! \n"
    )
  }
  
  # Data type for each abundance index
  # (1: total selectivity; 2:retention*selectivity)
  DatOut[["Sv_type"]] <- Sv_type
  # Number of rows of index data
  DatOut[["Nrows_SvDF"]] <- 1
  
  DatOut[["Surveys"]] <- list()
  colnamSurv <- c(
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
  )
  for (n in 1:DatOut$N_SurveyDF) {
    # n <- 1
    DatOut$Surveys[[n]] <-
      matrix(-1, nrow = 1, ncol = length(colnamSurv))
    colnames(DatOut$Surveys[[n]]) <- colnamSurv
    names(DatOut[["Surveys"]])[n] <- "Initial_Empty_Survey"
  }
  # -------------------------------------------------------------------------
  
  
  # SIZE COMPOSITION DATA FOR ALL FLEETS
  # -------------------------------------------------------------------------
  #Number of length frequency matrix
  if (sum(N_SizeFreq_df) > 0) {
    DatOut[["Nrows_SiseFreqDF"]] <-  rep(1, sum(N_SizeFreq_df))
    DatOut[["Nbins_SiseFreq"]] <-  rep(DatOut$N_sizeC, sum(N_SizeFreq_df))
    colnamSizeFreq <- c(
      "year",
      "seas",
      "fleet",
      "sex",
      "Type",
      "Shell",
      "Maturity",
      "Nsamp",
      rep("", DatOut[["Nbins_SiseFreq"]][n]))
    
    DatOut[["SizeFreq"]] <- list() # Length frequency matrices
    for (n in 1:sum(N_SizeFreq_df)) {
      DatOut[["SizeFreq"]][[n]] <-
        matrix(-1, nrow = 1, ncol = length(colnamSizeFreq))
      names(DatOut[["SizeFreq"]])[n] <-
        paste0("Initial_Empty_SizeFrqDF", n)
    }
  }
  DatOut[["N_SizeFreq_df"]] <- sum(N_SizeFreq_df)
  # -------------------------------------------------------------------------
  
  # GROWTH DATA
  # -------------------------------------------------------------------------
  # Type of observation (increment or change in size-class)
  DatOut[["GrowthObsType"]] <- GrowthObsType
  
  if (DatOut[["GrowthObsType"]] > 3)
    stop("GrowthObsType can only be 0,1, 2 or 3")
  
  if (DatOut[["GrowthObsType"]] > 0) {
    # Number of observation (lines of the df)
    DatOut[["NGrowthObs"]] <- 0
    
    if (DatOut[["NGrowthObs"]] > 0) {
      namGrowthObs <- base::switch(
        .ac(DatOut$GrowthObsType),
        "1" = c("Premolt", "Sex", "Molt_Inc", "CV"),
        "2" = c("Size_rel", "Size_Recap", "T_at_sea"),
        "3" = c(
          "Size_rel",
          "sex",
          "Size_Recap",
          "T_at_sea",
          "fleet",
          "Recap_Year",
          "Number"
        )
      )
      DatOut[["GrowthData"]] <-
        matrix(-1, nrow = 1, ncol = length(namGrowthObs))
      colnames(DatOut[["GrowthData"]]) <- namGrowthObs
    }
  }
  # -------------------------------------------------------------------------
  
  # ENVIRONMENTAL DATA
  # -------------------------------------------------------------------------
  # Number of environmental indices
  DatOut[["NenvIndics"]] <- NenvIndics
  
  if (DatOut[["NenvIndics"]] > 0) {
    DatOut[["EnvYrs"]] <- EnvYrs
    colnamEnv <- c("Index", "Year", "Value")
    DatOut[["EnvData"]] <-
      matrix(-1, nrow = 1, ncol = length(colnamEnv))
  }
  # ------------------------------------------------------------------------------
  
  # Write the new GMACS data input file
  if(pathFrom == file.path(dir_GMACS, "model", "build",
                           "assessment_files", fsep = fsep))
    overwrite <- FALSE
  
  if (overwrite && is.null(pathTo)) {
    Dir <- pathFrom
    FileName <- filename
  } else {
    Dir <- pathTo
    if (!is.null(filenameTo)) {
      FileName <- filenameTo
    } else {
      FileName <- filename
    }
  }
  if (verbose)
    cat("=> Writing the input data file in the following directory:\n",
        Dir,
        "\n")
  writeGmacsdatfile(
    Dir = Dir,
    FileName = FileName,
    overwrite = overwrite,
    DatFile = DatOut,
    stock = "Snow crab",
    model_name = model_name,
    Ass_Year = ""
  )
}


#' Function to write the base input control file for a MSE
#'
#' @param pathfrom (character string)- Path to the directory of the input control
#' file to be modified.
#' @param pathTo (character string)- path to the directory where to save the new
#' input control file.
#' @param filename (character string)- Name of the input control file to be modified
#' @param model_name (character string)- Name for the model/MSE. This will be use
#' to create the folder (if it does not already exist) where GMACS outputs will be 
#' housed.
#' @param overwrite (logical)- Do you want to overwrite the \code{filename} file 
#' in the \code{model_name} folder.
#' @param filenameTo (character string)- Name of the new input control file to write.
#' @param verbose (logical)- Print statement.
#'
#' This function write the input data file for a given configuration of a MSE
set_GMACSctlFile <- function(pathFrom = NULL,
                             pathTo = NULL,
                             filename = NULL,
                             model_name = NULL,
                             overwrite = NULL,
                             filenameTo = NULL,
                             verbose = NULL) {
  
  # Local declarations
  fsep <- .Platform$file.sep

  # Check Input
  if(pathFrom == pathTo){
    if(!overwrite && (is.null(filenameTo) || filenameTo == filename)){
      message("=> Something is wrong with the set up to write the control input file.\n
You can either:
\t1) provide a different name for the 'pathTo' directory (i.e., another model name)
\t2) allow overwriting the input control file in the 'pathFrom' directory, or 
\t3) privide a name (or a different one if already provided) for 'filenameTo'.")
      stop()
    }
    if(overwrite){
      cat("=> The input control file in the 'pathFrom' directory will be overwritten !\n")
    }
  }
  if(is.null(pathTo) && !overwrite){
    message("=> Please provide the 'pathTo' directory or 
allow overwriting the input control file in the 'pathFrom' directory.")
    stop()
  }
  if(is.null(pathTo) && overwrite){
    cat("=> The input control file in the 'pathFrom' directory will be overwritten !\n")
  }
  
  # Create pathTo if needed
  if(!is.null(pathTo) & !dir.exists(pathTo)){
    cat("=> Creating the folowing directory:\n",pathTo,"\n")
    dir.create(path = pathTo)
  }
  
  # Updating file or creating from the Assessment files
  FromAssessFile <- FALSE
  if(stringr::str_detect(string = pathFrom, pattern = "assessment_files"))
    FromAssessFile <- TRUE
  
  # Read the GMACS.dat file to get file names 
  gmacsFiles <-
    readGMACS.dat(path = file.path(pathFrom, "gmacs.dat", fsep = fsep),
                  verbose = FALSE)
  Dat <-
    readGMACSdat(
      FileName = file.path(pathFrom, gmacsFiles$DatFileName, fsep = fsep),
      verbose = FALSE
    )
  
  # Load the current data file for snow crab
  Ctl <- readGMACSctl(
    FileName = file.path(pathFrom, filename),
    DatFile = Dat,
    nyrRetro = 0,
    verbose = verbose
  )
  
  CtlOut <- Ctl
  
  # Read the gmacsall.out file if useEstParamVals == TRUE
  # =============================================== #
  # to get the estimated values of the parameters from the last assessment
  if ((length(unique(useEstParamVals)) > 1 &&
       unique(useEstParamVals)[1]) ||
      unique(useEstParamVals)) {
    cat("=> reading Gmacsall.out to get the estimates of the parameters
...... \n")
    AllOut <-
      readGMACSallOUT(
        FileName = file.path(pathFrom, "Gmacsall.out", fsep = fsep),
        verbose = FALSE,
        DatFile = Dat,
        CtlFile = Ctl,
        GmacsFile = gmacsFiles,
        nyrRetro = gmacsFiles$N_Year_Retro
      )
  }
  # =============================================== #
  
  # Define the number of theta parameters
  
  # Default number of parameters
  Defntheta <- nsex + 3 + 2 * nsex + 3
  # Depends on the Initial unfished conditions
  if (InitializeUnfished == 1) # FISHEDEQN
    ntheta <- Defntheta + nfleet
  if (InitializeUnfished == 2) # FREEPARS
    ntheta <- Defntheta + n_p*nsex*nMaturity*nShell
  if (InitializeUnfished == 3) # FREEPARSSCALED
    ntheta <- Defntheta + (n_p*nsex*nMaturity*nShell-1);
  
  # Key parameter control
  # -------------------------------------------------------------------------
  # Number of key parameters control (without size-class deviations)
  # CtlOut[["ntheta"]] <- ifelse(nsex == 1, yes = 9, no = 12)
  CtlOut[["ntheta"]] <- ntheta
  
  # Is the file created based on the snow crab assessment control file
  if (FromAssessFile) {
    
    # Add 3 to ntheta because the assessment is 2 sex
    if(Dat$N_sexes == 2)
      ntheta <- ntheta +3 
    
    # Natural mortality per sex
    # logR0
    # logRini, to estimate if NOT initialized at unfished (n68)
    # logRbar, to estimate if NOT initialized at unfished      #1
    # recruitment expected value (males or combined)
    # recruitment scale (variance component) (males or combined)
    # recruitment expected value (females)
    # recruitment scale (variance component) (females)
    # ln(sigma_R)
    # steepness
    # recruitment autocorrelation
    # Deviation for size-class - RELATED TO THE INITIALIZATION OF UNFISHED CONDITION
      # mature males/ immature males / mature females/ immature females
    
    # Matrix of key parameters control
    # Init_val | Lower_Bd | Upper_Bd | Phase | Prior_type | p1 | p2
    # if (nsex == 1) {
    #   theta <- CtlOut[["theta_control"]][c(1, 3:7, 10:12),]
    # } else {
    #   theta <- CtlOut[["theta_control"]][c(1:12),]
    # }
    if (nsex == 1) {
      theta <- CtlOut[["theta_control"]][c(1, 3:7, 10:ntheta),]
    } else {
      theta <- CtlOut[["theta_control"]][c(1:ntheta),]
    }
    
    # Fill in the initial values and turn off estimation if
    # useEstParamVals == TRUE
    # if (useEstParamVals[1]) {
    #   if (nsex == 1) {
    #     theta[, "Init_val"] <-
    #       AllOut$Param$theta[c(1, 3:7, 10:12), "Estimate"]
    #   } else {
    #     theta[, "Init_val"] <- AllOut$Param$theta[1:12, "Estimate"]
    #   }
    # }
    if (useEstParamVals[1]) {
      if (nsex == 1) {
        theta[, "Init_val"] <-
          AllOut$Param$theta[c(1, 3:7, 10:ntheta), "Estimate"]
      } else {
        theta[, "Init_val"] <- AllOut$Param$theta[1:ntheta, "Estimate"]
      }
    }
    
    if (TurnOffest[1]) {
      theta <- theta %>%
        dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
                                        Phase < 0 ~ Phase))
    }
    CtlOut[["theta_control"]] <- theta
    # -------------------------------------------------------------------------
    
    # Selectivity parameter controls
    # -------------------------------------------------------------------------
    # Selectivity
    CtlOut[["slx_nsel_period_in"]] <-
      .an(c(1, 1)) # Number of selectivity time period per fleet
    CtlOut[["slx_bsex_in"]] <-
      .an(c(0, 0)) # Male only selectivity
    CtlOut[["slx_type_in"]] <-
      .an(c(2, 0)) # logistic selectivity
    CtlOut[["slx_include_in"]] <-
      .an(c(0, 0)) # Insertion of fleet in another
    CtlOut[["slx_extra_in"]] <-
      .an(c(0, 0)) # Extra parameter for each pattern
    
    # Retention
    CtlOut[["ret_nret_period_in"]] <-
      .an(c(1, 1)) # Number of retention time period per fleet
    CtlOut[["ret_bsex_in"]] <-
      .an(c(0, 0)) # Male only retention
    CtlOut[["ret_type_in"]] <-
      .an(c(2, 6)) # Males & females retention types
    CtlOut[["slx_nret"]] <-
      .an(c(1, 0)) # boolean for retention/discard
    CtlOut[["ret_extra_in"]] <-
      .an(c(0, 0)) # Extra parameter for each pattern
    CtlOut[["slx_max_at_1_in"]] <-
      .an(c(1, 0)) # Selectivity for the maximum size class if forced to be 1?
    
    # Vulnerability controls
    # ==================================================== #
    
    # Selectivity parameters control
    # Fishing fleet
    if (nsex == 1) {
      tmpSelF <- CtlOut[["Selex_control"]][1:2, ]
    } else {
      tmpSelF <- CtlOut[["Selex_control"]][1:4, ]
    }
    
    # For survey - we are based on the selectivity of the BSFRF
    # Gear 5 in the snow crab assessment
    tmpSelSurv <- CtlOut[["Selex_control"]] %>%
      dplyr::filter(Fleet == 5)
    Ind1 <- tmpSelSurv[1, "Index"]
    
    if (nsex == 1) {
      Ind2 <- Ind1 + Dat$N_sizeC - 1
    } else {
      Ind2 <- Ind1 + 2 * Dat$N_sizeC - 1
    }
    tmpSelSurv <- tmpSelSurv %>%
      dplyr::filter(Index %in% Ind1:Ind2)  %>%
      dplyr::mutate(Fleet = 2)
    
    tmpSel <- rbind(tmpSelF, tmpSelSurv)
    rm(tmpSelF, tmpSelSurv)
    
    # Get the estimated values
    if (useEstParamVals[4]) {
      estSel <- AllOut$Param$Vul %>%
        dplyr::filter(Parameter %in% paste0("log_slx_pars_", tmpSel$Index))
      
      tmpSel <- tmpSel %>%
        mutate(Init_val = exp(estSel$Estimate))
    }
    # Turn of the estimation
    if (TurnOffest[4]) {
      tmpSel <- tmpSel %>%
        dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
                                        Phase < 0 ~ Phase))
    }
    CtlOut[["Selex_control"]] <- tmpSel
    rm(tmpSel)
    
    # Retention parameters control
    # fishing fleet
    tmpRetF <- CtlOut[["Ret_control"]] %>%
      dplyr::filter(Fleet == -1)
    
    if (nsex == 1) {
      tmpRetF <- tmpRetF[1:2, ]
    } else {
      tmpRetF <- tmpRetF[1:3, ]
      # Because retention if flat for females (equals to 0) - 1 param
    }
    
    # Survey
    # Retention null - One parameter
    tmpRetSurv <- CtlOut[["Ret_control"]] %>%
      dplyr::filter(Fleet == -5)  %>%
      dplyr::mutate(Fleet = -2)
    
    tmpRet <- rbind(tmpRetF, tmpRetSurv)
    rm(tmpRetF, tmpRetSurv)
    
    # Get the estimated values
    if (useEstParamVals[5]) {
      estRet <- AllOut$Param$Vul %>%
        dplyr::filter(Parameter %in% paste0("log_slx_pars_", tmpRet$Index))
      
      tmpRet <- tmpRet %>%
        mutate(Init_val = exp(estRet$Estimate))
    }
    # Turn of the estimation
    if (TurnOffest[5]) {
      tmpRet <- tmpRet %>%
        dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
                                        Phase < 0 ~ Phase))
    }
    CtlOut[["Ret_control"]] <- tmpRet
    rm(tmpRet)
    
    # Reindexing the parameters
    
    CtlOut$Selex_control$Index <- 1:nrow(CtlOut$Selex_control)
    CtlOut$Ret_control$Index <- (nrow(CtlOut$Selex_control) + 1):(nrow(CtlOut$Selex_control) +
                                                                    nrow(CtlOut$Ret_control))
    # ==================================================== #
    
    # Number of asymptotic selectivity parameter
    CtlOut[["NumAsympRet"]] <- 0
    # Asymptotic parameter control
    # CtlOut[["AsympSel_control"]] <-
    
    # Environmental parameters and random walk selectivity definition
    nslx_envpars <-
      c(CtlOut[["Selex_control"]]$Env_Link, CtlOut[["Ret_control"]]$Env_Link)
    nslx_envpars <- length(which(nslx_envpars > 0))
    
    CtlOut[["nslx_envpars"]] <- nslx_envpars
    
    # Env parameters if(nslx_envpars > 0)
    # CtlOut[["SlxEnvPar"]] <-
      
    # Estimation phase for the deviation parameter (random walk)
    CtlOut[["devParPhase"]] <- -1
    # -------------------------------------------------------------------------
    
    
    # Priors for catchabilities of surveys
    # -------------------------------------------------------------------------
    
    # Catchability control param
    # Based on BSFRF_TRWL MALES (fixed to 1)
    CtlOut[["q_controls"]] <-  CtlOut[["q_controls"]][5, ]
    # -------------------------------------------------------------------------
    
    # Additional survey CV control
    # -------------------------------------------------------------------------
    # Additional CV control param
    CtlOut[["add_cv_controls"]] <- CtlOut[["add_cv_controls"]][1, ]
    
    if (N_SurveyDF == 1) {
      # Additional variance control for each survey (0: ignore; >0 use)
      CtlOut[["add_cv_links"]] <- 0
    } else if (N_SurveyDF > 1) {
      CtlOut[["add_cv_links"]] <- rep(0, N_SurveyDF)
    } else {
      CtlOut[["add_cv_links"]] <- ""
    }
    # -------------------------------------------------------------------------
    
    # Penalties for fishing mortality control
    # -------------------------------------------------------------------------
    CtlOut[["f_controls"]] <- CtlOut[["f_controls"]][c(1, 5), ]
    # -------------------------------------------------------------------------
    
    # Size composition data control
    # -------------------------------------------------------------------------
    # Size composition likelihood type
    CtlOut[["nAgeCompType"]] <- c(2,2,2, 2)
    # Auto tail compression
    CtlOut[["bTailCompression"]] <- c(0,0,0, 0)
    # Initial value for effective sample size
    CtlOut[["nvn_ival"]] <- c(1,1,1, 1)
    # Phase for effective sample size
    CtlOut[["nvn_phz"]] <- c(-4,-4,-4, -4)
    # Should data be aggregated
    CtlOut[["iCompAggregator"]] <- c(1,1,1, 2)
    # 2:Survey-like predictions; 1:catch-like predictions
    CtlOut[["lf_catch"]] <- c(1,1,1, 2)
    # Lambda for effect N
    CtlOut[["lf_lambda"]] <- c(1,1,1, 1)
    # Weight for likelihood
    CtlOut[["lf_emphasis"]] <- c(1,1,1, 1)
    # -------------------------------------------------------------------------
    
    # Natural mortality control
    # -------------------------------------------------------------------------
    # Type of M specification
    # Constant natural mortality
    CtlOut[["m_type"]] <- 0
    
    # Is female M relative to M male?
    # 0: No absolute value; 1 relative
    CtlOut[["MrelFem"]] <- NULL
    
    # Phase of estimation for M
    CtlOut[["Mdev_phz_def"]] <- -4
    
    # standard deviation in M deviations
    CtlOut[["m_stdev"]] <- 0
    
    # Number of nodes for cubic spline or number of step-changes for option 3
    if (nsex == 1) {
      CtlOut[["m_nNodes_sex"]] <- 0
    } else {
      CtlOut[["m_nNodes_sex"]] <- matrix(0, nrow = 2, ncol = 1)
    }
    CtlOut[["m_nodeyear_sex"]] <- ""
    
    # Number of breakpoints in M by size class
    CtlOut[["nSizeDevs"]] <- 0
    # Size positions of breakpoints in M by size class
    CtlOut[["m_size_nodeyear"]] <- NULL
    # Specific initial value for natural mortality deviations
    CtlOut[["Init_Mdev"]] <- 0
    # Natural mortality deviations controls
    CtlOut[["Mdev_controls"]] <- NULL
    # -------------------------------------------------------------------------
    
    # Tagging control
    # -------------------------------------------------------------------------
    # Emphasis (likelihood weight) on tagging
    CtlOut[["tag_emphasis"]] <- 1
    # -------------------------------------------------------------------------
    
    # Other (additional) controls
    # -------------------------------------------------------------------------
    CtlOut[["Term_molt"]] <- 1 # Consider terminal molting?
    CtlOut[["rdv_phz"]] <- 1 # Phase for recruitment estimation
    CtlOut[["rec_prop_phz"]] <-
      2 # Phase for recruitment sex-ratio estimation
    CtlOut[["init_sex_ratio"]] <-
      0.5 # Initial value for expected sex-ratio
    CtlOut[["rec_ini_phz"]] <-
      -3 # Phase for initial recruitment estimation
    CtlOut[["verbose"]] <-
      1 # Verbose flag (0: off; 1: on; 2: objective function; 3: diagnostics)
    CtlOut[["bInitializeUnfished"]] <-
      3 # Initial conditions (1: unfished, 2: steady-state, 3: free params, 4: free params revised)
    CtlOut[["spr_lambda"]] <-
      1 # Proportion of mature male biomass for SPR reference points
    CtlOut[["nSRR_flag"]] <-
      0 # Stock-Recruit-Relationship (0 = none, 1 = Beverton-Holt)
    CtlOut[["TurnOffPhase"]] <-
      10 # Maximum phase (stop the estimation after this phase)
    CtlOut[["StopAfterFnCall"]] <-
      -1# Maximum number of function calls
    CtlOut[["CalcRefPoints"]] <-
      1 # Calculate reference points (0:no, 1: yes)
    CtlOut[["BRP_rec_sexR"]] <-
      0 # Use years specified to computed average sex ratio in the calculation of average recruitment for reference points (0 = off -i.e. Rec based on End year, 1 = on)
    CtlOut[["NyrEquil"]] <- 200 # Years to compute equilibria
    # -------------------------------------------------------------------------
    
    # Emphasis factor (weights for likelihood) controls
    # -------------------------------------------------------------------------
    # # Weights on catches for the likelihood component
    # CtlOut[["catch_emphasis"]] <- rep(1,3)
    # Penalties on deviations
    CtlOut[["Penalty_fdevs"]] <-
      matrix(
        c(1, 1, 0, 0, 0, 0, 0, 0),
        nrow = 2,
        ncol = 4,
        byrow = TRUE
      )
    
    # priors
    namEmph <- c(
      "Log_fdevs",
      "meanF",
      "Mdevs",
      "Rec_devs",
      "Initial_devs",
      "Fst_dif_dev",
      "Mean_sex-Ratio",
      "Molt_prob",
      "Free_selectivity",
      "Init_n_at_len",
      "Fdevs",
      "Fdovs",
      "Vul_devs"
    )
    
    for (n in namEmph) {
      if (n != "Log_fdevs")
        eval(parse(text = paste0(
          "CtlOut[['Penalty_emphasis']]$'", n, "' <- 0"
        )))
    }
    # -------------------------------------------------------------------------
  }
  
  # Update parameter values with the set up from the user input
  # Key parameter control
  # -------------------------------------------------------------------------
  # Fill in the initial values and turn off estimation if
  # useEstParamVals == TRUE
  if (!FromAssessFile) {
    if (useEstParamVals[1]) {
      CtlOut[["theta_control"]][, "Init_val"] <-
        AllOut$Param$theta[, "Estimate"]
    }
    if (TurnOffest[1]) {
      CtlOut[["theta_control"]] <- CtlOut[["theta_control"]] %>%
        dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
                                        Phase < 0 ~ Phase))
    }
  }
  # If parameters are define in the main script, then get those values
  # -------------------------------------------------------------------------
  
  # 1. Natural mortality ----
  # ================================ #
  # => GMACS: set a parameter for MATURE males and Females and use offset for
  # immature individuals i.e, set Maturity specific M
  CtlOut$theta_control[1, "Init_val"] <-
    mat_male_M        # M (mature male)
  if (nsex == 2) {
    CtlOut$theta_control[2, "Init_val"] <-
      mat_fem_M       # M (mature female)
  }
  # ================================
  
  # 2. Length-weight relationship type ----
  # ================================ #
  CtlOut[["lw_type"]] <- lw_type
  
  if (CtlOut[["lw_type"]] < 1 || CtlOut[["lw_type"]] > 3)
    stop("length-weight type can only be 1,2 or 3")
  
  if (CtlOut[["lw_type"]] == 1) {
    CtlOut[["lw_alfa"]] <-
      base::switch (.ac(nsex),
                    # alpha Length-weight relationship
                    "1" = weight_a_u,
                    "2" = c(weight_a_m, weight_a_f))
    CtlOut[["lw_beta"]] <-
      base::switch (.ac(nsex),
                    # beta Length-weight relationship
                    "1" = weight_b_u,
                    "2" = c(weight_b_m, weight_b_f))
    CtlOut[["mean_wt_in"]] <- ""
  } else {
    CtlOut[["lw_alfa"]] <- ""
    CtlOut[["lw_beta"]] <- ""
    
    if (nsex == 1) {
      # Input weight at size
      CtlOut[["mean_wt_in"]] <- wt_at_len_U
    } else {
      CtlOut[["mean_wt_in"]] <- wt_at_len
    }
  }

  # Fecundity for MMB/MMA calculation
  CtlOut[["maturity"]] <- CtlOut[["maturity"]][1:nsex,]
  
  # Proportion of mature at size by sex
  CtlOut[["legal_maturity"]] <- CtlOut[["legal_maturity"]][1:nsex,]
  # ================================
  
  
  # 3. Growth parameter controls ----
  # ================================ #
  # Options for the growth matrix
  CtlOut[["bUseCustomGrowthMatrix"]] <- CustomGrowthMatrix
  if (CtlOut[["bUseCustomGrowthMatrix"]] < 1 ||
      CtlOut[["bUseCustomGrowthMatrix"]] > 8)
    stop("Growth matrix type can only be 1-8")
  
  # Options for the growth increment model
  CtlOut[["bUseGrowthIncrementModel"]] <- GrowthIncrementModel
  
  # molt probability function
  CtlOut[["bUseCustomMoltProbability"]] <- CustomMoltProbability
  
  # Maximum sizes-class for recruitment
  rec_sizes <- rep(rec_sizes, nsex)
  if (nsex == 2 & length(rec_sizes) == 1) {
    cat("The maximum sizes-class for recruitment is set as the same for both sexes. \n")
  }
  CtlOut[["nSizeClassRec"]] <- rec_sizes
  
  # Number of size increment periods
  CtlOut[["nSizeIncVaries"]] <- nSizeIncBlock
  # Year(s) with changes in growth matrix - size increment (blank if no change)
  tmpSizeIncVaries <- NULL
  for (s in 1:nsex)
    tmpSizeIncVaries <-
    c(tmpSizeIncVaries, CtlOut[["nSizeIncVaries"]][s] - 1)
  
  if (max(tmpSizeIncVaries) > 0) {
    CtlOut[["iYrsSizeIncChanges"]] <- SizeIncBlock_Yrs
  } else {
    CtlOut[["iYrsSizeIncChanges"]] <- ""
  }
  
  # Number of molt periods
  CtlOut[["nMoltVaries"]] <- nMoltBlock
  # Year(s) molt period changes (blank if no change)
  tmpMoltVaries <- NULL
  for (s in 1:nsex)
    tmpMoltVaries <-
    c(tmpMoltVaries, CtlOut[["nMoltVaries"]][s] - 1)
  
  if (max(tmpMoltVaries) > 0) {
    CtlOut[["iYrsMoltChanges"]] <- MoltBlock_Yrs
  } else {
    CtlOut[["iYrsMoltChanges"]] <- ""
  }
  
  # Beta parameters are relative (0:No; 1:Yes)
  CtlOut[["BetaParRelative"]] <- 0
  
  # Model parameter control
  nGrwth <- 0
  nSizeIncPar <- 0
  for (s in 1:nsex)
  {
    if (CtlOut[["bUseGrowthIncrementModel"]] == 1)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * 3 # LINEAR_GROWTHMODEL
    if (CtlOut[["bUseGrowthIncrementModel"]] == 2)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * (nclass + 1) # INDIVIDUAL_GROWTHMODEL1
    if (CtlOut[["bUseGrowthIncrementModel"]] == 3)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * (nclass + 1) # INDIVIDUAL_GROWTHMODEL2
    if (CtlOut[["bUseGrowthIncrementModel"]] == 5)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * 3 # GROWTH_VARYK
    if (CtlOut[["bUseGrowthIncrementModel"]] == 6)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * 3 # GROWTH_VARYLINF
    if (CtlOut[["bUseGrowthIncrementModel"]] == 7)
      nGrwth <-
        nGrwth + CtlOut[["nSizeIncVaries"]][s] * 4 # GROWTH_VARYKLINF
    if (CtlOut[["bUseCustomMoltProbability"]] == 2)
      nSizeIncPar <-
        nSizeIncPar + CtlOut[["nMoltVaries"]][s] * 2 # LOGISTIC_PROB_MOLT
    if (CtlOut[["bUseCustomMoltProbability"]] == 3)
      nSizeIncPar <-
        nSizeIncPar + CtlOut[["nMoltVaries"]][s] * (nclass) # FREE_PROB_MOLT
  }
  # nGrwth <- nGrwth + nSizeIncPar
  CtlOut[["nGrwth"]] <- nGrwth
  CtlOut[["nSizeIncPar"]] <- nSizeIncPar
  
  # Growth Increment model parameters control
  if (GrowthIncrementModel == 1 || GrowthIncrementModel == 2) {
    if (GrowthIncrementModel == 1) {
      cat(
        "\t!! Considering a mean value between mature and immature individuals for the
growth increment parameters.\n"
      )
      
      # Fill in the parameter values if specified
      if (nsex == 1) {
        CtlOut[["Grwth_control"]] <- CtlOut[["Grwth_control"]][1:3,]
      }
      # Males
      if (!is.null(alpha_grow_m_imm) && !is.null(alpha_grow_m_mat))
        CtlOut[["Grwth_control"]][1, "Init_val"] <-
          mean(alpha_grow_m_imm, alpha_grow_m_mat)
      if (!is.null(beta_grow_m_imm) && !is.null(beta_grow_m_mat))
        CtlOut[["Grwth_control"]][2, "Init_val"] <-
          mean(beta_grow_m_imm, beta_grow_m_mat)
      if (nsex == 2) {
        # Females
        if (!is.null(alpha_grow_f_imm) && !is.null(alpha_grow_f_mat))
          CtlOut[["Grwth_control"]][4, "Init_val"] <-
            mean(alpha_grow_f_imm, alpha_grow_f_mat)
        if (!is.null(beta_grow_f_imm) && !is.null(beta_grow_f_mat))
          CtlOut[["Grwth_control"]][5, "Init_val"] <-
            mean(beta_grow_f_imm, beta_grow_f_mat)
      }
      
      # Use the estimated values
      if (useEstParamVals[2]) {
        if (nsex == 1) {
          CtlOut[["Grwth_control"]] <- CtlOut[["Grwth_control"]][1:3,] %>%
            mutate(Init_val = AllOut$Param$Grwth[1:3, "Estimate"])
        } else {
          CtlOut[["Grwth_control"]] <- CtlOut[["Grwth_control"]][1:6,] %>%
            mutate(Init_val = AllOut$Param$Grwth[1:6, "Estimate"])
        }
      }
      # Turn of the estimation
      if (TurnOffest[2]) {
        CtlOut[["Grwth_control"]] <- CtlOut[["Grwth_control"]] %>%
          dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
                                          Phase < 0 ~ Phase))
      }
    }
  }
  
  # Molting probability parameters
  if (CustomMoltProbability == 2) {
    # Logisitc parameters
    if (nsex == 1) {
      CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]][1:2, ]
    } else {
      CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]][1:4, ]
    }
    # Use the estimated values
    if (useEstParamVals[3]) {
      if (nsex == 1) {
        # CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]][1:2, ] %>%
        #   mutate(Init_val = AllOut$Param$Grwth[3:4, "Estimate"])
      } else {
        # CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]][1:4, ] %>%
        #   mutate(Init_val = AllOut$Param$Grwth[7:12, "Estimate"])
      }
    }
    # Turn of the estimation
    if (TurnOffest[3]) {
      # CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]] %>%
      #   dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
      #                                   Phase < 0 ~ Phase))
    }
  } else if (CustomMoltProbability == 3) {
    # free parameters
    
    # Use the estimated values
    if (useEstParamVals[3]) {
      if (nsex == 1) {
        # CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]][1:2, ] %>%
        #   mutate(Init_val = AllOut$Param$Grwth[1:2, "Estimate"])
      } else {
        CtlOut[["MoltProb_control"]] <-
          CtlOut[["MoltProb_control"]][1:4,] %>%
          mutate(Init_val = AllOut$Param$Grwth[1:4, "Estimate"])
      }
    }
    # Turn of the estimation
    if (TurnOffest[3]) {
      # CtlOut[["MoltProb_control"]] <- CtlOut[["MoltProb_control"]] %>%
      #   dplyr::mutate(Phase = case_when(Phase > 0 ~ -1 * Phase,
      #                                   Phase < 0 ~ Phase))
    }
  } else {
    CtlOut[["MoltProb_control"]] <- NULL
  }
  
  # Custom growth-increment matrix or size-transition matrix
  # GROWTH_FIXEDGROWTHTRANS || GROWTH_FIXEDSIZETRANS
  if (CtlOut[["bUseCustomGrowthMatrix"]] == 1 ||
      CtlOut[["bUseCustomGrowthMatrix"]] == 2) {
    maxSizeIncVaries = max(CtlOut[["nSizeIncVaries"]])

    # GROWTH_FIXEDSIZETRANS
    if (CtlOut[["bUseCustomGrowthMatrix"]] == 2)
      cat("\t!! Considering the male size transition matrix for all individuals.\n")
    CtlOut[["CustomGrowthMatrix"]] <- matrix(-1, 
                                             nrow = Dat$N_sizeC * nsex * maxSizeIncVaries,
                                             ncol = Dat$N_sizeC)
  } else {
  CtlOut[["CustomGrowthMatrix"]] <- NULL
  }
  
  # Custom molt probability matrix
  if (CtlOut[["bUseCustomMoltProbability"]] == 0) {
    # Fixed Molt probability
    # CtlOut[["CustomMoltProbabilityMatrix"]] <-
  } else {
    CtlOut[["CustomMoltProbabilityMatrix"]] <- NULL
  }
  # ================================
  
  
  # Immature/mature natural mortality ----
  # ================================ #
  
  # maturity specific natural mortality?
  CtlOut[["m_maturity"]] <- 1
  # (1:yes; 0:no - only for use if nmature > 1)
  
  if (nsex == 1) {
    CtlOut[["m_mat_controls"]] <- CtlOut[["m_mat_controls"]][1, ]
    # Fill in with the parameters
    # Natural mortality Immature Males
    CtlOut[["m_mat_controls"]][1, "Init_val"] <-
      log(imm_male_M / mat_male_M)
  } else {
    CtlOut[["m_mat_controls"]] <- CtlOut[["m_mat_controls"]][1:2, ]
    # Fill in with the parameters
    CtlOut[["m_mat_controls"]][1, "Init_val"] <-
      log(imm_male_M / mat_male_M)
    # Natural mortality Immature Females
    CtlOut[["m_mat_controls"]][2, "Init_val"] <-
      log(imm_fem_M / mat_fem_M)
  }
  # ================================
  # -------------------------------------------------------------------------

  # Write the new GMACS control input file
  if(pathFrom == file.path(dir_GMACS, "model", "build",
                           "assessment_files", fsep = fsep)){
    overwrite <- FALSE
  }
  
  if(overwrite){
    Dir <- pathFrom
    FileName <- filename
  } else {
    Dir <- pathTo
    if(!is.null(filenameTo)){
      FileName <- filenameTo
    } else {
      FileName <- filename
    }
  }
  
  if (FromAssessFile){
    # Read the GMACS.dat file to get file names 
    gmacsFiles <-
      readGMACS.dat(path = file.path(pathTo, "gmacs.dat", fsep = fsep),
                    verbose = FALSE)
    Dat <-
      readGMACSdat(
        FileName = file.path(pathTo, gmacsFiles$DatFileName, fsep = fsep),
        verbose = FALSE
      )
  }
  
  if(verbose)
    cat("=> Writing the input control file in the following directory:\n", Dir, "\n")
  writeGmacsctlfile(Dir = Dir, 
                    FileName = FileName,
                    DatFile = Dat,
                    CtlFile = CtlOut,
                    stock = "Snow crab", 
                    model_name = model_name, 
                    Ass_Year = "")
}

#' Function to write the base input data file for a MSE
#'
#' @param pathfrom (character string)- Path to the directory of the input data
#' file to be modified.
#' @param pathTo (character string)- path to the directory where to save the new
#' input control file.
#' @param filename (character string)- Name of the input data file to be modified
#' @param model_name (character string)- Name for the model/MSE. This will be use
#' to create the folder (if it does not already exist) where GMACS outputs will be
#' housed.
#' @param overwrite (logical)- Do you want to overwrite the \code{filename} file.
#' @param filenameTo (character string)- Name of the new input data file to write.
#' @param verbose (logical)- Print statement.
#'
#' This function write the input data file for a given configuration of a MSE
set_GMACSprojFile <- function(pathFrom = NULL,
                              pathTo = NULL,
                              model_name = NULL,
                              filename = NULL,
                              overwrite = NULL,
                              filenameTo = NULL,
                              verbose = NULL) {
  # Local declarations
  fsep <- .Platform$file.sep
  
  # Check Input
  if (pathFrom == pathTo) {
    if (!overwrite && (is.null(filenameTo) || filenameTo == filename)) {
      message(
        "=> Something is wrong with the set up to write the data input file.\n
You can either:
\t1) provide a different name for the 'pathTo' directory
\t2) allow overwriting the input data file in the 'pathFrom' directory, or
\t3) privide a name (or a different one if already provided) for 'filenameTo'."
      )
      stop()
    }
    if (overwrite) {
      cat("=> The input data file in the 'pathFrom' directory will be overwritten !\n")
    }
  }
  if (is.null(pathTo) && !overwrite) {
    message(
      "=> Please provide the 'pathTo' directory or
allow overwriting the input data file in the 'pathFrom' directory."
    )
    stop()
  }
  if (is.null(pathTo) && overwrite) {
    cat("=> The input data file in the 'pathFrom' directory will be overwritten !\n")
  }
  
  # Create pathTo if needed
  if (!is.null(pathTo) & !dir.exists(pathTo)) {
    if (verbose)
      cat("=> Creating the folowing directory:\n", pathTo, "\n")
    dir.create(path = pathTo)
  }
  
  # Read a given data file to be modified
  Prj <- readGMACSprj(FileName = file.path(pathFrom, filename),
                      verbose = verbose)
  PrjOut <- Prj
  rm(Prj)
  
  # Check the F applied to each fleet
  # 0 to apply F35%; 1 if F is to be fixed
  # Here we apply F35% to the fishery and F to the survey
  PrjOut[["Ffixed"]] <- c(0, 1)
  
  
  # Check settings for the First years in the projection 
  # ======================================================================= #
  # -------------------------------------------------------------------------
  WhichMod <- names(PrjOut[grepl(pattern = "_syr", x = names(PrjOut))])
  for(w in WhichMod){
    if(PrjOut[[w]] < StartY_Ass && PrjOut[[w]]!=0)
      eval(parse(text = paste0(
        "PrjOut[['",w,"']] <- StartY_Ass + 2"
      )))
  }
  
  # ------------------------------------------------------------------------------
  
  # Write the new GMACS projection input file
  if(pathFrom == file.path(dir_GMACS, "model", "build",
                           "assessment_files", fsep = fsep))
    overwrite <- FALSE
  
  if (overwrite && is.null(pathTo)) {
    Dir <- pathFrom
    FileName <- filename
  } else {
    Dir <- pathTo
    if (!is.null(filenameTo)) {
      FileName <- filenameTo
    } else {
      FileName <- filename
    }
  }
  if (verbose)
    cat("=> Writing the input projection file in the following directory:\n",
        Dir,
        "\n")
  writeGmacsprjfile(
    Dir = Dir,
    FileName = FileName,
    PrjFile = PrjOut,
    stock = "Snow crab",
    model_name = model_name,
    Ass_Year = ""
  )
}
