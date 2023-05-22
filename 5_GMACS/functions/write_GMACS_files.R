# Script to write the data and control files
# when GMACS is used as the assessment model in the MSE

# This script is sourced before the loop and set up the base
# input file. Catch, survey and other data are updated in the loop
# by sourcing the "update_GMACS.R" script

# 1. Save the current environment
var.to.save <- ls()

# Local declarations ----
fsep <- .Platform$file.sep
dir_GMACS <- file.path(here::here(), "5_GMACS", fsep = fsep)

# 2. Source the user Input for GMACS
# source(file = file.path(here::here(), 
#                         "5_GMACS", "scripts_mse", 
#                         paste0(Gmacs_Input_file, ".R")))
source(file = file.path(dir_GMACS, "scripts_mse", 
                        paste0(Gmacs_Input_file, ".R"), fsep = fsep))


# 3. load the functions
source(file.path(
  dir_GMACS,
  "functions",
  "Get_GMACS_dat_ctl_files.R",
  fsep = fsep
))

# 3. Get the gmacs.dat and proj files from the 
# pathFrom directory if it is different from the pathTo one

pathTo <- file.path(dir_GMACS, "model", "build",
                    model_name, fsep = fsep)

# Create pathTo if needed
if (!is.null(pathTo) & !dir.exists(pathTo)) {
  if (verbose)
    cat("=> Creating the folowing directory:\n", pathTo, "\n")
  dir.create(path = pathTo)
}

# 4. Read in the gmacs.dat to get the names of the data and control files
gmacsFiles <-
  readGMACS.dat(path = file.path(pathFrom, "gmacs.dat", fsep = fsep),
                verbose = FALSE)

# Copy files
if(pathFrom != pathTo){
  files <- file.path(pathFrom, c("gmacs.dat", gmacsFiles$PrjFileName), fsep = fsep)
  file.copy(from = files, 
            to = file.path(pathTo, c("gmacs.dat", gmacsFiles$PrjFileName), fsep = fsep),
            overwrite = TRUE,
            copy.date = TRUE)
}

# 5. Write the new data file
set_GMACSdataFile(
  pathFrom = pathFrom,
  pathTo = pathTo,
  model_name = model_name,
  filename = gmacsFiles$DatFileName,
  overwrite = overwrite,
  filenameTo = datfilenameTo,
  verbose = verbose
)

# 6. Write the new control file
set_GMACSctlFile(
  pathFrom = pathFrom,
  pathTo = pathTo,
  model_name = model_name,
  filename = gmacsFiles$CtlFileName,
  overwrite = overwrite,
  filenameTo = ctlfilenameTo,
  verbose = verbose
)
# pathFrom = pathFrom
# pathTo = pathTo
# model_name = model_name
# filename = gmacsFiles$CtlFileName
# overwrite = overwrite
# filenameTo = ctlfilenameTo
# verbose = verbose

# 7. Write the new projection file
set_GMACSprojFile(
  pathFrom = pathFrom,
  pathTo = pathTo,
  model_name = model_name,
  filename = gmacsFiles$PrjFileName,
  overwrite = overwrite,
  filenameTo = "snow.prj",
  verbose = verbose
)
# pathFrom = pathFrom
# pathTo = pathTo
# model_name = model_name
# filename = gmacsFiles$PrjFileName
# overwrite = overwrite
# filenameTo = "snow.prj"
# verbose = verbose

# Restore the environment as it was but saved some stuff needed in the loop
var.to.save <- c(var.to.save,
                 'dir_GMACS', 'model_name',
                 'Split_CatchDF',
                 'Nsam_SizeC')
rm(list = setdiff(ls(), var.to.save))
var.to.save <- ls()

# ===============
