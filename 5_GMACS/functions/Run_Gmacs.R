runGMACS <- function(){
  
  # Check for the existence of the .exe
  dir_model <- file.path(dir_GMACS, "model", fsep = fsep)
  
  gmacs_exe <- ifelse(isWindowsOS(),"gmacs.exe","gmacs")
  ADMBpaths <- ifelse(.Platform$OS.type == "windows",
                      "ADpaths_Windows.txt",
                      "ADpaths_MacOS.txt")
  
  if(!file.exists(file.path(dir_model, gmacs_exe, fsep = fsep))){
    ADMBpaths <- file.path(dir_GMACS, ADMBpaths, fsep = fsep)
    createGmacsExe(
      vv = 1,
      Dir = dir_model,
      verbose = FALSE,
      ADMBpaths = ADMBpaths
    )
  }
  
  # Run GMACS
  GMACS(
    Spc = model_name,
    GMACS_version = "MSE",
    Dir = dir_model,
    ASS = FALSE,
    compile = 0,
    run = TRUE,
    LastAssDat = FALSE,
    ADMBpaths = ADMBpaths,
    make.comp = FALSE,
    verbose = TRUE,
    cleanOut = FALSE
  )

}


