## Demographic settings
#----------------------
size_class_settings <- "fine"
# either "fine" (5 cm per 5 cm like GMACS)
# or "rough" (4 size classes accordingly to the IPM)
# In spatial IPM, size classes are: 0 < size1 =<40 / 40 < size 2 =< 78 / 78<size3 =<101 / 101<size4

# Original settings
if(size_class_settings == "fine"){
  binsize  <- 5
  sizes		 <-seq(27.5,132.5,binsize)
  binclass <- c(sizes - binsize / 2, sizes[length(sizes)] + binsize / 2)
  n_p <- as.numeric(length(sizes))
}

if(size_class_settings == "rough"){
  binclass <- c(0,40,78,101,132.5)
  sizes		 <- (binclass[1:(length(binclass) - 1)] + binclass[2:length(binclass)]) / 2
  n_p <- as.numeric(length(sizes))
}

ad_size = 45 # size at which crabs are considered mature

## Sex and maturity
#------------------
sexN = 2

imm_fem_M = 0.32
imm_male_M = 0.32
mat_fem_M = 0.26
mat_male_M = 0.28

morta_sd_imm = 0.1
morta_sd_mat = 0.1

morta_model = "cody_model"
if(morta_model == "max_model"){
  
  #==Maxime's morta parameterization
  ## Load morta settings for spatially varying parameters
  G_spatial <- F
  Morta_spatial <- F # are life history parameters spatial?
  source("4_full_MSE/LHP/morta_settings.R")
  
  ## Load function for morta
  source("4_full_MSE/LHP/morta.R")
  
}


## Recruitment
#-------------
rec_sizes    <-5
prop_rec     <-c(.25,.4,.25,.075,0.025)

if(size_class_settings == "rough"){
  rec_sizes    <-1
  prop_rec     <-1
}

if(size_class_settings == "fine"){
  rec_sizes    <-5
  prop_rec     <-c(.25,.4,.25,.075,0.025)
}
