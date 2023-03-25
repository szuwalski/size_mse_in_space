# Snow crab GMACS model output
# This is the gmacs.rep file

rm(list=ls())

library(magrittr)
library(dplyr)

# Directories
load(file = "4_full_MSE/data/Snow_GMACS_repfile.Rdata")

# Access names of the repfile
names(repfile)

# 1. Abundances
N_total <- repfile$N_total
N_females <- repfile$N_females
N_females_new <- repfile$N_females_new
N_females_old <- repfile$N_females_old
N_females_mature <- repfile$N_females_mature
N_males <- repfile$N_males
N_males_new <- repfile$N_males_new
N_males_old <- repfile$N_males_old
N_males_mature <- repfile$N_males_mature

# 2. Recruitment
rec_sdd <- repfile$rec_sdd
rec_dev <- repfile$rec_dev
recruits <- repfile$recruits

# 3. Selectivities 

# Let's take a look at the fleets
repfile$fleetname
  # Interest for 'Pot_Fishery' and 'NMFS_Trawl_1989' (probably - longest data serie)
# Capture
slx_capture <- repfile$slx_capture %>% 
  dplyr::filter(fleet == 1)
# Retained 
slx_retained <- repfile$slx_retaind %>% 
  dplyr::filter(fleet == 1)

# 4. Molt probability
molt_prob <- repfile$molt_probability

# 5. Molt increment
molt_Inc <- repfile$molt_increment

# 6. Growth
growth_transition <- repfile$growth_matrix
growth_param <- repfile$Growth_param  %>% 
  dplyr::mutate(model = "")

growth_param[1:6,] <- growth_param[1:6,] %>% 
  dplyr::mutate(Parameter = c(paste("Male", c("alpha", "beta", "scale"), sep = "_"),
                paste("Female", c("alpha", "beta", "scale"), sep = "_")))  %>% 
  dplyr::mutate(model = "Increment")

growth_param[7:dim(growth_param)[1],] <- growth_param[7:dim(growth_param)[1],] %>% 
  dplyr::mutate(Parameter = c(paste(rep(c("Male_SC", "female_SC"), each = length(repfile$mid_points)),
      rep(1:length(repfile$mid_points), 2), sep="_")))  %>% 
  dplyr::mutate(model = "molt prob")


# Selectivity
selec_curve = repfile$selectivity$Start_Y %>% 
  filter(fleet == 1) %>% 
  dplyr::select_at(vars(starts_with("SizeC_")))
