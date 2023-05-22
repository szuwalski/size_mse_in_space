## Demographic settings
#----------------------

##------------------------------------ Size class settings ---------------------------------------------

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


##------------------------------------ Mortality ---------------------------------------------

imm_fem_M = 0.32
imm_male_M = 0.32
mat_fem_M = 0.26
mat_male_M = 0.28
morta_sd_imm = 0.1
morta_sd_mat = 0.1
morta_model = "coldpool" # "cody_model", "max_model", "coldpool"
if(morta_model == "max_model"){
  
  #==Maxime's morta parameterization
  ## Load morta settings for spatially varying parameters
  G_spatial <- F
  Morta_spatial <- F # are life history parameters spatial?
  source("4_full_MSE/LHP/morta_settings.R")
  
  ## Load function for morta
  source("4_full_MSE/LHP/morta.R")
  
}else if(morta_model == "coldpool"){
  
  mod<-gam(data=m_gam_dat,Imm_mort~s(coldpool,k=3),family=tw())
  # summary(mod)
  # plot(mod,pages=1)
  
}


##------------------------------------ Recruitment ---------------------------------------------

rec_sizes <- 5
prop_rec <- c(.25,.4,.25,.075,0.025)

if(size_class_settings == "rough"){
  rec_sizes <- 1
  prop_rec <-1
}

if(size_class_settings == "fine"){
  rec_sizes <-5
  prop_rec <-c(.25,.4,.25,.075,0.025)
}

##------------------------------------ Growth ---------------------------------------------

# Parameterize with GMACS
growth_transition <- repfile$growth_matrix
growth_param <- repfile$Growth_param %>% 
  dplyr::mutate(model = "")

growth_param[1:6,] <- growth_param[1:6,] %>% 
  dplyr::mutate(Parameter = c(paste("Male", c("alpha", "beta", "scale"), sep = "_"),
                              paste("Female", c("alpha", "beta", "scale"), sep = "_"))) %>% 
  dplyr::mutate(model = "Increment")

growth_param[7:dim(growth_param)[1],] <- growth_param[7:dim(growth_param)[1],] %>% 
  dplyr::mutate(Parameter = c(paste(rep(c("Male_SC", "female_SC"),
                                        each = length(repfile$mid_points)),
                                    rep(1:length(repfile$mid_points), 2), sep="_")))  %>% 
  dplyr::mutate(model = "molt prob")

# Type of model (Cody's model --> non-spatial life-history parameters, Maxime's model --> spatially varying life-history parameters)
growth_model <- "cody_model" # "max_model" "cody_model"

growth_param_est = growth_param[1:6,] %>% 
  dplyr::select(Estimate) %>% 
  unlist

#==Cody's growth parameterization for generation of non spatially varying growth parameters
#==GMACS outputs
alpha_grow_f_imm <- growth_param_est[4]
alpha_grow_m_imm <- growth_param_est[1]
beta_grow_f_imm <- - growth_param_est[5]
beta_grow_m_imm <- - growth_param_est[2]

alpha_grow_f_mat <- growth_param_est[4]
alpha_grow_m_mat <- growth_param_est[1]
beta_grow_f_mat <- - growth_param_est[5]
beta_grow_m_mat <- - growth_param_est[2]

scale_f_imm = 1 # growth_param_est[6]
scale_f_mat = 1 # growth_param_est[6]
scale_m_imm = 1 # growth_param_est[3]
scale_m_mat = 1 # growth_param_est[3]

growth_sd_imm<-c(5,4) # Where to find better values?
growth_sd_mat<-c(5,4)

#==plug temp into growth curve
f_postmolt_imm<-sizes + (alpha_grow_f_imm + beta_grow_f_imm*sizes)/scale_f_imm
m_postmolt_imm<-sizes + (alpha_grow_m_imm + beta_grow_m_imm*sizes)/scale_m_imm
f_postmolt_mat<-sizes + (alpha_grow_f_mat + beta_grow_f_mat*sizes)/scale_f_mat
m_postmolt_mat<-sizes + (alpha_grow_m_mat + beta_grow_m_mat*sizes)/scale_m_mat

if(growth_model == "max_model"){
  
  #==Maxime's growth parameterization
  ## Load growth settings for spatially varying parameters
  G_spatial <- F 
  Growth_spatial <- F # are life history parameters spatial?
  source("4_full_MSE/LHP/growth_settings.R")
  
  ## Load function for growth
  source("2_Max_spatial_projection/LHP_functions/growth.R")
  
  ## Dimensions note in the same order as Cody's codes --> fix?
  growth_m_imm <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
  growth_m_mat <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
  growth_f_imm <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
  growth_f_mat <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
  
}


##------------------------------------ Molting probability ---------------------------------------------

terminal_molt<-1
if(size_class_settings == "fine") term_molt_prob<-c(0,0,0,.1,.2,.3,.4,.4,.4,.4,.4,.75,.9,1,1,1,1,1,1,1,1,1)
if(size_class_settings == "rough") term_molt_prob<-c(0.01, 0.3,0.58,0.9)
# plot(term_molt_prob~sizes)

##------------------------------------ Recruitment ---------------------------------------------

recruit_female = unlist(Snow_Out$Recruitment$recruits[1,])
meanlog_recruit_female = mean(log(recruit_female))
sdlog_female = sd(log(recruit_female))

recruit_male = unlist(Snow_Out$Recruitment$recruits[2,])
meanlog_recruit_male = mean(log(recruit_male))
sdlog_male = sd(log(recruit_male))

# Only mean recruitment is affected by growth
sd_meanlog_recruit_female = recruit_female / 0.1
sd_meanlog_recruit_male = recruit_male / 0.1

recruit_model = "cody_model"
if(recruit_model == "max_model"){
  
  #==Maxime's recruit parameterization
  ## Load recruit settings for spatially varying parameters
  G_spatial <- F 
  Recruit_spatial <- F # are life history parameters spatial?
  source("4_full_MSE/LHP/recruit_settings.R")
  
  ## Load function for recruitment
  source("4_full_MSE/LHP/recruitment.R")
  
}

##------------------------------------ Fishery pars ---------------------------------------------

fish_50<-95
fish_95<-101
fish_sel<-1/(1+exp(-log(19)*(sizes-fish_50)/(fish_95-fish_50)))


##--------------------------- weight at length parameters ---------------------------------------

weight_a_f<-0.001
weight_b_f<-3

weight_a_m<-0.0012
weight_b_m<-3  

wt_at_len<-rbind(weight_a_f*sizes^weight_b_f,weight_a_m*sizes^weight_b_m)  


##------------------------------------ Movement -------------------------------------------------

compute_movement_matrix = F
if(compute_movement_matrix){
  
  # Compute adjacency matrix
  A_gg = dnearneigh(point_expand[,c("Var1","Var2")],d1=0.45,d2=0.65)
  A_gg_mat = nb2mat(A_gg)
  A_gg_mat[which(A_gg_mat > 0)] = 1
  # plot(A_gg_mat[1:100,1:100])
  
  ## Diffusion
  #-----------
  # D * DeltaT / A
  # with D: diffusion coefficient, here btwn 0.1 and 1.1 km(^2?) per day (Cf. Olmos et al. SM)
  # DeltaT: time interval btwn time steps, 
  # A: area of grid cells
  # par(mfrow = c(3,2))
  D = 0.5
  DeltaT = (1/30) # convert day in month
  A = mean(cell_area)
  diffusion_coefficient = D * DeltaT / A
  diffusion_coefficient = 2 ^ 2
  diffusion_gg = A_gg_mat * diffusion_coefficient
  diag(diffusion_gg) = -1 * colSums(diffusion_gg)
  
  ## Taxis for juveniles
  #---------------------
  taxis_coef_juv = 10^2 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations
  preference_g_juv = (init_juv - mean(init_juv)) / sd(init_juv)
  # plot(t(preference_g_juv), breaks=20)
  preference_g_juv = as.vector(t(preference_g_juv))
  # # check
  # test = matrix(preference_g,nrow = 40,ncol = 40,byrow = T)
  # plot(test)
  
  taxis_gg_juv = A_gg_mat *  taxis_coef_juv * DeltaT / A * exp(outer(preference_g_juv, preference_g_juv, "-") * DeltaT / sqrt(A))
  diag(taxis_gg_juv) = -1 * colSums(taxis_gg_juv)
  
  # Total
  mrate_gg_juv = taxis_gg_juv # + diffusion_gg
  if( any(mrate_gg_juv-diag(diag(mrate_gg_juv))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
  mfraction_gg_juv = Matrix::expm(mrate_gg_juv)
  stationary_g_juv = eigen(mfraction_gg_juv)$vectors[,1]
  stationary_g_juv = stationary_g_juv / sum(stationary_g_juv)
  # test = matrix(stationary_g_juv,nrow=40,ncol=40,byrow = T)
  # plot(test)
  
  ## Taxis for adults
  #------------------
  taxis_coef_ad = 10^2 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations
  preference_g_ad = (init_adult - mean(init_adult)) / sd(init_adult)
  preference_g_ad = as.vector(t(preference_g_ad))
  # # Check
  # test = matrix(preference_g,nrow = 40,ncol = 40,byrow = T)
  # plot(test)
  
  taxis_gg_ad = A_gg_mat *  taxis_coef_ad * DeltaT / A * exp(outer(preference_g_ad, preference_g_ad, "-") * DeltaT / sqrt(A))
  diag(taxis_gg_ad) = -1 * colSums(taxis_gg_ad)
  
  # Total
  mrate_gg_ad = taxis_gg_ad # + diffusion_gg
  if( any(mrate_gg_ad-diag(diag(mrate_gg_ad))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
  
  mfraction_gg_ad = Matrix::expm(mrate_gg_ad)
  stationary_g_ad = eigen(mfraction_gg_ad)$vectors[,1]
  stationary_g_ad = stationary_g_ad / sum(stationary_g_ad)
  # test = matrix(stationary_g_ad,nrow=40,ncol=40,byrow = T)
  # x11();plot(test)
  
  save(data=mfraction_gg_ad,file="4_full_MSE/data/mfraction_gg_ad.RData")
  save(data=mfraction_gg_juv,file="4_full_MSE/data/mfraction_gg_juv.RData")
  
}else{
  
  load("4_full_MSE/data/mfraction_gg_ad.RData")
  load("4_full_MSE/data/mfraction_gg_juv.RData")
  
}

