rm(list=ls())

###########################################################################
# -------------------------------------------------------------------------
# STEP 1 SIMULATE LHP
# -------------------------------------------------------------------------
###########################################################################

setwd("C:/Users/Maxime/Documents/Git/size_mse_in_space")
source("spatial_projection/LHP_functions/libraries.R")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 0- Model Variable/parameters/integers------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# 0-a- Output file --------------------------------------------------------
# -------------------------------------------------------------------------
DateFile = paste0(getwd(), "/spatial_projection/Outputs/")
dir.create(DateFile)  

# 0-b- Spatial resolution -------------------------------------------------
# -------------------------------------------------------------------------
source("spatial_projection/LHP_functions/spatial_grid.R")
#-- Define grid to predict the Field
lat		   <- round(seq(70,51.5,length.out=40),2)
lon <-  round(seq(-179,-155,length.out=40),2)
n_lat <- length(lat)
n_lon <- length(lon)

spatial_grid <- spatial_grid(lon,lat)
attach(spatial_grid)

# 0-c-Size bins -----------------------------------------------------------
# -------------------------------------------------------------------------
binsize  <-5
sizes		 <-seq(27.5,132.5,binsize)
n_p = as.numeric(length(sizes))

# 0-d- Years --------------------------------------------------------------
# -------------------------------------------------------------------------
year_step <-12 # months
year_n	 <-30 # nber of years
proj_period	<-seq(1,year_step*year_n)
n_t <- length(proj_period)

#-- Years for climate scenarios
inits_year =2022
Years_climsc_temp <- ((rep(c(inits_year:(inits_year-1+year_n)),year_step)))
Years_climsc <- Years_climsc_temp[ order((rep(c(inits_year:(inits_year-1+year_n)),year_step)))]
month <- seq(1,12,1)

#-- Binary vectors related to period that determine when life events happen
#--July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun
survey_time   <-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),year_n)
fish_time		  <-rep(c(0,0,0,0,0,0,0,0,1,1,0,0),year_n)
recruit_time	<-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)
move_time		  <-rep(c(1,1,1,1,1,1,1,1,1,1,1,1),year_n)
molt_time_month_m <- c(0,0,0,0,0,0,0,0,0,1,0,0)
molt_time_month_f <- c(0,0,0,0,0,0,0,0,0,1,0,0)
molt_time_m <- rep(molt_time_month_m,year_n)
molt_time_f <- rep(molt_time_month_f,year_n)
molt_time  	  <-rbind(molt_time_m,molt_time_f) # females first, males second
mate_time 	  <-rep(c(0,0,0,0,0,0,0,1,0,1,0,0),year_n)

# 0-e- Sex ----------------------------------------------------------------
# -------------------------------------------------------------------------
sexN		 <- 2
n_n <- sexN


# 0-f- Climate scenarios --------------------------------------------------
# -------------------------------------------------------------------------
clim_sc =c("rcp45") # climate scenario to test
#n_cs <- length(clim_sc)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 1- Simulate Natural mortality -------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

imm_fem_M    <-0.32
imm_male_M   <-0.32
mat_fem_M    <-0.26
mat_male_M   <-0.28

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 2- Simulate Recruitment -------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# 2-a- Spatial  ----------------------------------------------------------
# -------------------------------------------------------------------------
R_spatial = FALSE

# 2-b Parameters ----------------------------------------------------------
rec_sizes    <- 5
prop_rec     <-c(.25,.4,.25,.075,0.025)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 3- Simulate growth ------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
LHP = "growth"

# 3-a- Spatial  ----------------------------------------------------------
# -------------------------------------------------------------------------

# - Is growth spatially defined ? -----------------------------------------
Growth_spatial = TRUE

# - Parameters for spatial field ------------------------------------------
scale_g <- 1.5 # scale of the RF
x_omega = 0.1 # how much of the total variance represent omega (in %)
x_epsilon = 0.05 # how much of the total variance represent epsilon (in %)


# 3-b- Parameters ---------------------------------------------------------
# -------------------------------------------------------------------------
# -- Growth : mean, sd and scale parameter for growth transition matrix
growth_par_beta <- array(0,c(2,3,4))
n_grpar_growth = dim(growth_par_beta)[2] # number of parameter to define growth 

# --- Growth is define by 2 parameters : initial and last size increment 
# ---- Female immature
growth_par_beta[1,1,1] <- 3 # initial size increment  
growth_par_beta[1,2,1] <- 15 #last size increment : mean
growth_par_beta[1,3,1] <-  0.3 #  scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,1] <- 0.1 # initial size increment  : sd
growth_par_beta[2,2,1] <- 0.1 #last size increment : sd
growth_par_beta[2,3,1] <-  NA 

# ---- Female mature
growth_par_beta[1,1,2] <- 3 # initial size increment  
growth_par_beta[1,2,2] <- 15 #last size increment : mean
growth_par_beta[1,3,2] <-  0.3 #  scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,2] <- 0.1 # initial size increment  : sd
growth_par_beta[2,2,2] <- 0.1 #last size increment : sd
growth_par_beta[2,3,2] <-  NA 

# ---- male immature
growth_par_beta[1,1,3] <- 3 # initial size increment  
growth_par_beta[1,2,3] <- 15 #last size increment : mean
growth_par_beta[1,3,3] <-  0.3 #  scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,3] <- 0.1 # initial size increment  : sd
growth_par_beta[2,2,3] <- 0.1 #last size increment : sd
growth_par_beta[2,3,3] <-  NA 

# ---- male mature
growth_par_beta[1,1,4] <- 3 # initial size increment  
growth_par_beta[1,2,4] <- 15 #last size increment : mean
growth_par_beta[1,3,4] <-  0.3 #  scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,4] <- 0.1 # initial size increment  : sd
growth_par_beta[2,2,4] <- 0.1 #last size increment : sd
growth_par_beta[2,3,4] <-  NA 


# 3-c- Uncertainty/sensitivity --------------------------------------------
# -------------------------------------------------------------------------
# -- 3-c-i Multiplicative or additive effect
ad_eff_growth <- TRUE # additive or multiplicative effect

# -- 3-c-ii Which preferential habitat function 
pref_hab_growth= c("quadr") # peferential habitat function (logi, quad, lin)
pars_pref_hab_growth = c(4,2)



# 3-d- Growth Settings ----------------------------------------------------
# -------------------------------------------------------------------------
source("spatial_projection/LHP_functions/pars_LHP_setting.R")

pars_Growth_setting_m_imm <- pars_LHP_setting(pref_hab_growth,
                                              pars_pref_hab_growth,
                                              x_omega,
                                              x_epsilon,
                                              scale_g,
                                              LHP,
                                              ad_eff_growth ,
                                              Growth_spatial,
                                              growth_par_beta[,,1],
                                              n_grpar_growth,
                                              Years_climsc)

pars_Growth_setting_m_mat <- pars_LHP_setting(pref_hab_growth,
                                              pars_pref_hab_growth,
                                              x_omega,
                                              x_epsilon,
                                              scale_g,
                                              LHP,
                                              ad_eff_growth ,
                                              Growth_spatial,
                                              growth_par_beta[,,2],
                                              n_grpar_growth,
                                              Years_climsc)

pars_Growth_setting_f_imm <- pars_LHP_setting(pref_hab_growth,
                                              pars_pref_hab_growth,
                                              x_omega,
                                              x_epsilon,
                                              scale_g,
                                              LHP,
                                              ad_eff_growth ,
                                              Growth_spatial,
                                              growth_par_beta[,,3],
                                              n_grpar_growth,
                                              Years_climsc)

pars_Growth_setting_f_mat <- pars_LHP_setting(pref_hab_growth,
                                              pars_pref_hab_growth,
                                              x_omega,
                                              x_epsilon,
                                              scale_g,
                                              LHP,
                                              ad_eff_growth ,
                                              Growth_spatial,
                                              growth_par_beta[,,4],
                                              n_grpar_growth,
                                              Years_climsc)
# 3-e - Generates Growth --------------------------------------------------
# -------------------------------------------------------------------------
source("spatial_projection/LHP_functions/growth.R")

growth_m_imm <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
growth_m_mat <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
growth_f_imm <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))
growth_f_mat <- array(0,dim=c(n_t,length(lon),length(lat),n_p,n_p))

#pb <- txtProgressBar(min = 1, max =  n_t, style = 3)
pb <- txtProgressBar(min = 1, max =  n_t, style = 3)
for(t in 1:n_t){
  #k= (c-1)*n_t + t
  if(molt_time_m[t] == 1){
    
    pars_Growth_setting_m_imm$year_LHP <- pars_Growth_setting_m_imm$Years_climsc[t]
    #pars_Growth_setting_m_imm$clim_sc <- pars_Growth_setting_m_imm$clim_sc_test[c]
    
    pars_Growth_setting_m_mat$year_LHP <- pars_Growth_setting_m_mat$Years_climsc[t]
    #pars_Growth_setting_m_mat$clim_sc <- pars_Growth_setting_m_mat$clim_sc_test[c]
    
    pars_Growth_setting_f_imm$year_LHP <- pars_Growth_setting_f_imm$Years_climsc[t]
    #pars_Growth_setting_f_imm$clim_sc <- pars_Growth_setting_f_imm$clim_sc_test[c]
    
    pars_Growth_setting_f_mat$year_LHP <- pars_Growth_setting_f_mat$Years_climsc[t]
    #pars_Growth_setting_f_mat$clim_sc <- pars_Growth_setting_f_mat$clim_sc_test[c]
    
    growth_m_imm[t,,,,] <- growth(sizes,
                                  binsize,
                                  pars_Growth_setting_m_imm,
                                  n_s,
                                  n_p,
                                  plot = FALSE)
    
    growth_m_mat[t,,,,] <- growth(sizes,
                                    binsize,
                                    pars_Growth_setting_m_mat,
                                    n_s,
                                    n_p,
                                    plot = FALSE)
    
    growth_f_imm[t,,,,] <- growth(sizes,
                                    binsize,
                                    pars_Growth_setting_f_imm ,
                                    n_s,
                                    n_p,
                                    plot = FALSE)
    
    
    growth_f_mat[t,,,,] <- growth(sizes,
                                    binsize,
                                    pars_Growth_setting_f_mat,
                                    n_s,
                                    n_p,
                                    plot = FALSE)
    Sys.sleep(0.1)
    setTxtProgressBar(pb,t)
    
    
  }
}

dim(growth_m_imm)
size_transition_mat_f_imm  <- growth_f_imm
size_transition_mat_f_mat  <- growth_f_mat
size_transition_mat_m_imm  <- growth_m_imm
size_transition_mat_m_mat  <- growth_m_mat  


# ------------------------------------------------------------------------- 
# -------------------------------------------------------------------------
# 4- Simulating molt ------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

Molt_spatial = FALSE

# molt
terminal_molt<-1
term_molt_prob<-c(0,0,0,.1,.2,.3,.4,.4,.4,.4,.4,.75,.9,1,1,1,1,1,1,1,1,1)
plot(term_molt_prob~sizes)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 5- Simulate fishery Selectivity -----------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

S_spatial = FALSE

#-- fishery pars
fish_50<-95
fish_95<-101
fish_sel<-1/(1+exp(-log(19)*(sizes-fish_50)/(fish_95-fish_50)))

#-- weight at length parameters
weight_a_f<-0.001
weight_b_f<-3

weight_a_m<-0.0012
weight_b_m<-3  

wt_at_len<-rbind(weight_a_f*sizes^weight_b_f,weight_a_m*sizes^weight_b_m)  

#-- fishery selectivity
fish_sel_50_f<-NA
fish_sel_95_f<-NA
fish_sel_50_m<-99
fish_sel_95_m<-101

fish_sel<-rbind(1/(1+exp(-log(19)*(sizes-fish_sel_50_f)/(fish_sel_95_f-fish_sel_50_f))),
                1/(1+exp(-log(19)*(sizes-fish_sel_50_m)/(fish_sel_95_m-fish_sel_50_m))))
fish_sel[is.na(fish_sel)]<-0

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 6 -Simulate movement ----------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#==movement at size
move_len_50<-60
move_len_95<-70

#==this function sets a dispersal kernel for every cell          
#==sd is 'one grid square'  
#source("movArray.R")
#  movement_dispersal<-movArray(SpaceR=length(lat),SpaceC=length(lon),sdx=1,sdy=1)        


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 7- Fisheries processes --------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# 7-a- Calculate costs to patch from the port for costs   ----------------
# -------------------------------------------------------------------------
#-- dutch harbor
port_lat<-  54
port_lon<- -166.54


# 7-b- calculate costs to fish --------------------------------------------
# -------------------------------------------------------------------------
#-- this should be related to the amount of fish in a patch
cost_fish<-10
cost_travel<-10000
price<-1.5
fishers<-2
quota<-rep(10000000,year_n)



fake_dist_data<-0

###########################################################################
# -------------------------------------------------------------------------
# STEP 2 RUN POP DYN
# -------------------------------------------------------------------------
###########################################################################

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 1- Initial states
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
initial_Ab <- initial_state(fake_dist_data, n_s, n_n, n_p, n_t,plot=FALSE) 


initial_Ab$imm_N_at_Len

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 2- Costs to fish
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# choose the location of the port 
# dutch harbor
port <- NULL
port$lat <- 54
port$lon<- -166.54

#==this should be related to the amount of fish in a patch
cost_fish<-10

cost_travel<-10000
cost_patch<-cost_travel*distance_map + cost_fish
price <- 1.5


# Here we should use output from assessment but for now we simulate -------
N_quota <- 10000000
fishers<-4
quota <- matrix(N_quota,fishers,year_n)


cost_to_fish_distance <- costs_to_fish(land_mask,
                                       port_lat,
                                       port_lon,
                                       lat,
                                       lon)











