###################################################
## Spatial MSE loop for Snow Crab of the Bering Sea
###################################################
## B. Alglave, M. Olmos, M. Veron, C. Szuwalski
rm(list=ls())

source("4_full_MSE/source/settings/libraries.R")

## Load and shape data for parameterzing and doing projection (GMACS, climatic projections)
#------------------------------------------------------------------------------------------
source("4_full_MSE/source/settings/load_data.r")

#--------------------------------------------------------------------------------------
#--------------------------------------- Settings -------------------------------------
#--------------------------------------------------------------------------------------

## Spatial extent
#----------------
source("4_full_MSE/source/settings/spatial_extent.r")

## Time series characteristics
#-----------------------------
# Time serie limit, seasonality 
# climatic scenarios
source("4_full_MSE/source/settings/time_series_settings.r")

## Demographic settings
#----------------------
source("4_full_MSE/source/settings/demographic_param.r")

## Fishery settings
#------------------
source("4_full_MSE/source/settings/fishery_param.R")

## Stock assessment
#------------------
source("4_full_MSE/source/settings/stock_assessment.R")

## Management settings
#---------------------
source("4_full_MSE/source/settings/proj_to_fish.R")

HCR = "cody"
fmsy_proxy = 1
bmsy_proxy = 183.1*1e9 / 1e3
b_hcr = 2
a_hcr = 0

## Matrices of abundance at size (initialization)
#------------------------------------------------
source("4_full_MSE/source/settings/ab_at_size.r")

## Dataframe for survey and commercial catch data
#------------------------------------------------
# From Max codes so that these are in good format to fit IPM model
source("4_full_MSE/source/settings/survey_commercial_data.r")

#--------------------------------------------------------------------------------------
#--------------------------------------- Project --------------------------------------
#--------------------------------------------------------------------------------------

## Stocking and debugging param for projections
source("4_full_MSE/source/settings/stocking_and_debug.r")

for(clim_sc in c("rcp45","rcp85","ssp126","ssp585")){
  
  R_stochastic = F
  simu_name = paste0("Clim_scen_",clim_sc,"-Nat_morta_",morta_model)
  
  #==indices: lat,lon,sex,size,time
  for(t in 1:(length(proj_period)-1))
    # for(t in 1:(length(proj_period)-1))
  {
    
    print(paste0(simu_name," | ",t))
    #==create a 'working' array for a given time step of both mature and immature critters
    temp_imm_N<-imm_N_at_Len[,,,,t]
    temp_mat_N<-mat_N_at_Len[,,,,t]
    
    x11()
    par(mfrow = c(4,6))
    for(i in 1:22){
      
      plot(mat_N_at_Len[,,2,i,t] * land_mask_na,
           main = paste0("Age ",i),asp = 1)
      
    }
    
    #==========================
    #== INIT FOR EACH YEAR ====
    #==========================
    # The environmental covariates are defined at yearly time step,
    # so at the moment we parameter life history traits at a yearly 
    # scale
    if(january_month[t] == T){
      
      source("4_full_MSE/source/projection/year_init.r")
      
    }
    
    if(print_distrib){plot(temp_mat_N[,,2,5],main="FIRST")}
    
    ## Follow up chunk (only male at the moment)
    follow_res = follow_ab_f(follow_ab_df,
                             simu_name,
                             phase = "init",
                             ab_imm_matrix = temp_imm_N[,,2,],
                             ab_mat_matrix = temp_mat_N[,,2,],
                             follow_ab_iter = follow_ab_iter,
                             t = t,
                             follow_ab = follow_ab)
    follow_ab_df = follow_res[[1]]
    follow_ab_iter = follow_res[[2]]
    
    #==========================
    #== FISHERY OCCURS ========
    #==========================
    if(fish_time[t]==1 & sum(quota) > 0)
    {
      for(f in 1:fishers)
      {
        
        source("4_full_MSE/source/projection/fishery.r")
        
      }
      
      if(print_distrib){plot(temp_mat_N[,,2,5],main="FISHERY OCCURED")}
      
      ## Follow up chunk (only male at the moment)
      follow_res = follow_ab_f(follow_ab_df,
                               simu_name,
                               phase = "post fishery",
                               ab_imm_matrix = temp_imm_N[,,2,],
                               ab_mat_matrix = temp_mat_N[,,2,],
                               follow_ab_iter = follow_ab_iter,
                               t = t,
                               follow_ab = follow_ab)
      follow_ab_df = follow_res[[1]]
      follow_ab_iter = follow_res[[2]]

    }
    
    
    
    #==========================
    #== GROWTH OCCURS =========
    #==========================
    if( molt_time[1,t]==1 | molt_time[2,t]==1 )
    {
      
      source("4_full_MSE/source/projection/growth.r")
      if(print_distrib){plot(temp_mat_N[,,2,5],main="GROWTH OCCURED")}
      
    }
    
    
    #==========================
    #== MOVEMENT OCCURS =======
    #==========================
    if(move_time[t]==1)
    {
      
      source("4_full_MSE/source/projection/movement.r")
      
      ## Follow up chunk (only male at the moment)
      follow_res = follow_ab_f(follow_ab_df,
                               simu_name,
                               phase = "post monthly movement",
                               ab_imm_matrix = temp_imm_N[,,2,],
                               ab_mat_matrix = temp_mat_N[,,2,],
                               follow_ab_iter = follow_ab_iter,
                               t = t,
                               follow_ab = follow_ab)
      follow_ab_df = follow_res[[1]]
      follow_ab_iter = follow_res[[2]]

    }
    
    #==========================
    #== SPAWNING OCCURS =======
    #==========================      
    #==this makes a map of spawning biomass to be used with transition matrices for recruitment
    if(mate_time[t]==1 )
    {
      #==aggregate spawning biomass
      #==just count female biomass?
      #==include sperm reserves?
      #==include biennial spawning?
      
      #==THIS IS JUST NUMBERS RIGHT NOW...
      spbiom<-apply(temp_mat_N[,,1,],c(1,2),sum,na.rm=T)
      plot_spb<-spbiom
      plot_spb[plot_spb==0]<-NA
      # filled.contour(x=lon,y=rev(lat),g(plot_spb),plot.axes=map(add=TRUE,fill=T,col='grey') )
    }
    
    
    #==========================
    #== RECRUITMENT OCCURS ====
    #========================== 
    #==how do we determine which bins they drop into?
    #==will temperature determine the size they reach in the time before they settle?
    if(recruit_time[t]==1)
    {
      
      source("4_full_MSE/source/projection/recruitment.r")
      
      if(print_distrib){plot(temp_mat_N[,,2,5],main="RECRUIT OCCURED")}
      
      ## Follow up chunk (only male at the moment)
      follow_res = follow_ab_f(follow_ab_df,
                               simu_name,
                               phase = "post recruitment",
                               ab_imm_matrix = temp_imm_N[,,2,],
                               ab_mat_matrix = temp_mat_N[,,2,],
                               follow_ab_iter = follow_ab_iter,
                               t = t,
                               follow_ab = follow_ab)
      follow_ab_df = follow_res[[1]]
      follow_ab_iter = follow_res[[2]]

    }
    
    
    #==========================
    #== Natural mortality =====
    #==========================
    source("4_full_MSE/source/projection/mortality.R")
    
    ## Follow up chunk (only male at the moment)
    follow_res = follow_ab_f(follow_ab_df,
                             simu_name,
                             phase = "post natural mortality",
                             ab_imm_matrix = imm_N_at_Len[,,2,,t+1],
                             ab_mat_matrix = mat_N_at_Len[,,2,,t+1],
                             follow_ab_iter = follow_ab_iter,
                             t = t,
                             follow_ab = follow_ab)
    follow_ab_df = follow_res[[1]]
    follow_ab_iter = follow_res[[2]]
    
    
    #==generate scientific (1 sample per cell grid - could be something else)
    if(survey_time[t] == 1){
      
      source("4_full_MSE/source/projection/sample_scientific.r")
      
    }
    
    
    #========================== 
    #==Stock assessment
    #========================== 
    # "spatialIPM": spatially-explicit model IPM
    # "nonspatialIPM": non spatial model IPM
    # "GMACS": standard stock assessment model
    # "none": no feedback loop
    if(SA_time[t] == 1){
      
      source("4_full_MSE/source/projection/make_commercial.r")
      
      source("4_full_MSE/source/projection/stock_assessment.r")
      
      source("4_full_MSE/source/projection/hcr.r")
      
    }
    
  }
  
  ## Save runs
  source("4_full_MSE/source/projection/save.r")
  
}

## Plot results
# source("4_full_MSE/source/projection/plot.r")

