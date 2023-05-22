###################################################
## spatial MSE loop for Snow Crab of the Bering Sea
###################################################
## B. Alglave, M. Olmos, M. Veron, C. Szuwalski
rm(list=ls())

## Load packages and make paths
#------------------------------
source("2_Max_spatial_projection/LHP_functions/libraries.R")
print_messages <- F

## Load and shape data for parameterzing and doing projection (GMACS, climatic projections)
#------------------------------------------------------------------------------------------
source("4_full_MSE/source/settings/load_data.r")

#--------------------------------------------------------------------------------------
#--------------------------------------- Settings -------------------------------------
#--------------------------------------------------------------------------------------

## Spatial extent
#----------------
source("4_full_MSE/source/settings/spatial_extent.r")

## Climate projection
#--------------------
source("4_full_MSE/source/settings/proj_param.r")

## Seasonality
#-------------
#==Binary vectors related to period that determine when life events happen
#==July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun
survey_time   <-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),year_n)
fish_time		  <-rep(c(0,0,0,0,0,0,0,1,1,1,0,0),year_n)
recruit_time	<-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)
move_time		  <-rep(c(1,1,1,1,1,1,1,1,1,1,1,1),year_n)
molt_time  	  <-rbind(rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n),
                      rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)) # females first, males second
mate_time 	  <-rep(c(0,0,0,0,0,0,0,1,0,1,0,0),year_n)
SA_time       <-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)
january_month <-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),year_n)

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
source("4_full_MSE/Management/proj_to_fish.R")

HCR = "cody"
fmsy_proxy = 0.5
bmsy_proxy = 183.1*1e9
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
list_it = 0

# Abundance
imm_N_at_Len_full = list()
mat_N_at_Len_full = list()

# Exploitation
total_spatial_catch_full = list()
catch_by_fisher_full = list()
profit_by_fisher_full = list()
cost_by_fisher_full = list()
all_net_benefit_patch_full = list()
chosen_patch_full = list()
all_chosen_patch_full = list()

list_info = data.frame(it = 0,
                       cost_travel = 0)

max_quota_it = 100 # maximum iteration for filling the quota
use_harv = 1

print_distrib = F

for(cost_travel in 1000){ # c(0,1000,1000*2)
  
  total_spatial_catch<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
  total_spatial_catch_nb<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
  catch_by_fisher<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period),fishers))
  profit_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  cost_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  all_net_benefit_patch=array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  all_chosen_patch=array(0,dim=c(2,length(proj_period),fishers))
  
  #==indices: lat,lon,sex,size,time
  for(t in 1:(length(proj_period)-1))  
  {
    
    print(t)
    #==create a 'working' array for a given time step of both mature and immature critters
    temp_imm_N<-imm_N_at_Len[,,,,t]
    temp_mat_N<-mat_N_at_Len[,,,,t]
    
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
    #==movement input as a .csv?
    #==movement constant?
    #==movement follows gradient?
    
    if(move_time[t]==1)
    {
      
      source("4_full_MSE/source/projection/movement.r")

    }
    
    #==========================
    #== SPAWNING OCCURS =======
    #==========================      
    #==this makes a map of spawning biomass to be used with transition matrices for recruitment
    
    if(mate_time[t]==1 )
    {
      #==aggregate spawnign biomass
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
      
    }
    
    
    #==========================
    #== Natural mortality =====
    #==========================
    source("4_full_MSE/source/projection/mortality.R")
    
    
    #==generate scientific (1 sample per cell grid - could be something else)
    # Data_Geostat
    if(survey_time[t] == 1){
      
      source("4_full_MSE/source/projection/sample_scientific.r")
      
    }
    
    
    #==Stock assessment
    # "spatialIPM": spatially-explicit model IPM
    # "nonspatialIPM": non spatial model IPM
    # "GMACS": standard stock assessment model
    if(SA_time[t] == 1){
      
      source("4_full_MSE/source/projection/make_commercial.r")
      
      source("4_full_MSE/source/projection/stock_assessment.r")
      
      source("4_full_MSE/source/projection/hcr.r")
      
    }
    
  }
  
  list_it = list_it + 1
  list_info[list_it,"it"] = list_it
  list_info[list_it,"cost_travel"] = cost_travel
  
  # Abundance
  imm_N_at_Len_full[[list_it]] = imm_N_at_Len
  mat_N_at_Len_full[[list_it]] = mat_N_at_Len
  
  # Catch list
  total_spatial_catch_full[[list_it]] = total_spatial_catch
  catch_by_fisher_full[[list_it]] = catch_by_fisher
  profit_by_fisher_full[[list_it]] = profit_by_fisher
  cost_by_fisher_full[[list_it]] = cost_by_fisher
  all_net_benefit_patch_full[[list_it]] = all_net_benefit_patch
  all_chosen_patch_full[[list_it]] = all_chosen_patch
  
}

source("4_full_MSE/source/projection/plot.r")
