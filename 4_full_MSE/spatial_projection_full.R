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

## Climate projection
#--------------------
source("4_full_MSE/source/settings/proj_param.r")

## Seasonality
#-------------
#==Binary vectors related to period that determine when life events happen
#==July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun

# Easier to manage for GMACS ----
survey_time_Yr   <-c(1,0,0,0,0,0,0,0,0,0,0,0)
fish_time_Yr		 <-c(0,0,0,0,0,0,0,1,1,1,0,0)
recruit_time_Yr	 <-c(0,0,0,0,0,0,0,0,0,1,0,0)
move_time_Yr		 <-c(1,1,1,1,1,1,1,1,1,1,1,1)
Fem_molt_time_Yr <- c(0,0,0,0,0,0,0,0,0,1,0,0)
Mal_molt_time_Yr <- c(0,0,0,0,0,0,0,0,0,1,0,0)
mate_time_Yr 	   <-c(0,0,0,0,0,0,0,1,0,1,0,0)
SA_time_Yr       <-c(0,0,0,0,0,0,0,0,0,1,0,0)
january_month_Yr <-c(1,0,0,0,0,0,0,0,0,0,0,0)

# ==========

survey_time   <-rep(survey_time_Yr,year_n)
fish_time		  <-rep(fish_time_Yr,year_n)
recruit_time	<-rep(recruit_time_Yr,year_n)
move_time		  <-rep(move_time_Yr,year_n)
molt_time  	  <-rbind(rep(Fem_molt_time_Yr,year_n),rep(Mal_molt_time_Yr,year_n)) # females first, males second
mate_time 	  <-rep(mate_time_Yr,year_n)
SA_time       <-rep(SA_time_Yr,year_n)
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

list_it = 0

## Objects to stock variables
#----------------------------
# Abundance list
imm_N_at_Len_full = list()
mat_N_at_Len_full = list()

# Fishery list
total_spatial_catch_full = list()
catch_by_fisher_full = list()
profit_by_fisher_full = list()
cost_by_fisher_full = list()
all_net_benefit_patch_full = list()
chosen_patch_full = list()
all_chosen_patch_full = list()

# Time-steps array
total_spatial_catch<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
total_spatial_catch_nb<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
catch_by_fisher<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period),fishers))
profit_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
cost_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
all_net_benefit_patch=array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
all_chosen_patch=array(0,dim=c(2,length(proj_period),fishers))

## For debugging
#---------------
## Fix some values so that the loop works fine
max_quota_it = 100 # maximum iteration for filling the quota (so that we are not blocked in a while loop)
use_harv = 1

## Following abundance throughout the loop
follow_ab = T
follow_ab_iter = 1
follow_ab_df = data.frame(phase = "0",
                          size = 0,
                          ab_male_mat = 0,
                          ab_male_imm = 0,
                          t = 1)

## Print distribution throughout the loop
print_distrib = F

## Print messages for identifying where are bugs
print_messages <- F


#==indices: lat,lon,sex,size,time
for(t in 1:(length(proj_period)-1))
  # for(t in 1:(length(proj_period)-1))
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
  
  ## Follow up chunk
  if(follow_ab){
    
    for(size in 1:length(sizes)){
      
      follow_ab_df[follow_ab_iter,"phase"] = "init"
      follow_ab_df[follow_ab_iter,"size"] = size
      follow_ab_df[follow_ab_iter,"ab_male_imm"] = sum(temp_imm_N[,,2,size])
      follow_ab_df[follow_ab_iter,"ab_male_mat"] = sum(temp_mat_N[,,2,size])
      follow_ab_df[follow_ab_iter,"t"] = t
      follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
      follow_ab_iter=follow_ab_iter+1
    }
    
  }
  
  
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
    
    ## Follow up chunk
    if(follow_ab){
      
      for(size in 1:length(sizes)){
        
        follow_ab_df[follow_ab_iter,"phase"] = "post fishery"
        follow_ab_df[follow_ab_iter,"size"] = size
        follow_ab_df[follow_ab_iter,"ab_male_imm"] = sum(temp_imm_N[,,2,size])
        follow_ab_df[follow_ab_iter,"ab_male_mat"] = sum(temp_mat_N[,,2,size])
        follow_ab_df[follow_ab_iter,"t"] = t
        follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
        follow_ab_iter=follow_ab_iter+1
        
      }
      
    }
    
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
    
    ## Follow up chunk
    if(follow_ab){
      
      for(size in 1:length(sizes)){
        
        follow_ab_df[follow_ab_iter,"phase"] = "post monthly movement"
        follow_ab_df[follow_ab_iter,"size"] = size
        follow_ab_df[follow_ab_iter,"ab_male_imm"] = sum(temp_imm_N[,,2,size])
        follow_ab_df[follow_ab_iter,"ab_male_mat"] = sum(temp_mat_N[,,2,size])
        follow_ab_df[follow_ab_iter,"t"] = t
        follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
        follow_ab_iter=follow_ab_iter+1
      }
      
    }
    
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
    
    ## Follow up chunk
    if(follow_ab){
      
      for(size in 1:length(sizes)){
        
        follow_ab_df[follow_ab_iter,"phase"] = "post recruitment"
        follow_ab_df[follow_ab_iter,"size"] = size
        follow_ab_df[follow_ab_iter,"ab_male_imm"] = sum(temp_imm_N[,,2,size])
        follow_ab_df[follow_ab_iter,"ab_male_mat"] = sum(temp_mat_N[,,2,size])
        follow_ab_df[follow_ab_iter,"t"] = t
        follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
        follow_ab_iter=follow_ab_iter+1
        
      }
      
    }
    
  }
  
  
  #==========================
  #== Natural mortality =====
  #==========================
  source("4_full_MSE/source/projection/mortality.R")

  ## Follow up chunk
  if(follow_ab){
    
    for(size in 1:length(sizes)){
      
      follow_ab_df[follow_ab_iter,"phase"] = "post natural mortality"
      follow_ab_df[follow_ab_iter,"size"] = size
      follow_ab_df[follow_ab_iter,"ab_male_imm"] = sum(imm_N_at_Len[,,2,size,t+1],na.rm=T)
      follow_ab_df[follow_ab_iter,"ab_male_mat"] = sum(mat_N_at_Len[,,2,size,t+1],na.rm=T)
      follow_ab_df[follow_ab_iter,"t"] = t
      follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
      follow_ab_iter=follow_ab_iter+1
      
    }
    
  }
  
  
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
  if(SA_time[t] == 1){
    
    source("4_full_MSE/source/projection/make_commercial.r")
    
    source("4_full_MSE/source/projection/stock_assessment.r")
    
    source("4_full_MSE/source/projection/hcr.r")
    
  }
  
}

## Plot results
source("4_full_MSE/source/projection/plot.r")

## Save codes
