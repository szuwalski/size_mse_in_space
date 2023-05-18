###################################################
## spatial MSE loop for Snow Crab of the Bering Sea
###################################################
## B. Alglave and C. Szuwalski
rm(list=ls())

## Load packages and make paths
#------------------------------
source("2_Max_spatial_projection/LHP_functions/libraries.R")
print_messages <- F

## Load required data for parameterzing and doing projection (GMACS, climatic projections)
#-----------------------------------------------------------------------------------------
load("4_full_MSE/data/Snow_GMACS_repfile.Rdata") # Outputs provided by Mathieu from GMACS
load("4_full_MSE/data/Abundance_Recruitment_Assessment.Rdata")
load("4_full_MSE/data/NRS_vars_wide_op.Rdata") # Environmental covariates --> all variables
load("2_Max_spatial_projection/Climate_scenarios/mn_var_all.Rdata") # --> bottom temperature

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
    #==INIT FOR EACH YEAR
    #==========================
    # The environmental covariates are defined at yearly time step,
    # so at the moment we parameter life history traits at a yearly 
    # scale
    if(january_month[t] == T){
      
      ## Temperature data frame
      En.cond <- mn_var_all %>% filter(simulation %in% clim_sc, year %in% Years_climsc[t] ) %>% dplyr::select(year,simulation, latitude,longitude, val)
      En.cond <- En.cond %>% mutate(Lon = ifelse(longitude <= 180, longitude, longitude -360),Lat=latitude)
      dat_cond <-  data.frame(coords=En.cond[,c("Lon","Lat")],Year=En.cond$year, Temp=En.cond$val, Climate_sc = En.cond$simulation)
      colnames(dat_cond) <- c( "Lon","Lat","Year","Temp","Climate_sc")
      crs_LL = CRS(proj4string(wrld_simpl))
      dat_cond_sf <- st_as_sf(dat_cond, coords=c("Lon","Lat"), crs=crs_LL)
      
      ## Cold pool extent
      
      
      if(growth_model == "max_model"){
        
        print("init growth")
        source("4_full_MSE/LHP/growth_t.R")
        
      }
      
      if(morta_model == "max_model"){
        
        print("init morta")
        source("4_full_MSE/LHP/morta_t.R")
        
      }
      
      if(recruit_model == "max_model"){
        
        print("init recruit")
        source("4_full_MSE/LHP/recruit_t.R")
        
      }
      
    }
    
    
    
    if(print_distrib){plot(temp_mat_N[,,2,5],main="FIRST")}

    #==========================
    #==FISHERY OCCURS
    #==========================
    #==indices: lat,lon,sex,size,time
    if(fish_time[t]==1 & sum(quota) > 0)
    {
      for(f in 1:fishers)
      {
        quota_remaining <- quota[f]
        
        #==THISNEEDS TO BE FIXED====
        #==calculate net benefits by patch
        #==this needs to be based on the amount of quota available
        # if there is a patch nearby they can get their quota filled, they will go there
        #==there should be some relationship between cost and biomass in a cell?
        
        temp_catch<-array(dim=c(length(lat),length(lon),length(sizes)))
        for(sex in 1:2)
          for(x in 1:length(sizes))
          {
            temp_catch[,,x]<-temp_imm_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
              temp_mat_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
          }
        
        ## Choose fishing location
        catch_patch<-apply(temp_catch,c(1,2),sum)
        catch_patch[catch_patch>quota[f]]<-quota[f] # this makes it so they don't travel a long way if they can get it close < # We could modify this
        net_benefit_patch<-catch_patch*price-cost_patch
        all_net_benefit_patch[,,t,f] = net_benefit_patch
        
        if(fishing_process == "max_benefit_min_dist"){
          
          max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
          chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
          
        }else if(fishing_process == "stochastic"){
          
          prob_net = net_benefit_patch[which(!is.na(net_benefit_patch) & net_benefit_patch > 0)]
          chosen_patch<-which(net_benefit_patch == sample(prob_net,size=1,prob=prob_net),arr.ind=T)
          chosen_patch=chosen_patch[sample(1:nrow(chosen_patch),1),]
          
        }
        
        #========================================================
        #==subtract catch from locations while quota is remaining
        quota_it = 1
        while(quota_it < max_quota_it & quota_remaining>0.1 & net_benefit_patch[chosen_patch[1],chosen_patch[2]]>0)
        {
          
          quota_it = quota_it + 1
          
          if(fishing_process == "max_benefit_min_dist"){
            
            max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
            chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
            
          }else if(fishing_process == "stochastic"){
            
            prob_net = net_benefit_patch[which(!is.na(net_benefit_patch) & net_benefit_patch > 0)]
            chosen_patch<-which(net_benefit_patch == sample(prob_net,size=1,prob=),arr.ind=T)
            chosen_patch=chosen_patch[sample(1:nrow(chosen_patch),1),]
            
          }
          
          if(is.array(chosen_patch)) all_chosen_patch[,t,f] = chosen_patch[1,]
          if(is.integer(chosen_patch) & length(chosen_patch) == 2) all_chosen_patch[,t,f] = chosen_patch
          
          #==calculate total potential catch in a patch
          potential_catch<-0
          for(sex in 1:2)
            for(x in 1:length(sizes))
            {
              if(print_messages) print(paste0("sex=",sex," | sizes=",x," | ",
                                             " | temp_imm_N=",temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x],
                                             " | temp_mat_N=",temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x],
                                             " | potential catch= ",potential_catch))
              potential_catch<-potential_catch + temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
            }
          
          #==patch has less than needed to fill quota
          if(potential_catch<=quota_remaining)
          {
            
            if(print_messages) print(paste0("potential_catch<=quota_remaining | potential_catch=",potential_catch,"| quota_remaining=",quota_remaining,"| use_harv=",use_harv))
            
            for(sex in 1:2)
              for(x in 1:length(sizes))
              {
                
                total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                  temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
                
                total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                  temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x] + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
                
                temp_catch <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
                
                catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
                
              }
            cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]   <- cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f] + cost_patch[chosen_patch[1],chosen_patch[2]]        
            quota_remaining<-quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
            #==update temp array of n at len
            for(sex in 1:2)
              for(x in 1:length(sizes))
              {
                temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
              }
            
          }
          
          #==patch has more than needed to fill quota
          if(potential_catch>quota_remaining)
          {
            
            #==find harvest rate that would fill quota
            maxHarv<-1
            minHarv<-.0000001
            for(o in 1:25)
            {
              use_harv<-(maxHarv+minHarv)/2
              temp_cat<-0
              for(sex in 1:2)
                for(x in 1:length(sizes))
                {
                  temp_cat<-temp_cat+
                    temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
                    temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
                }
              if(temp_cat<quota_remaining)
                minHarv<-use_harv
              if(temp_cat>quota_remaining)
                maxHarv<-use_harv
            }
            
            if(print_messages) print(paste0("potential_catch<=quota_remaining | potential_catch=",potential_catch,"| quota_remaining=",quota_remaining,"| use_harv=",use_harv))
            
            temp_catch<-0
            for(sex in 1:2)
              for(x in 1:length(sizes))
              {
                total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                  temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv  

                total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                  temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x] + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
                
                temp_catch<- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
                  temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
                
                catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
              }
            #sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
            quota_remaining <- quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
            cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]   <- cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f] + cost_patch[chosen_patch[1],chosen_patch[2]]
            
            #==update temp array of n at len
            for(sex in 1:2)
              for(x in 1:length(sizes))
              {
                temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv
              }
          }
          
          #===update the net benefits in while loop
          temp_catch<-array(dim=c(length(lat),length(lon),length(sizes)))
          for(sex in 1:2)
            for(x in 1:length(sizes))
              temp_catch[,,x]<-temp_imm_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + temp_mat_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
          
          catch_patch<-apply(catch_patch,c(1,2),sum)
          catch_patch[catch_patch>quota_remaining]<-quota_remaining
          net_benefit_patch<-catch_patch*price - cost_patch
        }
        
        profit_by_fisher[chosen_patch[1],chosen_patch[2],t,f] <- sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])*price - cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]
        
      }
      
      if(print_distrib){plot(temp_mat_N[,,2,5],main="FISHERY OCCURED")}
      
    }
    
    
    #==========================
    #==GROWTH OCCURS
    #==========================
    if( molt_time[1,t]==1 | molt_time[2,t]==1 )
    {
      
      # bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
      for(x in 1:nrow(imm_N_at_Len[,,,,t])){
        for(y in 1:ncol(imm_N_at_Len[,,,,t]))
        {
          
          ## Cody's growth model
          #---------------------
          if(growth_model == "cody_model"){
            
            if(land_mask[x,y]!=0)
            {
              
              #======================================
              #==make size transition matrix immature
              size_transition_mat_m_imm<-matrix(ncol=length(sizes),nrow=length(sizes))
              size_transition_mat_f_imm<-matrix(ncol=length(sizes),nrow=length(sizes))
              
              for(n in 1:nrow(size_transition_mat_m_imm))
              {
                size_transition_mat_m_imm[n,]<-dnorm(x=sizes,mean=m_postmolt_imm[n],sd=growth_sd_imm[2])
                size_transition_mat_f_imm[n,]<-dnorm(x=sizes,mean=f_postmolt_imm[n],sd=growth_sd_imm[1])
              }
              
              #==ensure no crab shrinks in size after molting
              size_transition_mat_f_imm[lower.tri(size_transition_mat_f_imm)]<-0
              size_transition_mat_m_imm[lower.tri(size_transition_mat_m_imm)]<-0
              
              #==standardize to sum to 1
              for(z in 1:nrow(size_transition_mat_f_imm))
              {
                size_transition_mat_f_imm[z,]<-size_transition_mat_f_imm[z,]/sum(size_transition_mat_f_imm[z,],na.rm=T)
                size_transition_mat_m_imm[z,]<-size_transition_mat_m_imm[z,]/sum(size_transition_mat_m_imm[z,],na.rm=T)
              }
              
              #==immature crab molt, some mature, some remain immature
              if(molt_time[1,t]==1)
              {
                tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm
                temp_imm_N[x,y,1,]<-tmp_molt*term_molt_prob
                temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (1-term_molt_prob)*tmp_molt
                
              }
              if(molt_time[2,t]==1)
              {
                tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm
                temp_imm_N[x,y,2,]<-tmp_molt*term_molt_prob
                temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (1-term_molt_prob)*tmp_molt
              }
              
              if(terminal_molt==0)
              { 
                #======================================
                #==make size transition matrix mature
                size_transition_mat_m_mat<-matrix(ncol=length(sizes),nrow=length(sizes))
                size_transition_mat_f_mat<-matrix(ncol=length(sizes),nrow=length(sizes))
                
                for(n in 1:nrow(size_transition_mat_m_mat))
                {
                  size_transition_mat_m_mat[n,]<-dnorm(x=sizes,mean=m_postmolt_mat[n],sd=growth_sd_mat[2])
                  size_transition_mat_f_mat[n,]<-dnorm(x=sizes,mean=f_postmolt_mat[n],sd=growth_sd_mat[1])
                }
                
                #==ensure no crab shrinks in size after molting
                size_transition_mat_f_mat[lower.tri(size_transition_mat_f_mat)]<-0
                size_transition_mat_m_mat[lower.tri(size_transition_mat_m_mat)]<-0
                
                #==normalize
                for(z in 1:nrow(size_transition_mat_m_mat))
                {
                  size_transition_mat_f_mat[z,]<-size_transition_mat_f_mat[z,]/sum(size_transition_mat_f_mat[z,],na.rm=T)
                  size_transition_mat_m_mat[z,]<-size_transition_mat_m_mat[z,]/sum(size_transition_mat_m_mat[z,],na.rm=T)
                }
                if(!is.na(match(molt_time[1,t],t)) ) 
                  temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat
                if(!is.na(match(molt_time[2,t],t)))
                  temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat     
              }
              
            }
            
          }
          
          
          ## Maximes' growth model
          #-----------------------
          if(growth_model == "max_model"){
            
            if(land_mask[x,y]!=0)
            {
              
              #==immature crab molt, some mature, some remain immature
              if(molt_time[1,t]==1)
              {
                
                tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm[t,x,y,,] ##### Is the indexing right ??????
                temp_imm_N[x,y,1,]<-tmp_molt*term_molt_prob
                temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (1-term_molt_prob)*tmp_molt
                
              }
              
              if(molt_time[2,t]==1)
              {
                tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm[t,x,y,,] ##### Is the indexing right ??????
                temp_imm_N[x,y,2,]<-tmp_molt*term_molt_prob
                temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (1-term_molt_prob)*tmp_molt
              }
              
              if(terminal_molt==0)
              {
                
                if(!is.na(match(molt_time[1,t],t)))
                  temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat[t,x,y,,] ##### Is the indexing right ??????
                if(!is.na(match(molt_time[2,t],t)))
                  temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat[t,x,y,,] ##### Is the indexing right ??????
              }
              
            }
            
          }
          
        }
        
      }
      
      ## Ontogenic movement
      #--------------------
      # that's rough, but cannot make better at the moment
      # the best would be to make some kind of transition 
      # movement from adult to juvenile with Jim's model 
      
      which_ad_size = which(sizes < ad_size)
      which_ad_size = which_ad_size[length(which_ad_size)] + 1
      
      if(length(dim(size_transition_mat_m_imm)) == 2){
        last_size_mig = which(size_transition_mat_m_imm[which_ad_size,] < 0.001)[which_ad_size]
      }else if(length(dim(size_transition_mat_m_imm)) == 5){
        last_size_mig2 = c()
        for(x in 1:length(lat)){
          for(y in 1:length(lon)){
            if(init_juv[x,y] > 0){
              last_size_mig2 = c(last_size_mig2,which(size_transition_mat_m_imm[t,x,y,which_ad_size,] < 0.001)[which_ad_size])
            }
          }
        }
        last_size_mig = max(last_size_mig2)
      }
      
      for(s in which_ad_size:last_size_mig){
        print(s)
        ontog_mig_N_imm_female = 0
        ontog_mig_N_imm_male = 0
        ontog_mig_N_mat_female = 0
        ontog_mig_N_mat_male = 0
        for(x in 1:length(lon)){
          for(y in 1:length(lat)){
            if(init_juv[x,y] > 0){
              
              ontog_mig_N_imm_female = ontog_mig_N_imm_female + temp_imm_N[x,y,1,s]
              ontog_mig_N_imm_male = ontog_mig_N_imm_male + temp_imm_N[x,y,2,s]
              ontog_mig_N_mat_female = ontog_mig_N_mat_female + temp_mat_N[x,y,1,s]
              ontog_mig_N_mat_male = ontog_mig_N_mat_male + temp_mat_N[x,y,2,s]
              
              temp_imm_N[x,y,1,s] = 0
              temp_imm_N[x,y,2,s] = 0
              temp_mat_N[x,y,1,s] = 0
              temp_mat_N[x,y,2,s] = 0
              
            }
          }
        }
        
        Ab_ontog_mig_N_imm_female = sum(ontog_mig_N_imm_female)
        Ab_ontog_mig_N_imm_male = sum(ontog_mig_N_imm_male)
        Ab_ontog_mig_N_mat_female = sum(ontog_mig_N_mat_female)
        Ab_ontog_mig_N_mat_male = sum(ontog_mig_N_mat_male)
        
        ontog_mig_N_imm_female = Ab_ontog_mig_N_imm_female * init_adult
        ontog_mig_N_imm_male = Ab_ontog_mig_N_imm_male * init_adult
        ontog_mig_N_mat_female = Ab_ontog_mig_N_mat_female * init_adult
        ontog_mig_N_mat_male = Ab_ontog_mig_N_mat_male * init_adult
        
        temp_imm_N[,,1,s] = ontog_mig_N_imm_female
        temp_imm_N[,,2,s] = ontog_mig_N_imm_male
        temp_mat_N[,,1,s] = ontog_mig_N_mat_female
        temp_mat_N[,,2,s] = ontog_mig_N_mat_male
        
      }
      
      if(print_distrib){plot(temp_mat_N[,,2,5],main="GROWTH OCCURED")}
      
    }
    
    
    #==========================
    #==MOVEMENT OCCURS
    #==========================  
    #==movement input as a .csv?
    #==movement constant?
    #==movement follows gradient?
    
    if(move_time[t]==1)
    {
      #==two options: gaussian and temperature mediated
      #==gaussian
      #==create disperal kernel for each space
      for(sex in 1:2)
        for(x in 1:length(sizes))
        {
          
          if(sizes[x] < ad_size){
            
            mov_imm_temp = temp_imm_N[,,sex,x]
            mov_mat_temp = temp_mat_N[,,sex,x]
            
            mov_imm_v = mfraction_gg_juv %*% as.vector(t(mov_imm_temp))
            mov_mat_v = mfraction_gg_juv %*% as.vector(t(mov_imm_temp))
            
            mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
            mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
            
          }else{
            
            mov_imm_temp = temp_imm_N[,,sex,x]
            mov_mat_temp = temp_mat_N[,,sex,x]
            
            mov_imm_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp))
            mov_mat_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp))
            
            mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
            mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
            
          }
          
          temp_imm_N[,,sex,x] = mov_imm_temp_2
          temp_mat_N[,,sex,x] = mov_mat_temp_2
          
          # ## Check that movement happens
          # x11()
          # par(mfrow = c(4,4))
          # 
          # mov_imm_temp = temp_imm_N[,,sex,x]
          # mov_mat_temp = temp_mat_N[,,sex,x]
          # 
          # plot(mov_imm_temp,main=paste0("t = 0"),asp = 1)
          # 
          # for(t in 1:15){
          # 
          #   mov_imm_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp)) # mrate_gg %*%
          #   mov_mat_v = mfraction_gg_ad %*% as.vector(t(mov_mat_temp))  # mrate_gg %*%
          # 
          #   mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
          #   mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
          # 
          #   # print(which((mov_imm_temp != mov_mat_temp_2)))
          # 
          #   mov_imm_temp = mov_imm_temp_2
          #   mov_mat_temp = mov_mat_temp_2
          # 
          #   plot(mov_mat_temp,main=paste0("t = ",t),asp = 1)
          # 
          # }
          # 
          # # and compare with
          # test = matrix(stationary_g_ad,nrow=40,ncol=40,byrow = T)
          # plot(test)
          # 
          # # 
          # ## Or
          # land_mask_na = land_mask
          # land_mask_na[which(land_mask_na == 0)] = NA
          # for(t in 1:(12*3)){
          #   if(t %in% c(1,13,25)){
          #     x11()
          #     par(mfrow = c(4,3))
          #   }
          #   
          #   plot(mat_N_at_Len[,,2,5,t] * land_mask_na,main=paste0(t),asp = 1)
          # }
          
        }
      
    }
    
    #==========================
    #==SPAWNING OCCURS
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
    #==RECRUITMENT OCCURS
    #========================== 
    #==how do we determine which bins they drop into?
    #==will temperature determine the size they reach in the time before they settle?
    if(recruit_time[t]==1)
    {
      #==this is a dumb temporary fix
      #==this implants the original recruitment with some error
      #==ultimately we need an algorithm to determine location and intensity of recruitment
      #==teleconnections postdoc will hopefully help with this
      
      print('Recruit')
      
      if(recruit_model == "cody_model"){
        
        meanlog_recruit_male_2 = meanlog_recruit_male
        meanlog_recruit_female_2 = meanlog_recruit_female
        
      }else if(recruit_model == "max_model"){
        
        source("4_full_MSE/LHP/recruit_t.R")
        
      }
      
      aggreg_rec_1 = rlnorm(n = 1, mean = meanlog_recruit_female_2, sdlog = sdlog_female)
      aggreg_rec_2 = rlnorm(n = 1, mean = meanlog_recruit_male_2, sdlog = sdlog_male)
      
      tmp_rec_1 <- init_juv * aggreg_rec_1
      tmp_rec_2 <- init_juv * aggreg_rec_2

      for(r in 1:rec_sizes)
      {
        temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + imm_N_at_Len[,,1,r,1]*tmp_rec_1*prop_rec[r]
        temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + imm_N_at_Len[,,2,r,1]*tmp_rec_2*prop_rec[r]
      }
      
      if(print_distrib){plot(temp_mat_N[,,2,5],main="RECRUIT OCCURED")}
      
    }
    
    #==update dynamics
    if(morta_model == "cody_model"){
      
      print("Cody morta")
      imm_N_at_Len[,,1,,t+1] <-  temp_imm_N[,,1,]*exp(-imm_fem_M*1/year_step)
      imm_N_at_Len[,,2,,t+1] <-  temp_imm_N[,,2,]*exp(-imm_male_M*1/year_step)
      mat_N_at_Len[,,1,,t+1] <-  temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
      mat_N_at_Len[,,2,,t+1] <-  temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
      
    }else if(morta_model == "max_model"){
      
      print("Max morta")
      imm_N_at_Len[,,1,,t+1] <- temp_imm_N[,,1,]*exp(-imm_fem_M_size*1/year_step)
      imm_N_at_Len[,,2,,t+1] <- temp_imm_N[,,2,]*exp(-imm_male_M_size*1/year_step)
      mat_N_at_Len[,,1,,t+1] <- temp_mat_N[,,1,]*exp(-mat_fem_M_size*1/year_step)
      mat_N_at_Len[,,2,,t+1] <- temp_mat_N[,,2,]*exp(-mat_male_size*1/year_step)
      
    }else if(morta_model == "coldpool"){
      
      df = data.frame(coldpool = coldpool_range_t)
      
      imm_fem_M_size_cp = exp(predict(mod,df))
      imm_male_M_size_cp = exp(predict(mod,df))
      
      imm_N_at_Len[,,1,,t+1] <- temp_imm_N[,,1,]*exp(-imm_fem_M_size_cp*1/year_step)
      imm_N_at_Len[,,2,,t+1] <- temp_imm_N[,,2,]*exp(-imm_male_M_size_cp*1/year_step)
      mat_N_at_Len[,,1,,t+1] <- temp_mat_N[,,1,]*exp(-mat_fem_M_size*1/year_step)
      mat_N_at_Len[,,2,,t+1] <- temp_mat_N[,,2,]*exp(-mat_male_size*1/year_step)
      
    }
    
    #==generate scientific (1 sample per cell grid - could be something else)
    # Data_Geostat
    if(survey_time[t] == 1){
      
      if(print_messages) print("Simulate scientific data")
      ## Implemented through matrix because much faster
      sci_it = 0
      for(x in 1:length(lat))
      {
        for(y in 1:length(lon))
        {
          for(i in 1:length(sizes))
          {
            
            include_loc = T
            if(sp_domain == "EBS" & EBS_mask_matrix[x,y] == 1) include_loc = F
            
            if(land_mask[x,y]!=0 & include_loc)
            {
              
              # Need to implement observation pdf
              # print(paste0("x:",x,"|y:",y,"|i:",i))
              if(sci_it == 0){
                
                sci_mat = matrix(data = c(i, # size_class
                                          Years_climsc[t], # year
                                          imm_N_at_Len[x,y,2,i,t] + mat_N_at_Len[x,y,2,i,t], # Catch_N: only male
                                          cell_area[x,y], # AreaSwept_km2
                                          0, # Vessel
                                          lat[x], # Lat
                                          lon[y]), # Lon
                                 nrow = 1)
                
                sci_it = sci_it + 1
                
              }else{
                
                line_vec = matrix(data = c(i, # size_class
                                           Years_climsc[t], # year
                                           imm_N_at_Len[x,y,2,i,t] + mat_N_at_Len[x,y,2,i,t], # Catch_N: only male
                                           cell_area[x,y], # AreaSwept_km2
                                           0, # Vessel
                                           lat[x], # Lat
                                           lon[y]), # Lon
                                  nrow = 1)
                sci_mat = rbind(sci_mat,line_vec)
                
              }
            }
          }
        }
      }
      
      sci_df = as.data.frame(sci_mat)
      colnames(sci_df) = colnames(Data_Geostat)
      sci_df$size_class = as.factor(sci_df$size_class)
      Data_Geostat = rbind(Data_Geostat,sci_df)
      
    }
    
    
    #==Stock assessment
    # "spatialIPM": spatially-explicit model IPM
    # "nonspatialIPM": non spatial model IPM
    # "GMACS": standard stock assessment model
    if(SA_time[t] == 1){
      
      if(print_messages) print("Make commercial data")
      
      com_it = 0
      for(x in 1:length(lat))
      {
        for(y in 1:length(lon))
        {
          for(i in 1:length(sizes))
          {
            
            include_loc = T
            if(sp_domain == "EBS" & EBS_mask_matrix[x,y] == 1) include_loc = F
            
            if(land_mask[x,y]!=0 & include_loc)
            {
              
              #------------------------------------------------------------------------------------------------------------------------
              # Need to implement observation pdf
              # print(paste0("x:",x,"|y:",y,"|i:",i))
              if(com_it==0){
                
                catch_mat = matrix(data = c(Years_climsc[t], # Year
                                            i, # size_bin
                                            sum(total_spatial_catch[x,y,2,i,which(Years_climsc == Years_climsc[t])]), # Catches_N: only male
                                            lon[y], # Lon
                                            lat[x]), # Lat
                                   nrow = 1)
                com_it = com_it + 1
                
              }else{
                
                line_vec = matrix(data = c(Years_climsc[t], # Year
                                           i, # size_bin
                                           sum(total_spatial_catch[x,y,2,i,which(Years_climsc == Years_climsc[t])]), # Catches_N: only male
                                           lon[y], # Lon
                                           lat[x]), # Lat
                                  nrow = 1)
                
                catch_mat = rbind(catch_mat,line_vec)
                
              }
              
            }
          }
        }
      }
      
      catch_df = as.data.frame(catch_mat)
      colnames(catch_df) = colnames(catch_N)
      catch_df$Year = as.character(catch_df$Year)
      catch_df$size_bin = as.character(catch_df$size_bin)
      catch_N = rbind(catch_N,catch_df)
      
      catch_N %>% 
        filter(Year > 2018) %>% 
        group_by(Year) %>% 
        dplyr::summarise(Catches_N = sum(Catches_N))
      
      if(print_messages) print("Make Stock assessment")
      
      if(SA == "spatialIPM"){
        
        if(display_print) print("spatial IPM")
        
        source(paste0(project_spatialIPM,"03_spatial_model/run_model_mse.R"))
        
      }
      
      if(SA == "nonspatialIPM"){
        
      }
      
      if(SA == "GMACS"){
        
      }
      
      if(SA == "none"){
        
        n_at_size = c()
        c_at_size = c()
        f_at_size = c()
        for(s in 1:length(sizes)){
          
          # c: catch, m: natural mortality, n: abundance
          # --> c = f / (f + m) * n * (1 - exp(-(f + m)))
          imm_male_M = 0.32 
          mat_male_M = 0.28
          m = 0.32
          
          c = sum(total_spatial_catch_nb[,,2,s,(t-9):t])
          
          ab_at_t = c()
          catch_at_t = c()
          ab_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
          catch_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
          for(t2 in (t-9):t){
            
            # Abundance
            ab_at_t0 = sum(imm_N_at_Len[,,2,s,t2] + mat_N_at_Len[,,2,s,t2])
            ab_at_t = c(ab_at_t,ab_at_t0)
            
            ab_at_t0_mat = imm_N_at_Len[,,2,s,t2] + mat_N_at_Len[,,2,s,t2]
            ab_at_t_mat = ab_at_t_mat + ab_at_t0_mat
            
            # Catch
            catch_at_t0 = sum(total_spatial_catch_nb[,,2,s,t2])
            catch_at_t = c(catch_at_t,catch_at_t0)
            
            catch_at_t0_mat = total_spatial_catch_nb[,,2,s,t2]
            catch_at_t_mat = catch_at_t_mat + catch_at_t0_mat
            
          }
          
          n = ab_at_t[1]
          n_mat = ab_at_t_mat / 9
          
          # Aggregated value
          n_plus_1 = (n - c * exp(m/2)) / exp(m)  # Pope's approximation (see https://www.fao.org/3/x9026e/x9026e06.htm)
          f = log(n / n_plus_1) - m
          
          # # Spatial value
          # n_plus_1_mat = (n_mat - catch_at_t_mat * exp(m/2)) / exp(m)  # Pope's approximation (see https://www.fao.org/3/x9026e/x9026e06.htm)
          # f = log(n_mat / n_plus_1_mat) - m
          
          n_at_size = c(n_at_size,n)
          c_at_size = c(c_at_size,c)
          f_at_size = c(f_at_size,f)
          
        }

      }
      
      ## Implement HCR for next year quota
      if(HCR == "cody"){
        
        f_mort = rep(0,length(sizes))
        temp_imm = apply(imm_N_at_Len[,,2,,t],c(3),sum,na.rm=T)
        temp_mat = apply(mat_N_at_Len[,,2,,t],c(3),sum,na.rm=T)
        fish_sel_hcr = fish_sel[2,]
        weight_at_size = wt_at_len[2,]
        functional_mat = 0
        mmb_def = 'morphometric'
        
        source("4_full_MSE/Management/hcr_cody.R")
        
      }
      
      
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

tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len,c(5),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len,c(5),sum,na.rm=T)

x11()
par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,.1,1,1))
plot(tot_imm,type='l',las=1,xaxt='n',ylim=c(0,max(tot_imm,tot_mat)))
lines(tot_mat,lty=2)
legend('topright',bty='n',lty=c(1,2),legend=c("Immature N","Mature N"))
# plot(tot_catch,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total catch"))
# plot(tot_cost,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total cost"))
# plot(tot_profit,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total profits"))
plot(tot_catch[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16,ylim=c(0,max(tot_catch)))
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit[(tot_profit>0)],las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total profits"))

x11()
par(mfrow = c(4,3))
for(x in 1:12){
  
  plot(imm_N_at_Len[,,2,x,600] * land_mask_na,
       main = paste0("Size: [",sizes[x],", ",sizes[x+1],"] mm"),
       asp=1)
  
}

