rm(list=ls())
#==what is the first size class to be modeled?
#==this depend on how size at maturity changes
#==do crab mature after a set number of molts?
#==do they molt no matter what, just different increments?
#==OR do they molt the same size increment, but fewer times?
#==size dependent molting?
#==or start the model at the point that they are already only molting once a year
#==and the size they enter the model change based on the temperature during the time period

source("2_Max_spatial_projection/LHP_functions/libraries.R")

library(gdistance)
library(maps)
library(maptools)
library(raster)
library(rnaturalearth)
library(reshape2)
library(sf)
library(tidyverse)

world_sf <- ne_countries(scale = "medium", returnclass = "sf")

print_messages <- F

# Paths
project_spatialIPM <- "../Spatial_snow_crab_2021/" # to go to working directory related to spatial IPM
DateFile = paste0(getwd(), "/2_Max_spatial_projection/Outputs/")
if(!dir.exists(DateFile)) dir.create(DateFile)

## Spatial extent
#----------------
lat		   <-seq(70,51.5,length.out=40)
lon		   <-seq(-179,-155,length.out=40)

source("2_Max_spatial_projection/LHP_functions/spatial_grid.R")
spatial_grid <- spatial_grid(lon,lat)
attach(spatial_grid)

# designate areas of potential habitat (i.e. not land)
data(wrld_simpl)
point_expand <- expand.grid(lon, lat)  
pts <- SpatialPoints(point_expand, proj4string=CRS(proj4string(wrld_simpl)))
proj4string(wrld_simpl)<-CRS(proj4string(pts))
## Find which points fall over land
land <- !is.na(over(pts, wrld_simpl)$FIPS)

# ## Check
# plot(wrld_simpl,xlim = c(min(lon), max(lon)), ylim = c(min(lat),max(lat)))
# points(pts, col=1+land, pch=16)


## Climate scenarios 
#-------------------
clim_sc =c("rcp45") # climate scenario to test


## Projection period
#-------------------
year_n	 <-50
year_step <-12
proj_period	<-seq(1,year_step*year_n)
Years_climsc <- rep(2022:(2022+year_n), each=year_step)
n_t <- length(proj_period)


## Size class settings
#---------------------
size_class_settings <- "rough"
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


## Sex and maturity
#------------------
sexN		 <-2
imm_fem_M    <-0.32
imm_male_M   <-0.32
mat_fem_M    <-0.26
mat_male_M   <-0.28
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

## Seasonality
#-------------
#==Binary vectors related to period that determine when life events happen
#==July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun
survey_time   <-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),year_n)
fish_time		  <-rep(c(0,0,0,0,0,0,0,1,1,1,0,0),year_n)
recruit_time	<-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)
move_time		  <-rep(c(1,1,1,1,1,1,1,1,1,1,1,1),year_n)
molt_time  	  <-rbind(rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n),rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)) # females first, males second
mate_time 	  <-rep(c(0,0,0,0,0,0,0,1,0,1,0,0),year_n)

## Load survey data for sampling locations
load(paste0(project_spatialIPM,'02_transformed_data/survey/Data_Geostat_4class.RData'))


## Stock assessment
#------------------
SA <- "spatialIPM"
# "spatialIPM": spatially-explicit model IPM
# "nonspatialIPM": non spatial model IPM
# "GMACS": standard stock assessment model


## Matrices of abundance at size
#-------------------------------
## need to track numbers at size by sex by maturity by location
imm_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
mat_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))

load(file="1_OM_Cody_version0/smooth_mat_N_at_Len_2017.RData") # smooth_mat
load(file="1_OM_Cody_version0/smooth_imm_N_at_Len_2017.RData") # smooth_imm

#==FIX THIS SO THAT THESE DISTRIBUTIONS ARE MADE TO TRUE DISTRIBUTIONS AND MULTIPLIED BY NUMBERS AT LENGTH FROM ASSESSMENT
if(size_class_settings=="fine"){
  imm_N_at_Len[,,,,1]<-exp(smooth_imm)
  mat_N_at_Len[,,,,1]<-exp(smooth_mat)
}

if(size_class_settings=="rough"){
  
  for(i in 1:4){
    
    imm_N_at_Len[,,1,i,1] <- exp(smooth_imm)[,,1,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    imm_N_at_Len[,,2,i,1] <- exp(smooth_imm)[,,2,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    mat_N_at_Len[,,1,i,1] <- exp(smooth_mat)[,,1,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    mat_N_at_Len[,,2,i,1] <- exp(smooth_mat)[,,2,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    
  }
  
  
}

imm_N_at_Len[imm_N_at_Len=="NaN"]<-0
mat_N_at_Len[mat_N_at_Len=="NaN"]<-0

imm_N_at_Len[imm_N_at_Len<0]<-0
mat_N_at_Len[mat_N_at_Len<0]<-0


#==This makes up a random distribution for the population
#==only used for testing
fake_dist_data<-0
if(fake_dist_data==1)
{
  #==set dummy initial distribution--take this from the survey for real
  #==this is only for getting mechanics down
  imm_N_at_Len[,,1,1,1]<-rnorm(dim(imm_N_at_Len[,,1,1,1])[1]*dim(imm_N_at_Len[,,1,1,1])[2],1000,100)
  imm_N_at_Len[,,2,1,1]<-rnorm(dim(imm_N_at_Len[,,1,1,1])[1]*dim(imm_N_at_Len[,,1,1,1])[2],1000,100)
  mat_N_at_Len[,,1,1,1]<-rnorm(dim(imm_N_at_Len[,,1,1,1])[1]*dim(imm_N_at_Len[,,1,1,1])[2],1000,100)
  mat_N_at_Len[,,2,1,1]<-rnorm(dim(imm_N_at_Len[,,1,1,1])[1]*dim(imm_N_at_Len[,,1,1,1])[2],1000,100)
  
  for(x in 2:length(sizes))
  {
    imm_N_at_Len[,,1,x,1]<-imm_N_at_Len[,,1,x-1,1]*exp(-M)
    imm_N_at_Len[,,2,x,1]<-imm_N_at_Len[,,2,x-1,1]*exp(-M)
    mat_N_at_Len[,,1,x,1]<-mat_N_at_Len[,,1,x-1,1]*exp(-M)
    mat_N_at_Len[,,2,x,1]<-mat_N_at_Len[,,2,x-1,1]*exp(-M)
  }
  
}

#==delete critters where there is land (this will be used for movement as well)
land_mask<-matrix(1,ncol=length(lon),nrow=length(lat),byrow=T)
for(x in 1:length(lat))
  for(y in 1:length(lon))
  {
    if(land[intersect(which(!is.na(match(pts$Var1,(lon)[y]))) ,which(!is.na(match( pts$Var2,(lat)[x]))))])
      land_mask[x,y]<-0
  }

#==ugh. this is dumb, but how to automate?
g<-function(m) t(m)[,nrow(m):1]

land_mask[25:40,20:40]
land_mask[31,32]<-0
land_mask[32,29]<-0
filled.contour(land_mask)
filled.contour(x=lon,y=rev(lat),g(imm_N_at_Len[,,1,1,1]))
#write.csv(land_mask,'landmask.csv')

#==ensures no critters on land
for(x in 1:2)
  for(y in 1:length(sizes))
  {
    imm_N_at_Len[,,x,y,1]<- imm_N_at_Len[,,x,y,1]*land_mask
    mat_N_at_Len[,,x,y,1]<- imm_N_at_Len[,,x,y,1]*land_mask
  }

#==SHOULD THERE ALSO BE A 'DEPTH MASK'???
#==MAYBE A 'FISHERY MASK'??  Depth should limit the fishery.
filled.contour(x=lon,y=rev(lat),g(imm_N_at_Len[,,1,1,1]))

#==create a file with bottom temperature for all time periods
avg_bot_tmp<-rep(c(1,2,3,3,3,2,1,0,0,0,0,1),year_n)
for(x in 1:length(proj_period))
{
  temp<-land_mask*rnorm(length(land_mask),avg_bot_tmp[x],.5)
  #write.csv(temp,paste("temp_data/bot_temp_",x,".csv",sep=""),row.names=FALSE)
}


#==========================
# Population processes 
#==========================
## Growth
#--------
growth_model <- "cody_model" # "max_model" "cody_model"
# Cody's model --> non-spatial life-history parameters
# Maxime's model --> spatially varying life-history parameters

if(growth_model == "cody_model"){
  
  #==Cody's growth parameterization for generation of non spatially varying growth parameters
  alpha_grow_f_imm<-4
  alpha_grow_m_imm<-7
  beta_grow_f_imm<-1.05
  beta_grow_m_imm<-1.1
  
  alpha_grow_f_mat<-4
  alpha_grow_m_mat<-7
  beta_grow_f_mat<-1.05
  beta_grow_m_mat<-1.1
  
  growth_sd_imm<-c(5,4)
  growth_sd_mat<-c(5,4)

}


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

terminal_molt<-1


## Molting probability
if(size_class_settings == "fine") term_molt_prob<-c(0,0,0,.1,.2,.3,.4,.4,.4,.4,.4,.75,.9,1,1,1,1,1,1,1,1,1)
if(size_class_settings == "rough") term_molt_prob<-c(0.01, 0.3,0.58,0.9)
plot(term_molt_prob~sizes)


## fishery pars
#--------------
fish_50<-95
fish_95<-101
fish_sel<-1/(1+exp(-log(19)*(sizes-fish_50)/(fish_95-fish_50)))

## weight at length parameters
#-----------------------------
weight_a_f<-0.001
weight_b_f<-3

weight_a_m<-0.0012
weight_b_m<-3  

wt_at_len<-rbind(weight_a_f*sizes^weight_b_f,weight_a_m*sizes^weight_b_m)  

## Movement at size
#------------------
move_len_50<-60
move_len_95<-70

## Fishery selectivity
#---------------------
fish_sel_50_f<-NA
fish_sel_95_f<-NA
fish_sel_50_m<-99
fish_sel_95_m<-101

fish_sel<-rbind(1/(1+exp(-log(19)*(sizes-fish_sel_50_f)/(fish_sel_95_f-fish_sel_50_f))),
                1/(1+exp(-log(19)*(sizes-fish_sel_50_m)/(fish_sel_95_m-fish_sel_50_m))))
fish_sel[is.na(fish_sel)]<-0

#==this function sets a dispersal kernel for every cell          
#==sd is 'one grid square'  
#source("movArray.R")
#  movement_dispersal<-movArray(SpaceR=length(lat),SpaceC=length(lon),sdx=1,sdy=1)        

#=====================================================  
#===CALCULATE COSTS TO A PATCH FROM THE PORT FOR COSTS  
#=====================================================

# dutch harbor
port_lat<-  54
port_lon<- -166.54

load("4_full_MSE/")

cost<-raster(nrow=length(lat), ncol=length(lon), 
<<<<<<< HEAD
             xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs="+proj=utm")

save(data = cost, "/4_full_MSE/")

=======
             xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
             crs="+proj=utm")
>>>>>>> 41b62322f3a6f98b77d6b408267b35741d5748a7
cost[]<-1
for(x in 1:nrow(cost))
  for(y in 1:ncol(cost))
    if(land_mask[x,y]==0) cost[x,y]<-100000

trCost <- transition(1/cost, min, directions=8)
trCostC <- geoCorrection(trCost, type="c")
trCostR <- geoCorrection(trCost, type="r")

## Create three points (representing three points in time series)
pnts <- cbind(x=c(port_lon, -155,-160), y=c(port_lat, 52,58))
costDistance(trCostC,pnts)

## Display results for one set of points
plot(cost)
plot(SpatialPoints(pnts), add=TRUE, pch=20, col="red")
plot(shortestPath(trCostC, pnts[1,], pnts[2,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCostC, pnts[1,], pnts[3,], output="SpatialLines"), add=TRUE)
plot(shortestPath(trCostC, pnts[2,], pnts[3,], output="SpatialLines"), add=TRUE)

#==find the distance from harbor to all points
distance_map<-matrix(ncol=length(lon),nrow=length(lat))
colnames(distance_map)<-lon
rownames(distance_map)<-(lat)
for(x in 1:length(lon))
  for(y in 1:length(lat))
  {
    if(land_mask[y,x]!=0)
    {
      pts <- cbind(x=c(port_lon, as.numeric(colnames(distance_map)[x])), y=c(port_lat,  as.numeric(rownames(distance_map)[y])))
      distance_map[y,x]<-costDistance(trCostC, pts[1,],pts[2,])
      if(distance_map[y,x]>10000)distance_map[y,x]<-NA
    }
  }
# filled.contour(x=lon,y=rev(lat),g(distance_map*land_mask),plot.axes=c(map(add=TRUE,fill=T,col='grey'),
#                                                                       points(y=port_lat,x=port_lon,pch=16,col='red')))
#write.csv(distance_map,'dist.csv')


#================================================
# calculate costs to fish
#=============================================
#==this should be related to the amount of fish in a patch
cost_fish<-10

cost_travel<-10000
cost_patch<-cost_travel*distance_map + cost_fish
price<-1.5

fishers<-4
quota<-rep(10000000,year_n)

#=============================================
# PROJJEEEECCCT
#============================================
total_spatial_catch<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
catch_by_fisher<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period),fishers))
profit_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
cost_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))

#==indices: lat,lon,sex,size,time
for(t in 1:(length(proj_period)-1))
  #for(t in 1:320)
{
  print(t)
  #==create a 'working' array for a given time step of both mature and immature critters
  
  temp_imm_N<-imm_N_at_Len[,,,,t]
  temp_mat_N<-mat_N_at_Len[,,,,t]
  #filled.contour(x=lon,y=rev(lat),g(temp_mat_N[,,1,5]),plot.axes=map(add=TRUE,fill=T,col='grey') )
  # if(survey_time[t]==1)
  #   collect_survey_data()
  
  #==========================
  #==FISHERY OCCURS
  #==========================
  #==indices: lat,lon,sex,size,time
  
  if(fish_time[t]==1)
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
      
      catch_patch<-apply(temp_catch,c(1,2),sum)
      catch_patch[catch_patch>quota[f]]<-quota[f] # this makes it so they don't travel a long way if they can get it close
      
      #filled.contour(x=lon,y=rev(lat),g(catch_patch),plot.axes=map(add=TRUE,fill=T,col='grey') )
      net_benefit_patch<-catch_patch*price-cost_patch
      
      #filled.contour(x=lon,y=rev(lat),g(net_benefit_patch),plot.axes=map(add=TRUE,fill=T,col='grey'),zlim=c(0,max(net_benefit_patch,na.rm=T)) )
      max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
      chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
      # filled.contour(x=lon,y=rev(lat),g(net_benefit_patch),
      #               plot.axes=c(map(add=TRUE,fill=T,col='grey'),
      #                          points(x=lon[chosen_patch[2]],y=lat[chosen_patch[1]],col=2,pch=16)),
      #                                zlim=c(0,max(net_benefit_patch,na.rm=T)) )
      # 
      #========================================================
      #==subtract catch from locations while quota is remaining
      print("Begin quota")
      while(quota_remaining>0.1 & net_benefit_patch[chosen_patch[1],chosen_patch[2]]>0)
      {
        #==find closest, highest value, fishable patch
        max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
        #==have to do this if there are two patches with identical net benefits
        chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
        
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
              
              
              temp_catch<- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
              
              catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
            }
          #sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
          quota_remaining<-quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
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
      print("End quota")
      
      profit_by_fisher[chosen_patch[1],chosen_patch[2],t,f] <- sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])*price - cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]
    }
    
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
    
  }
  
  #==========================
  #==GROWTH OCCURS
  #==========================
  
  ########################################################################################################################
  ######################################## This is where Maxime's work plug-in ###########################################
  ########################################################################################################################
  
  if( molt_time[1,t]==1 | molt_time[2,t]==1 )
  {
    
    if(growth_model == "max_model"){
      
      if(print_messages) print("growth Max")
      source("4_full_MSE/LHP/growth_t.R")
      
    }
    
    
    # bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
    for(x in 1:nrow(imm_N_at_Len[,,,,t]))
      for(y in 1:ncol(imm_N_at_Len[,,,,t]))
      {
        
        if(growth_model == "cody_model"){
          
          if(land_mask[x,y]!=0)
          {
            #==this is where Max's maps could be implemented
            #==plug temp into growth curve
            f_postmolt_imm<-alpha_grow_f_imm + beta_grow_f_imm*sizes
            m_postmolt_imm<-alpha_grow_m_imm + beta_grow_m_imm*sizes
            f_postmolt_mat<-alpha_grow_f_mat + beta_grow_f_mat*sizes
            m_postmolt_mat<-alpha_grow_m_mat + beta_grow_m_mat*sizes
            
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
        
        #==========================================================================================================
        
        if(growth_model == "max_model"){
          
          # bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
          for(x in 1:nrow(imm_N_at_Len[,,,,t]))
            for(y in 1:ncol(imm_N_at_Len[,,,,t]))
            {
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
                  
                  if(!is.na(match(molt_time[1,t],t)) ) 
                    temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat[t,x,y,,] ##### Is the indexing right ??????
                  if(!is.na(match(molt_time[2,t],t)))
                    temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat[t,x,y,,] ##### Is the indexing right ??????
                }
                
              }
              
            }
          
        }
        
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
    tmp_rec_1<- matrix(rnorm(length(imm_N_at_Len[,,1,1,1]),1,1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,1,1,1]))
    tmp_rec_1[tmp_rec_1<0]<-0
    tmp_rec_2<- matrix(rnorm(length(imm_N_at_Len[,,2,1,1]),1,1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,2,1,1]))
    tmp_rec_2[tmp_rec_2<0]<-0  
    
    for(r in 1:rec_sizes)
    {
      temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + imm_N_at_Len[,,1,r,1]*tmp_rec_1
      temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + imm_N_at_Len[,,2,r,1]*tmp_rec_2
    }
  }
  
  #==update dynamics
  imm_N_at_Len[,,1,,t+1] <-  temp_imm_N[,,1,]*exp(-imm_fem_M*1/year_step)
  imm_N_at_Len[,,2,,t+1] <-  temp_imm_N[,,2,]*exp(-imm_male_M*1/year_step)
  mat_N_at_Len[,,1,,t+1] <-  temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
  mat_N_at_Len[,,2,,t+1] <-  temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
  
  
  
  #==generate scientific data based on existing data
  # Data_Geostat
  
  
  #==Stock assessment
  # "spatialIPM": spatially-explicit model IPM
  # "nonspatialIPM": non spatial model IPM
  # "GMACS": standard stock assessment model
  
  
  if(SA == "spatialIPM"){
    
    source(paste0(project_spatialIPM,"03_spatial_model/run_model_mse.R"))
    
  }
  
  if(SA == "nonspatialIPM"){
    
  }
  
  if(SA == "GMACS"){
    
  }
  
  
}

tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len,c(5),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len,c(5),sum,na.rm=T)

x11()
par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,.1,1,1))
plot(tot_imm[100:length(tot_imm)],type='l',las=1,xaxt='n',ylim=c(0,9000000000))
lines(tot_mat[100:length(tot_imm)],lty=2)
legend('topright',bty='n',lty=c(1,2),legend=c("Immature N","Mature N"))
# plot(tot_catch,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total catch"))
# plot(tot_cost,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total cost"))
# plot(tot_profit,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total profits"))
plot(tot_catch[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16,ylim=c(0,60000000))
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit[(tot_profit>0)],las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total profits"))
