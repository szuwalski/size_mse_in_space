rm(list=ls())
#==what is the first size class to be modeled?
#==this depend on how size at maturity changes
#==do crab mature after a set number of molts?
#==do they molt no matter what, just different increments?
#==OR do they molt the same size increment, but fewer times?
#==size dependent molting?
#==or start the model at the point that they are already only molting once a year
#==and the size they enter the model change based on the temperature during the time period
library(gdistance)
library(raster)
library(maps)
library(maptools)
source("function_LHP/script_growth_inc.R")

#==Set model variables
lat		   <-seq(70,51.5,-.25)
lon		   <-seq(-179,-155,.5)
lat		   <-seq(70,51.5,length.out=40)
lon		   <-seq(-179,-155,length.out=40)
#lat		   <-seq(38,22,-.25)
#lon		   <-seq(115,135,.5)
#lat		   <-seq(31,29,-.05)
#lon		   <-seq(120,124,.05)

binsize  <-5
sizes		 <-seq(27.5,132.5,binsize)
year_n	 <-20
year_step <-12
proj_period	<-seq(1,year_step*year_n)
sexN		 <-2
imm_fem_M    <-0.32
imm_male_M   <-0.32
mat_fem_M    <-0.26
mat_male_M   <-0.28
rec_sizes    <-5
prop_rec     <-c(.25,.4,.25,.075,0.025)

inits_year =2022
Years_climsc_temp <- ((rep(c(inits_year:(inits_year-1+year_n)),year_step)))
Years_climsc <- Years_climsc_temp[ order((rep(c(inits_year:(inits_year-1+year_n)),year_step)))]

#==Binary vectors related to period that determine when life events happen
#==July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun
survey_time   <-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),year_n)
fish_time		  <-rep(c(0,0,0,0,0,0,0,1,0,0,0,0),year_n)
recruit_time	<-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)
move_time		  <-rep(c(1,1,1,1,1,1,1,1,1,1,1,1),year_n)
molt_time  	  <-rbind(rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n),rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)) # females first, males second
mate_time 	  <-rep(c(0,0,0,0,0,0,0,1,0,1,0,0),year_n)

#==designate areas of potential habitat (i.e. not land)

data(wrld_simpl)

## Create a SpatialPoints object
set.seed(0)
point_expand <- expand.grid(lon, lat)  
pts <- SpatialPoints(point_expand, proj4string=CRS(proj4string(wrld_simpl)))
proj4string(wrld_simpl)<-CRS(proj4string(pts))
## Find which points fall over land
land <- !is.na(over(pts, wrld_simpl)$FIPS)

## Check that it worked
plot(wrld_simpl,xlim = c(min(lon), max(lon)), ylim = c(min(lat),max(lat)))
points(pts, col=1+land, pch=16)

#==need to track numbers at size by sex by maturity by location
imm_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
mat_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))

load(file="smooth_mat_N_at_Len_2017.RData") # smooth_mat
load(file="smooth_imm_N_at_Len_2017.RData") # smooth_imm

#==FIX THIS SO THAT THESE DISTRIBUTIONS ARE MADE TO TRUE DISTRIBUTIONS AND MULTIPLIED BY NUMBERS AT LENGTH FROM ASSESSMENT
imm_N_at_Len[,,,,1]<-exp(smooth_imm)
mat_N_at_Len[,,,,1]<-exp(smooth_mat)

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
#==growth pars
# set up for growth
source("function_LHP/LHP_function.R")

G_spatial = TRUE
clim_sc=c("rcp45") # climate scenario
pref_hab= c("quadr") # peferential habitat function (logi, quad, lin)
pars =  c(4,2)# c(4,2)# parameters to define the preferential habitat function

x_om_init =0# 0.1 # if additive effect, % of the mean of the initial growth increment represented by omega (field driven by climate sc.)
x_om_last=0#0.1 # if additive effect, % of the mean of the last growth increment represented by omega (field driven by climate sc.)
x_eps_last=0#0.05 # if additive effect, % of the mean of the last growth increment represented by epsilon (random field)
x_eps_init=0#0.05 # if additive effect, % of the mean of the intial growth increment represented by epsilon (random field.)
beta_scale= c(0.3)# scale parameter of the gamma distribution used to generate the growth transition matric 
binsize # bin size
sizes # size classes
m_last = c(25,0.1) # last size increment : meam, sd
m_init = c(9,0.1) # initial size increment  : meam, sd
scale_g <- 1.5 # scale of the RF

# molt
terminal_molt<-1
term_molt_prob<-c(0,0,0,.1,.2,.3,.4,.4,.4,.4,.4,.75,.9,1,1,1,1,1,1,1,1,1)
plot(term_molt_prob~sizes)

#==fishery pars
fish_50<-95
fish_95<-101
fish_sel<-1/(1+exp(-log(19)*(sizes-fish_50)/(fish_95-fish_50)))

#==weight at length parameters
weight_a_f<-0.001
weight_b_f<-3

weight_a_m<-0.0012
weight_b_m<-3  

wt_at_len<-rbind(weight_a_f*sizes^weight_b_f,weight_a_m*sizes^weight_b_m)  

#==movement at size
move_len_50<-60
move_len_95<-70

#==fishery selectivity
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

cost<-raster(nrow=length(lat), ncol=length(lon), 
             xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs="+proj=utm")
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
filled.contour(x=lon,y=rev(lat),g(distance_map*land_mask),plot.axes=c(map(add=TRUE,fill=T,col='grey'),
                                                                      points(y=port_lat,x=port_lon,pch=16,col='red')))
#write.csv(distance_map,'dist.csv')


#================================================
# calculate costs to fish
#=============================================
#==this should be related to the amount of fish in a patch
cost_fish<-10

cost_travel<-10000
cost_patch<-cost_travel*distance_map + cost_fish
price<-1.5

fishers<-2
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
  years =  Years_climsc[t]
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
            potential_catch<-potential_catch + temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
              temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
          }
        
        #==patch has less than needed to fill quota
        if(potential_catch<=quota_remaining)
        { 
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
  if(molt_time[1,t]==1 | molt_time[2,t]==1 )
  {
    #======================================
    #==make size transition matrix immature
    # transition matrix 
    n_p <- length(sizes)
    
    if(!G_spatial){
      # female
      size_transition_mat_f_imm_temp <- growth_trans(sizes,m_last[1],
                                                     m_init[1],
                                                     sizes[1],
                                                     sizes[n_p],
                                                     binsize,
                                                     beta_scale)
      
      size_transition_mat_f_imm <- array(0, dim = c(length(lat),length(lon),n_p,n_p))
      size_transition_mat_f_imm[,,1:n_p,1:n_p] <- size_transition_mat_f_imm_temp[1:n_p,1:n_p]
      
      # male
      size_transition_mat_m_imm_temp <- growth_trans(sizes,m_last[1],
                                                     m_init[1],
                                                     sizes[1],
                                                     sizes[n_p],
                                                     binsize,
                                                     beta_scale)
      size_transition_mat_m_imm <- array(0, dim = c(length(lat),length(lon),n_p,n_p))
      size_transition_mat_m_imm[,,1:n_p,1:n_p] <- size_transition_mat_m_imm_temp[1:n_p,1:n_p]
      
      
    }else{
      size_transition_mat_f_imm <- growth_map(G_spatial,lon, lat,clim_sc, pref_hab, years,pars,
                                              x_om_init,x_om_last,x_eps_last,x_eps_init,
                                              beta_scale,binsize,sizes, m_last, m_init,scale_g)
      
      size_transition_mat_m_imm <- growth_map(G_spatial,lon, lat,clim_sc, pref_hab, years,pars,
                                              x_om_init,x_om_last,x_eps_last,x_eps_init,
                                              beta_scale,binsize,sizes, m_last, m_init,scale_g) 
    }
    
    if(terminal_molt==0)
    { 
      #==make size transition matrix mature
      if(!G_spatial){
        # female
        size_transition_mat_f_mat_temp <- growth_trans(sizes,m_last[1],
                                                       m_init[1],
                                                       sizes[1],
                                                       sizes[n_p],
                                                       binsize,
                                                       beta_scale)
        
        size_transition_mat_f_mat <- array(0, dim = c(length(lat),length(lon),n_p,n_p))
        size_transition_mat_f_mat[,,1:n_p,1:n_p] <- size_transition_mat_f_mat_temp[1:n_p,1:n_p]
        
        # male
        size_transition_mat_m_mat_temp <- growth_trans(sizes,m_last[1],
                                                       m_init[1],
                                                       sizes[1],
                                                       sizes[n_p],
                                                       binsize,
                                                       beta_scale)
        size_transition_mat_m_mat <- array(0, dim = c(length(lat),length(lon),n_p,n_p))
        size_transition_mat_m_mat[,,1:n_p,1:n_p] <- size_transition_mat_m_mat_temp[1:n_p,1:n_p]
        
        
      }else{
        size_transition_mat_m_mat<-growth_map(G_spatial,lon, lat,clim_sc, pref_hab, years,pars,
                                              x_om_init,x_om_last,x_eps_last,x_eps_init,
                                              beta_scale,binsize,sizes, m_last, m_init,scale_g) 
        size_transition_mat_f_mat<-growth_map(G_spatial,lon, lat,clim_sc, pref_hab, years,pars,
                                              x_om_init,x_om_last,x_eps_last,x_eps_init,
                                              beta_scale,binsize,sizes, m_last, m_init,scale_g) 
      }
    }
    
    #======================================
    
    # bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
    for(x in 1:nrow(imm_N_at_Len[,,,,t]))
      for(y in 1:ncol(imm_N_at_Len[,,,,t]))
      {
        if(land_mask[x,y]!=0)
        {
          
          #==immature crab molt, some mature, some remain immature
          if(molt_time[1,t]==1)
          {
            tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm[x,y,,]
            temp_imm_N[x,y,1,]<-tmp_molt*(1-term_molt_prob)
            temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (term_molt_prob)*tmp_molt
            
          }
          if(molt_time[2,t]==1)
          {
            tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm[x,y,,]
            temp_imm_N[x,y,2,]<-tmp_molt*(1-term_molt_prob)
            temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (term_molt_prob)*tmp_molt
          }
          
          #== mature crab 
          if(terminal_molt==0)
          { 
            if(!is.na(match(molt_time[1,t],t)) ) 
              temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat[x,y,,]
            if(!is.na(match(molt_time[2,t],t)))
              temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat[x,y,,]   
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
  
  
}

tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len,c(5),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len,c(5),sum,na.rm=T)

par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
plot(log(tot_imm),type='l',las=1,xaxt='n')
lines(log(tot_mat),lty=2)
legend('topright',bty='n',lty=c(1,2),legend=c("Immature N","Mature N"))
plot(tot_catch,xaxt='n',las=1)
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost,xaxt='n',las=1)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit,xaxt='n',las=1)
legend('right',bty='n',legend=c("Total profits"))

Save=NULL
Save=list("imm_N_at_Len"=imm_N_at_Len,"mat_N_at_Len"=mat_N_at_Len,
          "catch_by_fisher"=catch_by_fisher,
          "cost_by_fisher"=cost_by_fisher,
          "profit_by_fisher"=profit_by_fisher )
save(Save,file="Save.RData")
getwd()