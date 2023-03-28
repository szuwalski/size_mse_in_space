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
library(FishStatsUtils)
library(maps)
library(maptools)
library(raster)
library(rnaturalearth)
library(reshape2)
library(sf)
library(spdep)
library(stars)
library(tidyverse)
library(TMB)

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
sp_domain = "EBS"

source("2_Max_spatial_projection/LHP_functions/spatial_grid.R")
spatial_grid <- spatial_grid(lon,lat)
attach(spatial_grid)

# designate areas of potential habitat (i.e. not land)
data(wrld_simpl)
point_expand <- expand.grid(lon, lat)
point_expand$key = 1:nrow(point_expand)
# plot(point_expand[,1],point_expand[,2])
# text(point_expand[,1],point_expand[,2],labels = 1:nrow(point_expand))
pts <- SpatialPoints(point_expand, proj4string=CRS(proj4string(wrld_simpl)))
proj4string(wrld_simpl)<-CRS(proj4string(pts))
## Find which points fall over land
land <- !is.na(over(pts, wrld_simpl)$FIPS)

# Compute cell area
test = rasterFromXYZ(cbind(point_expand,rep(rnorm(nrow(point_expand)),nrow(point_expand))), crs=CRS(proj4string(wrld_simpl)), digits=5)
cell_area = raster::area(test) %>% as.matrix()

## EBS area
load("4_full_MSE/data/EBS.RData")
xys = st_as_sf(as.data.frame(Extrapolation_List$Data_Extrap), coords=c("Lon","Lat"),crs=CRS(proj4string(wrld_simpl)))
EBS_sf = xys %>% 
  group_by() %>% 
  summarise(Include = sum(Include)) %>% 
  st_cast("MULTIPOINT") %>% 
  st_cast("MULTILINESTRING") %>% 
  st_cast("MULTIPOLYGON")

EBS_sf = st_convex_hull(EBS_sf)
# plot(EBS_sf)
EBS_sp = as_Spatial(EBS_sf)
EBS_mask = is.na(over(pts, EBS_sp)$Include)

# ## Check
# plot(wrld_simpl,xlim = c(min(lon), max(lon)), ylim = c(min(lat),max(lat)))
# points(pts, col=1+land, pch=16)

## Climate scenarios 
#-------------------
clim_sc = c("rcp45") # climate scenario to test

## Projection period
#-------------------
year_n	 <- 50
year_step <- 12
proj_period	<- seq(1,year_step*year_n)
Years_climsc <- rep(2022:(2022+year_n), each=year_step)
n_t <- length(proj_period)

## Load GMACS outputs for parameterizing
#---------------------------------------
load(file = "4_full_MSE/data/Snow_GMACS_repfile.Rdata")

## Size class settings
#---------------------
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
SA_time       <-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),year_n)

## Survey data
#-------------
## Load survey data for sampling locations
load('4_full_MSE/data/Data_Geostat_4class.RData')

DF <- as_tibble(Data_Geostat)

Data_Geostat = cbind(
  "size_class" = DF[, "size_bin"],
  "year" = DF[, "Year"],
  "Catch_N" = DF[, "Catch"],
  "AreaSwept_km2" = DF[, "AreaSwept_km2"],
  "Vessel" = 0,
  "Lat" = DF[, "Lat"],
  "Lon" = DF[, "Lon"]
)

colnames(Data_Geostat) <-
  c("size_class",
    "year",
    "Catch_N",
    "AreaSwept_km2",
    "Vessel",
    "Lat",
    "Lon")

Data_Geostat = Data_Geostat %>%
  filter(year %in% 2015:2016)

## Catch data
#------------
load(paste0(project_spatialIPM,"02_transformed_data/movement/COG/COG_smoother_ADFG/catch_fishery_mov_intersect.RData"))
catch_N <- catch_fishery_mov_intersect

# Choose if we implement fisheries catches (Fisheries_catches <-TRUE) or not in the model 
Fisheries_catches <-TRUE
mov <- FALSE

# Choose if we implement simulated fisheries catches(Catches_sim_random <-TRUE) ro not in the model 
Catches_sim_random <- FALSE
Catches_sim_fractionAb <- FALSE

heterMov <- FALSE
homoMov <- FALSE

# Smoother
#smoother <- TRUE # smoother is with knots
#ADFG <- TRUE # smoother is with ADFG cells
#KNOT_i <- FALSE


## Stock assessment
#------------------
SA <- "GMACS"
# "spatialIPM": spatially-explicit model IPM
# "nonspatialIPM": non spatial model IPM
# "GMACS": standard stock assessment model

## Spatial settings
Data_Set <- 'Snow_crab'

## Number of knot/stations
n_x = c(30, 50, 75, 100, 150, 200, 300)[3]

# Output from Calc_Kmeans with knots for a triangulated mesh
# Calc_Kmeans determines the location for a set of knots for
# approximating spatial variation
# n_x: the number of knots to select
# nstart the number of times that the k-means algorithm is run while
# searching for the best solution (default=100)
# ter.max the number of iterations used per k-means algorithm
Kmeans_Config = list("randomseed" = 1,
                     "nstart" = 100,
                     "iter.max" = 1e3)

# Define studied region
# builds an object used to determine areas to extrapolation
# densities to when calculating indices
strata.limits <- data.frame('STRATA' = "All_areas")
Region = "Eastern_Bering_Sea"
Extrapolation_List = FishStatsUtils::make_extrapolation_info(Region = Region,
                                                             strata.limits = strata.limits)


## Matrices of abundance at size
#-------------------------------
## need to track numbers at size by sex by maturity by location
imm_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
mat_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))

load(file="1_OM_Cody_version0/smooth_mat_N_at_Len_2017.RData") # smooth_mat
load(file="1_OM_Cody_version0/smooth_imm_N_at_Len_2017.RData") # smooth_imm

## Initialization data
ebs_2019 <- read_csv("4_full_MSE/data/ebs_2019.csv")
nbs_2019 <- read_csv("4_full_MSE/data/nbs_2019.csv")

ebs_2019_2 = ebs_2019 %>% 
  mutate(area = "ebs") %>% 
  mutate(stage = ifelse(WIDTH < 45,"Juv",'Mat')) %>% 
  dplyr::group_by(HAUL,MID_LATITUDE,MID_LONGITUDE,stage,area) %>% 
  tally()

nbs_2019_2 = nbs_2019 %>% 
  mutate(area = "nbs") %>% 
  mutate(stage = ifelse(WIDTH < 45,"Juv",'Mat')) %>% 
  dplyr::group_by(HAUL,MID_LATITUDE,MID_LONGITUDE,stage,area) %>% 
  tally()

count_df = rbind(ebs_2019_2,nbs_2019_2)

lat_range = range(count_df$MID_LATITUDE)
lon_range = range(count_df$MID_LONGITUDE)

Juv_plot = ggplot(count_df[which(count_df$stage == "Juv"),])+
  geom_point(aes(x=MID_LONGITUDE,y=MID_LATITUDE,col=log(n)),size = 2)+
  scale_color_distiller(palette = "Spectral")+
  ggtitle("Juveniles")+
  xlim(lon_range)+ylim(lat_range) # + facet_wrap(.~area)

Mat_plot = ggplot(count_df[which(count_df$stage == "Mat"),])+
  geom_point(aes(x=MID_LONGITUDE,y=MID_LATITUDE,col=log(n)),size = 2)+
  scale_color_distiller(palette = "Spectral")+
  ggtitle("Mature")+
  xlim(lon_range)+ylim(lat_range) # + facet_wrap(.~area)

point_sf = st_as_sf(point_expand,coords = c("Var1","Var2"))
raster_dom = st_rasterize(point_sf %>% dplyr::select(key, geometry))
grid_sf = st_as_sf(raster_dom)

count_sf = st_as_sf(count_df,coords = c("MID_LONGITUDE","MID_LATITUDE"))
test = st_intersection(grid_sf,count_sf) %>%
  as.data.frame %>% 
  group_by(key,stage) %>%
  dplyr::summarise(n = mean(n)) %>% 
  full_join(grid_sf) %>%
  filter(!is.na(n)) %>% 
  filter(!is.na(stage)) %>% 
  st_as_sf

ggplot()+
  geom_sf(data=test,aes(fill=n))+
  scale_fill_distiller(palette="Spectral")+
  facet_wrap(.~stage)

init_adult<-array(0,dim=c(length(lat),length(lon)))
init_juv<-array(0,dim=c(length(lat),length(lon)))

for(x in 1:length(lat))
  for(y in 1:length(lon))
  {
    
    key = point_expand$key[which(point_expand$Var1 == lon[y] & point_expand$Var2 == lat[x])]
    test_mat = test$n[test$key == key & test$stage == "Mat"]
    test_immat = test$n[test$key == key & test$stage == "Juv"]
    if(length(test_mat) > 0) init_adult[x,y] = test$n[test$key == key & test$stage == "Mat"]
    if(length(test_immat) > 0) init_juv[x,y] = test$n[test$key == key & test$stage == "Juv"]

  }

init_adult = init_adult / sum(init_adult)
init_juv = init_juv / sum(init_juv)

# cowplot::plot_grid(Juv_plot,Mat_plot)

if(size_class_settings=="fine"){
  imm_N_at_Len[,,,,1]<-array(0,dim = c(length(lat),length(lon),2,length(sizes)))
  mat_N_at_Len[,,,,1]<-array(0,dim = c(length(lat),length(lon),2,length(sizes)))
}

if(size_class_settings=="rough"){
  
  for(i in 1:4){
    
    # imm_N_at_Len[,,1,i,1] <- exp(smooth_imm)[,,1,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    # imm_N_at_Len[,,2,i,1] <- exp(smooth_imm)[,,2,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    # mat_N_at_Len[,,1,i,1] <- exp(smooth_mat)[,,1,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))
    # mat_N_at_Len[,,2,i,1] <- exp(smooth_mat)[,,2,i] # matrix(rnorm(n = length(lat) * length(lon),mean = 1,sd = 1),nrow = length(lat), ncol = length(lon))

    model <- RMexp(scale = 20,var = 1e10)
    imm_N_at_Len[,,1,i,1] <- (as.matrix(RFsimulate(model, lon, rev(lat), grid=TRUE)))
    imm_N_at_Len[,,2,i,1] <- (as.matrix(RFsimulate(model, lon, rev(lat), grid=TRUE)))
    mat_N_at_Len[,,1,i,1] <- (as.matrix(RFsimulate(model, lon, rev(lat), grid=TRUE)))
    mat_N_at_Len[,,2,i,1] <- (as.matrix(RFsimulate(model, lon, rev(lat), grid=TRUE)))

  }
  
  
}

imm_N_at_Len[imm_N_at_Len=="NaN"]<-0
mat_N_at_Len[mat_N_at_Len=="NaN"]<-0

imm_N_at_Len[imm_N_at_Len<0]<-0
mat_N_at_Len[mat_N_at_Len<0]<-0


## Abundance
load("4_full_MSE/data/Abundance_Recruitment_Assessment.Rdata") # Outputs provided by Mathieu from GMACS
init_year = 37 # Time series: 1982-2022
Ab_males = unlist(Snow_Out$Abundance$N_males[init_year,])
Ab_males_mat = unlist(Snow_Out$Abundance$N_males_mature[init_year,])
Ab_females = unlist(Snow_Out$Abundance$N_females[init_year,])
Ab_females_mat = unlist(Snow_Out$Abundance$N_females_mature[init_year,])

#==This makes up a random distribution for the population
#==only used for testing
fake_dist_data<-1
if(fake_dist_data==1)
{
  #==set dummy initial distribution--take this from the survey for real
  #==this is only for getting mechanics down
  imm_N_at_Len[,,1,1,1]<-init_juv * (Ab_females[1] - Ab_females_mat[1])
  imm_N_at_Len[,,2,1,1]<-init_juv * (Ab_males[1] - Ab_males_mat[1])
  mat_N_at_Len[,,1,1,1]<-init_juv * Ab_females_mat[1]
  mat_N_at_Len[,,2,1,1]<-init_juv * Ab_males_mat[1]
  
  for(x in 2:length(sizes))
  {
    
    if(sizes[x] < ad_size){
      imm_N_at_Len[,,1,x,1]<-init_juv * (Ab_females[x] - Ab_females_mat[x])
      imm_N_at_Len[,,2,x,1]<-init_juv * (Ab_males[x] - Ab_males_mat[x])
      mat_N_at_Len[,,1,x,1]<-init_juv * Ab_females_mat[x]
      mat_N_at_Len[,,2,x,1]<-init_juv * Ab_males_mat[x]
    }else if(sizes[x] > 45){
      imm_N_at_Len[,,1,x,1]<-init_adult * (Ab_females[x] - Ab_females_mat[x])
      imm_N_at_Len[,,2,x,1]<-init_adult * (Ab_males[x] - Ab_males_mat[x])
      mat_N_at_Len[,,1,x,1]<-init_adult * Ab_females_mat[x]
      mat_N_at_Len[,,2,x,1]<-init_adult * Ab_males_mat[x]
    }
    
  }
  
}

#==delete critters where there is land (this will be used for movement as well)

if(size_class_settings == "rough") area_mask = "EBS_only"

land_mask<-matrix(1,ncol=length(lon),nrow=length(lat),byrow=T)
land_matrix = t(matrix(land,ncol=length(lon),nrow=length(lat)))
EBS_mask_matrix = t(matrix(EBS_mask,ncol=length(lon),nrow=length(lat)))
EBS_mask_matrix[which(EBS_mask_matrix == T)] = 1
EBS_mask_matrix[which(EBS_mask_matrix == F)] = 0

# plot(t(land_matrix))
# plot(t(EBS_mask_matrix))
for(x in 1:length(lat))
  for(y in 1:length(lon))
  {
    
    if(land_matrix[x,y] == 1){
      land_mask[x,y]<-0
    }
    
  }

#==ugh. this is dumb, but how to automate?
g<-function(m) t(m)[,nrow(m):1]

land_mask[25:40,20:40]
land_mask[31,32]<-0
land_mask[32,29]<-0
# filled.contour(land_mask)

# filled.contour(x=lon,y=rev(lat),g(imm_N_at_Len[,,1,1,1]))
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
# filled.contour(x=lon,y=rev(lat),g(imm_N_at_Len[,,1,1,1]))

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
# Parameterize with GMACS
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

# Type of model
growth_model <- "cody_model" # "max_model" "cody_model"
# Cody's model --> non-spatial life-history parameters
# Maxime's model --> spatially varying life-history parameters

growth_param_est = growth_param[1:6,] %>% 
  dplyr::select(Estimate) %>% 
  unlist

if(growth_model == "cody_model"){
  
  #==Cody's growth parameterization for generation of non spatially varying growth parameters
  alpha_grow_f_imm<-4 # growth_param_est[4]
  alpha_grow_m_imm<-7 # growth_param_est[1]
  beta_grow_f_imm<-1.05 # growth_param_est[5]
  beta_grow_m_imm<-1.1 # growth_param_est[2]
  
  alpha_grow_f_mat<-4 # growth_param_est[4]
  alpha_grow_m_mat<-7 # growth_param_est[1]
  beta_grow_f_mat<-1.05 # growth_param_est[5]
  beta_grow_m_mat<-1.1 # growth_param_est[2]
  
  growth_sd_imm<-c(5,4) # Where to find better values?
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


## Recruitment
#-------------
recruit_female = unlist(Snow_Out$Recruitment$recruits[1,])
meanlog_recruit_female = mean(log(recruit_female))
sdlog_female = sd(log(recruit_female))

recruit_male = unlist(Snow_Out$Recruitment$recruits[2,])
meanlog_recruit_male = mean(log(recruit_male))
sdlog_male = sd(log(recruit_male))


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


## Movement
#----------
## Movement at size
move_len_50<-60
move_len_95<-70

## Movement matrices
#-------------------
# --> might need to make different matrices for mature and immature
compute_movement_matrix = T

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
  
  # plot(diffusion_gg[1:100,1:100], breaks=20)
  # hist(diffusion_gg)
  
  ## Taxis for juveniles
  #---------------------
  taxis_coef_juv = 10^4 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations
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
  # plot(diffusion_gg[1:100,1:100],breaks=20)
  if( any(mrate_gg_juv-diag(diag(mrate_gg_juv))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
  mfraction_gg_juv = Matrix::expm(mrate_gg_juv)
  # test = matrix(mfraction_gg_juv@x,nrow=mfraction_gg_juv@Dim[1],ncol=mfraction_gg_juv@Dim[2])
  # plot(test[1:100,1:100])
  stationary_g_juv = eigen(mfraction_gg_juv)$vectors[,1]
  stationary_g_juv = stationary_g_juv / sum(stationary_g_juv)
  test = matrix(stationary_g_juv,nrow=40,ncol=40,byrow = T)
  plot(test)

  ## Taxis for adults
  #------------------
  taxis_coef_ad = 10^4 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations
  preference_g_ad = (init_adult - mean(init_adult)) / sd(init_adult)
  preference_g_ad = as.vector(t(preference_g_ad))
  # # check
  # test = matrix(preference_g,nrow = 40,ncol = 40,byrow = T)
  # plot(test)
  
  taxis_gg_ad = A_gg_mat *  taxis_coef_ad * DeltaT / A * exp(outer(preference_g_ad, preference_g_ad, "-") * DeltaT / sqrt(A))
  diag(taxis_gg_ad) = -1 * colSums(taxis_gg_ad)
  
  # Total
  mrate_gg_ad = taxis_gg_ad # + diffusion_gg
  if( any(mrate_gg_ad-diag(diag(mrate_gg_ad))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
  
  mfraction_gg_ad = Matrix::expm(mrate_gg_ad)
  # test = matrix(mfraction_gg@x,nrow=mfraction_gg@Dim[1],ncol=mfraction_gg@Dim[2])
  # plot(test[1:100,1:100])
  
  stationary_g_ad = eigen(mfraction_gg_ad)$vectors[,1]
  stationary_g_ad = stationary_g_ad / sum(stationary_g_ad)
  # test = matrix(stationary_g_ad,nrow=40,ncol=40,byrow = T)
  # plot(test)
  
}

## Fishery selectivity
#---------------------
# Parameterize with GMACS
selec_curve = repfile$selectivity$Start_Y %>% 
  filter(fleet == 1) %>% 
  dplyr::select_at(vars(starts_with("SizeC_")))

fish_sel = selec_curve

# fish_sel_50_f<-NA
# fish_sel_95_f<-NA
# fish_sel_50_m<-50
# fish_sel_95_m<-101
# 
# fish_sel<-rbind(1/(1+exp(-log(19)*(sizes-fish_sel_50_f)/(fish_sel_95_f-fish_sel_50_f))),
#                 1/(1+exp(-log(19)*(sizes-fish_sel_50_m)/(fish_sel_95_m-fish_sel_50_m))))
# fish_sel[is.na(fish_sel)]<-0

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
# filled.contour(x=lon,y=rev(lat),g(distance_map*land_mask),plot.axes=c(map(add=TRUE,fill=T,col='grey'),
#                                                                       points(y=port_lat,x=port_lon,pch=16,col='red')))
#write.csv(distance_map,'dist.csv')


#================================================
# calculate costs to fish
#=============================================
#==this should be related to the amount of fish in a patch
cost_fish <- 0

cost_travel <- 0
cost_patch <- cost_travel * distance_map + cost_fish
price <- 1.5

fishers<-20
quota<-rep(1000,fishers)

fishing_process="stochastic"
# fishing_process="max_benefit_min_dist"

#=============================================
# PROJJEEEECCCT
#============================================
list_it = 0

# Abundance
imm_N_at_Len_full = list()
mat_N_at_Len_full = list()

#  Exploitation
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
for(cost_travel in 1000){ # c(0,1000,1000*2)
  
  total_spatial_catch<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
  catch_by_fisher<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period),fishers))
  profit_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  cost_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  all_net_benefit_patch=array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
  all_chosen_patch=array(0,dim=c(2,length(proj_period),fishers))
  
  #==indices: lat,lon,sex,size,time
  for(t in 1:(length(proj_period)-1))
    #for(t in 1:320)
  {
    print(t)
    #==create a 'working' array for a given time step of both mature and immature critters
    
    temp_imm_N<-imm_N_at_Len[,,,,t]
    temp_mat_N<-mat_N_at_Len[,,,,t]
    # filled.contour(x=lon,y=rev(lat),g(temp_mat_N[,,1,5]),plot.axes=map(add=TRUE,fill=T,col='grey') )
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
          chosen_patch=chosen_patch[1,]
          
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
            chosen_patch=chosen_patch[1,]
            
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
          # par(mfrow = c(3,2))
          # 
          # mov_imm_temp = temp_imm_N[,,sex,x]
          # mov_mat_temp = temp_mat_N[,,sex,x]
          # 
          # plot(mov_imm_temp,main=paste0("t = 0"),asp = 1)
          # 
          # for(t in 1:5){
          # 
          #   mov_imm_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp)) # mrate_gg %*%
          #   mov_mat_v = mfraction_gg_ad %*% as.vector(t(mov_mat_temp))  # mrate_gg %*%
          # 
          #   mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
          #   mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
          # 
          #   print(which((mov_imm_temp != mov_mat_temp_2)))
          # 
          #   mov_imm_temp = mov_imm_temp_2
          #   mov_mat_temp = mov_mat_temp_2
          # 
          #   plot(mov_imm_temp,main=paste0("t = ",t),asp = 1)
          # 
          # }
          # 
          ## Or
          # 
          # land_mask_na = land_mask
          # land_mask_na[which(land_mask_na == 0)] = NA
          # for(t in 1:(12*3)){
          #   if(t %in% c(1,13,25)){
          #     x11()
          #     par(mfrow = c(4,3))
          #   }
          # 
          #   plot(mat_N_at_Len[,,2,7,t] * land_mask_na,main=paste0(t),asp = 1)
          # }

        }
      
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
      aggreg_rec_1 = rlnorm(n = 1, mean = meanlog_recruit_female, sdlog = sdlog_female)
      aggreg_rec_2 = rlnorm(n = 1, mean = meanlog_recruit_male, sdlog = sdlog_male)
      
      tmp_rec_1 <- init_juv * aggreg_rec_1
      # matrix(rep(aggreg_rec_1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,1,1,1]))
      tmp_rec_1[tmp_rec_1<0]<-0
      tmp_rec_2 <- init_juv * aggreg_rec_2
      # matrix(rnorm(length(imm_N_at_Len[,,2,1,1]),1,1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,2,1,1]))
      tmp_rec_2[tmp_rec_2<0]<-0
      
      for(r in 1:rec_sizes)
      {
        temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + imm_N_at_Len[,,1,r,1]*tmp_rec_1*prop_rec[r]
        temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + imm_N_at_Len[,,2,r,1]*tmp_rec_2*prop_rec[r]
      }
    }
    
    #==update dynamics
    imm_N_at_Len[,,1,,t+1] <-  temp_imm_N[,,1,]*exp(-imm_fem_M*1/year_step)
    imm_N_at_Len[,,2,,t+1] <-  temp_imm_N[,,2,]*exp(-imm_male_M*1/year_step)
    mat_N_at_Len[,,1,,t+1] <-  temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
    mat_N_at_Len[,,2,,t+1] <-  temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
    
    #==generate scientific (1 sample per cell grid - could be something else)
    # Data_Geostat
    
    ###############################
    ## Pass to do.call to go faster
    ###############################
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
      
      if(print_messages) print("Make Stock assessment")
      
      if(SA == "spatialIPM"){
        
        if(display_print) print("spatial IPM")
        
        source(paste0(project_spatialIPM,"03_spatial_model/run_model_mse.R"))
        
      }
      
      if(SA == "nonspatialIPM"){
        
      }
      
      if(SA == "GMACS"){
        
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
