#################
## Fleet dynamics
#################
# T. Wang, B alglave

library(ggplot2)
library(ggridges)
library(sf)
library(tidyverse)

#==== read in model outputs
load(file="4_full_MSE/data/fit_mnl_fleet_only_r_none.rda")

#==== read in the data
the_d = read.csv("4_full_MSE/data/cpue_alt_covariate_prewide_2023_0329.csv")

#==== get parameter estimates and se
parameters_fixed =  fit$par[which(names(fit$par) %in% c("CpueV_par","CpueF_par","Ice_par","Trad_par",
                                                 "Sort_par","Cong_par","Port_par","Surv_par",
                                                 "CvV_par",
                                                 "CvF_par"))]
# se=summary(sd_error_rep,"fixed")[,2]
# se_fixed = se[which(names(se) %in% c("CpueV_par","CpueF_par","Ice_par","Trad_par",
#                                      "Sort_par","Cong_par","Port_par","Surv_par",
#                                      "CvV_par",
#                                      "CvF_par"))]

#==== get area intercepts
AreaPar = fit$par[which(names(fit$par) %in% c("AreaPar"))]
AreaPar=append(AreaPar,0)

area_par_df = bind_cols(stat_area = sort(unique(the_d$stat_area)),AreaPar=AreaPar)
area_id = sort(unique(the_d$stat_area))
N_area = length(area_id)

#==== Get average of non-effort interacting covariates per stat_area (Sort, Port, Ice)
#== Ice
# how much fishing occurs each month?
# the_d %>% group_by(month) %>%
#   summarise(effort=sum(effort,na.rm = T)) %>%
#   mutate(prop_eff = effort/sum(effort))
# focus only on months 1,2,and 3

ice_area_grid <- read.csv("4_full_MSE/data/weekly_ice_area_per_grid_cell.csv") %>% dplyr::select(-X)
ice_area_grid <- ice_area_grid %>% rename(ice_date = date) %>% mutate(ice_date = as.Date(ice_date))

mean_ice = ice_area_grid %>% 
  filter(month %in% 1:3 & year %in% 2006:2020) %>%
  group_by(stat_area,month) %>%
  summarise(ice_prop_mean = mean(ice_prop, na.rm=T)) %>% filter(stat_area %in% unique(the_d$stat_area))%>%
 arrange(stat_area)

#== Sorting
sorting_season = read.csv("4_full_MSE/data/sorting_by_season.csv") %>%
  dplyr::select(-X) %>%
  dplyr::select(statarea,season,legal_prop)

mean_sort = sorting_season %>%
  group_by(statarea) %>%
  summarise(legal_prop_mean = mean(legal_prop, na.rm=T)) %>%
  filter(statarea %in% unique(the_d$stat_area))

#== Port distance
# since port of landing varies for each trip, I assume that vessels will deliver to nearest port
# issue: this won't deal with catcher-processor vessels and stationary floating processors
plants = sf::st_read("4_full_MSE/data/ebs_crab_plants/ebs_crab_plants.shp") %>% 
  st_drop_geometry() %>% dplyr::select(General_Ar,Latitude,Longitude) 

the_d_plants = the_d %>% dplyr::select(stat_area, lat,lon) %>% expand_grid(plants)

stat_area_df = the_d_plants %>% 
  group_by(stat_area,lon,lat) %>% 
  dplyr::slice(1) %>% 
  dplyr::select(stat_area,lon,lat)

library(geosphere)
the_d_plants$port_dist = distHaversine(the_d_plants[, c('lon', 'lat')], the_d_plants[, c('Longitude', 'Latitude')])/1000

stat_area_port_dist = the_d_plants %>% 
  group_by(stat_area) %>% 
  filter(port_dist==min(port_dist)) %>% 
  dplyr::select(stat_area,port_dist) %>% unique()

# set up months and how many fishing trips per month
month_seq=rep(1:3, each=2)
month_i_df = data.frame(run_i=1:length(month_seq),month=month_seq)

# set up data frames that keep track of fishing effort, CPUE, and coeff of variance of CPUE
df = data.frame(matrix(-99,nrow=length(month_seq),ncol=N_area)) %>% 
  magrittr::set_colnames(area_id) 

fishery_history <- replicate(3, df, simplify = FALSE)
names(fishery_history) <- c("effort","cpue","risk")


#==== initialize first fishing spatial distribution
for(i in 1:nrow(month_i_df)){
  area_util_ti = area_par_df$AreaPar+stat_area_port_dist$port_dist*parameters_fixed["Port_par"]+
    # ice cover varies per month
    mean_ice$ice_prop_mean[mean_ice$month == month_i_df[i,"month"]]*parameters_fixed["Ice_par"]+
    mean_sort$legal_prop_mean*parameters_fixed["Sort_par"]
  if (i>1){
   # time varying covariates for time steps after 1
  area_util_ti = area_util_ti + fishery_history$effort[i-1,]*parameters_fixed["Trad_par"]+
    fishery_history$cpue[i-1,]*parameters_fixed["CpueF_par"]+
    fishery_history$risk[i-1,]*parameters_fixed["CvF_par"]
  }
  
  # proportion of fishing each stat_area
  prob_ti = exp(area_util_ti)/sum(exp(area_util_ti))

  fishery_history$effort[i,] = prob_ti
  
  fishery_history$cpue[i,] = rexp(32, rate = 1/100)
  fishery_history$risk[i,] = rnorm(32,10,2)
  
}

stat_area_df$stat_area = as.character(stat_area_df$stat_area)

effort_df = data.frame(stat_area = colnames(fishery_history$effort),
                      effort = t(fishery_history$effort)) %>% 
  inner_join(stat_area_df) %>% 
  pivot_longer(effort.1:effort.6)

ggplot(effort_df)+
  geom_raster(aes(x=lon,y=lat,fill=value))+
  scale_fill_distiller(palette = "Spectral")+
  facet_wrap(.~name)

