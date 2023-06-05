## Initialization for each year
#------------------------------

## Temperature data frame
En.cond <- mn_var_all %>% filter(simulation %in% clim_sc, year %in% Years_climsc[t] ) %>% dplyr::select(year,simulation, latitude,longitude, val)
En.cond <- En.cond %>% mutate(Lon = ifelse(longitude <= 180, longitude, longitude -360),Lat=latitude)
dat_cond <-  data.frame(coords=En.cond[,c("Lon","Lat")],Year=En.cond$year, Temp=En.cond$val, Climate_sc = En.cond$simulation)
colnames(dat_cond) <- c( "Lon","Lat","Year","Temp","Climate_sc")
crs_LL = CRS(proj4string(wrld_simpl))
dat_cond_sf <- st_as_sf(dat_cond, coords=c("Lon","Lat"), crs=crs_LL)

## Cold pool extent
if(! Years_climsc[t] > 2021) proj_cold_pool_df_2 = hind_cold_pool_df %>% 
  filter(str_detect(GCM_scen,"gfdl") & str_detect(GCM_scen,clim_sc))

if(Years_climsc[t] > 2021) proj_cold_pool_df_2 = proj_cold_pool_df %>%
  filter(simulation %in% clim_sc)

# values for rescaling cold pool for GAM
mu_forc = mean(proj_cold_pool_df_2$fracbelow2)
sd_forc = sd(proj_cold_pool_df_2$fracbelow2)
mu_hind = mean(m_gam_dat$coldpool,na.rm=T)
sd_hind = sd(m_gam_dat$coldpool,na.rm=T)

proj_cold_pool_df_3 = proj_cold_pool_df_2 %>% 
  filter(year %in% Years_climsc[t]) %>% 
  mutate(fracbelow2_scaled = (fracbelow2 - mu_forc) * sd_hind / sd_forc + mu_hind)



if(growth_model == "max_model"){
  
  print("init growth")
  source("4_full_MSE/source/projection/growth_t.R")
  
}

if(morta_model == "max_model"){
  
  print("init morta")
  source("4_full_MSE/source/projection/morta_t.R")
  
}

if(recruit_model == "max_model"){
  
  print("init recruit")
  source("4_full_MSE/source/projection/recruit_t.R")
  
}

