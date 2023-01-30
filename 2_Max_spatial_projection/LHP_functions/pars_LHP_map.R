
pars_LHP_map <- function(pref_hab,
                         lon,
                         lat,
                         clim_sc,
                         years,
                         pars_pref_hab,
                         x_omega,
                         x_epsilon,
                         pars_LHP,
                         scale_g,
                         LHP,
                         ad_eff=TRUE,
                         plot=FALSE)
{
  
  suppressWarnings({ suppressMessages({ 


  # -------------------------------------------------------------------------
  # OMEGA
  # -------------------------------------------------------------------------
  
  # 1 - Simulate Field conditionned on climate scenario----------------------
  # -------------------------------------------------------------------------
  
  # Load Temperature scenarios  ---------------------------------------------
  load("2_Max_spatial_projection/Climate_scenarios/mn_var_all.RData")
  En.cond <- mn_var_all %>% filter(simulation %in% clim_sc , year %in% years ) %>% dplyr::select(year,simulation, latitude,longitude, val)
  En.cond <- En.cond %>% mutate(Lon = ifelse(longitude <= 180, longitude, longitude -360),Lat=latitude)
  dat_cond <-  data.frame(coords=En.cond[,c("Lon","Lat")],Year=En.cond$year, Temp=En.cond$val, Climate_sc = En.cond$simulation)
  colnames(dat_cond) <- c( "Lon","Lat","Year","Temp","Climate_sc")
  data(wrld_simpl)
  crs_LL = CRS(proj4string(wrld_simpl))
  dat_cond_sf <- st_as_sf(dat_cond, coords=c("Lon","Lat"), crs=crs_LL)
  
  
  # Extract year and scenarios ----------------------------------------------
  temp_sc_year_sf <- NULL
  En.cov_plot <- NULL
  temp_sc_year_sf <- (dat_cond_sf) %>% dplyr::filter(Year == years,Climate_sc == clim_sc)
  dat_cond_temp <- data.frame(coords=(st_coordinates(temp_sc_year_sf)),Temp=temp_sc_year_sf$Temp)
  
  # Simulate Field------------------------------------------------------------
  spde_cov <-  RMexp() #RMgauss()#var=0.1*0.1,scale=scale_g
  En.cov_temp <- RFsimulate(model = spde_cov, x=loc_x[,2], y=loc_x[,3],data= dat_cond_temp)
  En.cov_plot <- data.frame( Temp = En.cov_temp$Temp,Lon=loc_x[,2],Lat=loc_x[,3])
  
  if (plot==TRUE){
  world_sf <- st_as_sf(wrld_simpl,crs=CRS(proj4string(wrld_simpl)))
  p_Fcond <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))+
    geom_raster(data=En.cov_plot,aes(Lon, Lat,fill=Temp)) +
    scale_fill_viridis()  +
    theme_bw() + ylab("") +xlab("")

  
  ggsave((paste0(DateFile, '/Spatial_Field_env',LHP,'.png')),plot=p_Fcond,
         width = 27,
         height = 18,
         units = "cm")
  
  }
  
  # 2-Apply Preferential Habitat function to Field---------------------------
  # -------------------------------------------------------------------------
  X_s <- NULL
  
  X_s_temp <- NULL
  rast <- raster(ncol=150,nrow=150)
  loc_x_tmp <- loc_x[,c(2,3)]
  colnames(loc_x_tmp) <- c("x" ,"y")
  extent(rast) <- extent(as.data.frame(loc_x_tmp))
  
  P2 = SpatialPoints(loc_x_tmp)
  data_temp <-   En.cov_plot 
  P2$data =     data_temp$Temp
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  
  envir_stack <- stack(rast.temp)
  names(envir_stack) <- c('temp')
  
  # Apply the habitat preference function
  if (pref_hab=="logi") fun_hab <- c(fun='logisticFun', alpha=pars_pref_hab[1], beta=pars_pref_hab[2])
  if (pref_hab=="quadr") fun_hab <-  c(fun="dnorm",mean=pars_pref_hab[1],sd=pars_pref_hab[2])
  if (pref_hab=="line") fun_hab <-  c(fun="linearFun",a=pars_pref_hab[1],b=pars_pref_hab[2])
  
  parameters <- virtualspecies::formatFunctions(temp = fun_hab)
  envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)
  
  X_s <- data.frame(Lat=loc_x[,3],Lon=loc_x[,2], Temp_pref =   raster::extract(envirosuitability$suitab.raster,loc_x[,c(2,3)]) )   
  
  if (plot==TRUE){
  # Plot suitability function    
  png(paste('PH_fun',pref_hab,'_',LHP,'.png',sep=''), height = 7, width = 8, units = 'in', res=600)
  virtualspecies::plotResponse(envirosuitability) #plot response curves
  dev.off()
  }
  
  # 3- Rescale X_s to get Omega varying within a range ----------------------
  # -------------------------------------------------------------------------
  
  par_omega <-   rescale(X_s$Temp_pref, to = c(-x_omega*pars_LHP[1],    x_omega*pars_LHP[1]))
  
  if (plot==TRUE){
  # Plot Omega
  par_omega_plot <- data.frame( Omega = par_omega,lon=loc_x[,2],lat=loc_x[,3])
  
  p_omega <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))+
    geom_raster(data=par_omega_plot,aes(lon, lat,fill=Omega)) +
    scale_fill_viridis()  +
    theme_bw() + ylab("") +xlab("")
  
  
  ggsave((paste0(DateFile, '/Omega',LHP,'.png')),plot=p_omega,
         width = 27,
         height = 18,
         units = "cm")
  }
  
  # -------------------------------------------------------------------------
  #  EPSILON ----------------------------------------------------------------
  # -------------------------------------------------------------------------
  
  # Define variance and scale for RF
  # define model
  model_epsilon_par    <- RMgauss(var=pars_LHP[2]^2, scale=scale_g)

  # Simulate RF
  par_epsilon_RF = RFsimulate(model = model_epsilon_par, x=loc_x[,2], y=loc_x[,3])@data[,1]

  # Rescale Epsilon 
  par_epsilon <- rescale(par_epsilon_RF, to = c(-x_epsilon*pars_LHP[1],    x_epsilon*pars_LHP[1]))

  
  if (plot==TRUE){
  # Plot Epsilon
  par_epsilon_plot <- data.frame( Epsilon = par_epsilon,lon=loc_x[,2],lat=loc_x[,3])
  p_epsilon <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))+
    geom_raster(data=par_epsilon_plot,aes(lon, lat,fill=Epsilon)) +
    scale_fill_viridis()  +
    theme_bw() + ylab("") +xlab("")
  
  
  ggsave((paste0(DateFile, '/Epsilon',LHP,'.png')),plot=p_epsilon,
         width = 27,
         height = 18,
         units = "cm")
  }
  

 # -------------------------------------------------------------------------
 # Generate parameters ---------------------------------------------------
 # -------------------------------------------------------------------------
  if (ad_eff==TRUE){
  LHP_pars_map = pars_LHP[1] + par_epsilon +  par_omega
  }else{
  LHP_pars_map = exp(pars_LHP[1] + par_epsilon +  par_omega) 
  }
  
  if (plot==TRUE){
  LHP_pars_map_plot <- data.frame( LHP = LHP_pars_map,lon=loc_x[,2],lat=loc_x[,3])
  plot <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))+
    geom_raster(data=LHP_pars_map_plot,aes(lon, lat,fill=LHP)) +
    scale_fill_viridis()  +
    theme_bw() + ylab("") +xlab("")
  }
  
  LHP_pars <- list("LHP_pars_map"=LHP_pars_map,"par_epsilon"=par_epsilon, "par_omega"=par_omega)
  return(LHP_pars)
  })})
  
}

