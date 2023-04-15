########################
## Modelling recruitment
########################

recruitment <- function(E_recruit_par,
                        pars_LHP_setting,
                        n_s,
                        R_spatial=F,
                        plot=FALSE){
  
  #-- Map of growth matrix
  recruit_matrix = array (data = 0, dim = c(n_s))
  recruit_par = array (data = NA, dim = c(n_s,pars_LHP_setting$n_grpar))
  
  #-- no space
  if (!pars_LHP_setting$LHP_spatial) {
    for (k in 1:n_s) {
      recruit_matrix[k] = E_recruit_par
    }
    
  }else{
    
    # -- Space  
    source("2_Max_spatial_projection/LHP_functions/pars_LHP_map.R")
    
    recruit_par <- array(0, c(n_s,pars_LHP_setting$n_grpar))
    par <- pars_LHP_map(
      pars_LHP_setting$pref_hab,
      pars_LHP_setting$lon,
      pars_LHP_setting$lat,
      pars_LHP_setting$clim_sc,
      pars_LHP_setting$year_LHP,
      pars_LHP_setting$pars_pref_hab,
      pars_LHP_setting$x_omega,
      pars_LHP_setting$x_epsilon,
      pars_LHP_setting$pars_LHP,
      pars_LHP_setting$scale_g,
      pars_LHP_setting$LHP,
      pars_LHP_setting$ad_eff,
      plot=FALSE
    )
    
    recruit_matrix =  par$LHP_pars_map
    
  }
  
  # match lat/lon with cell number
  Recruit_temp <- data.frame(cell = 1:n_s,
                             recruit=recruit_matrix)
  
  Cell <- data.frame(cell = loc_x[,1], lat=loc_x[,3],lon=loc_x[,2])
  Recruit <- left_join(Recruit_temp,Cell,by = "cell")
  
  # Transform to appropriate format for OM input
  lon2 <- data.frame(lon=lon,lon_nb = c(1:length(lon)))
  lat2 <- data.frame(lat=lat,lat_nb = c(1:length(lat)))
  colnames(lon2) <- c("lon","cell")
  colnames(lat2) <- c("lat","cell")
  Recruit_temp2 <- left_join(Recruit, lon2,by = c("cell", "lon"))
  Recruit_temp3 <- left_join(Recruit_temp2, lat2,by = c("cell", "lat"))
  
  if(plot==TRUE){
    
    # Plot
    if(!G_spatial) {p_size_trans <- p }else{ p_size_trans <- p + facet_wrap(~cell)}
    
    p <- ggplot()+
      geom_raster(Recruit_temp2, mapping=aes(lon, lat, fill= recruit))+
      scale_fill_viridis()+
      theme_bw()+ 
      geom_sf(data = world_sf,fill="black",color=NA)+ 
      coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))
    
    plot(p)
    
  }
  
  # Format for pop. dyn model
  temp <- Recruit_temp2 %>% 
    mutate(value=recruit) %>% 
    dplyr::select(lat,lon,recruit) %>% 
    arrange(desc(lat))
  
  Recruit_OM <- t(g(acast(temp, lat~lon,value.var = "recruit")))
  Recruit_OM <- ifelse(is.na(Recruit_OM),0,Recruit_OM)
  
  # plot(t(g(acast(temp, lat~lon))))
  
  return(Recruit_OM)
  
}

