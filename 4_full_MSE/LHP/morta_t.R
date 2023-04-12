######################
## Modelling mortality
######################

mortality <- function(E_morta_par,
                      pars_LHP_setting,
                      n_s,n_p,
                      G_spatial=F,
                      plot=FALSE){
  
  #-- Map of growth matrix
  morta_matrix = array (data = 0, dim = c(n_s))
  morta_par = array (data = NA, dim = c(n_s,pars_LHP_setting$n_grpar))
  
  #-- no space
  if (!pars_LHP_setting$LHP_spatial) {
    for (k in 1:n_s) {
      morta_matrix[k] = E_morta_par
    }
    
  }else{
    
    # -- Space  
    source("2_Max_spatial_projection/LHP_functions/pars_LHP_map.R")
    
    growth_par <- array(0, c(n_s,pars_LHP_setting$n_grpar))
    par <- pars_LHP_map(
      pars_LHP_setting$pref_hab,
      pars_LHP_setting$lon,
      pars_LHP_setting$lat,
      pars_LHP_setting$clim_sc,
      pars_LHP_setting$year_LHP,
      pars_LHP_setting$pars_pref_hab,
      pars_LHP_setting$x_omega,
      pars_LHP_setting$x_epsilon,
      pars_LHP_setting$pars_LHP[, i],
      pars_LHP_setting$scale_g,
      pars_LHP_setting$LHP,
      pars_LHP_setting$ad_eff,
      plot=FALSE
    )
    
    morta_matrix =  par$LHP_pars_map
    
  }
  
  # match lat/lon with cell number
  Morta_temp <- melt(morta_matrix)
  colnames(Morta_temp) <- c("cell","morta")
  Morta_temp <- data.frame(Morta_temp)
  
  Morta_temp$Size_class_from <- as.factor(Morta_temp$Size_class_from)
  Morta_temp$Size_class_to <- as.factor(Morta_temp$Size_class_to)
  
  Cell <- data.frame(cell = loc_x[,1], lat=loc_x[,3],lon=loc_x[,2])
  Morta <- left_join(Morta_temp,Cell,by = "cell")
  
  # Transform to appropriate format for OM input
  lon2 <- data.frame(lon=lon,lon_nb = c(1:length(lon)))
  lat2 <- data.frame(lat=lat,lat_nb = c(1:length(lat)))
  colnames(lon2) <- c("lon","cell")
  colnames(lat2) <- c("lat","cell")
  Morta_temp2 <- left_join(Morta, lon2,by = c("cell", "lon"))
  Morta_temp3 <- left_join(Morta_temp2, lat2,by = c("cell", "lat"))
  
  if(plot==TRUE){
    
    # Plot
    if(!G_spatial) {p_size_trans <- p }else{ p_size_trans <- p + facet_wrap(~cell)}
    
    p <- ggplot()+
      geom_raster(Morta_temp2, mapping=aes(lon, lat, fill= morta))+
      scale_fill_viridis()+
      theme_bw()+ 
      geom_sf(data = world_sf,fill="black",color=NA)+ 
      coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))
    
    plot(p)
    
  }
  
  # Format for pop. dyn model
  temp <- Morta_temp2 %>% mutate(value=morta) %>% dplyr::select(lat,lon,morta)
  Morta_OM <- acast(temp, lat~lon)
  Morta_OM <- ifelse(is.na(Morta_OM),0,Morta_OM)
  
  return(Morta_OM)
  
}


