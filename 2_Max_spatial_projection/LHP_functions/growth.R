
# pars_LHP_setting = pars_Growth_setting_f_mat

growth <- function(sizes,
                   binclass,
                   pars_LHP_setting,
                   n_s,n_p,
                   G_spatial=F,
                   plot=FALSE
                   ){

#-- Map of growth matrix
growth_matrix = array (data = 0, dim = c(n_s,n_p,n_p))
growth_par = array (data = NA, dim = c(n_s,pars_LHP_setting$n_grpar))

#-- no space
source("2_Max_spatial_projection/LHP_functions/script_growth_inc.R")

if (!pars_LHP_setting$LHP_spatial) {
  for (k in 1:n_s) {
    growth_matrix[k, ,] = growth_trans(sizes,
                                       pars_LHP_setting$pars_LHP[1,2],
                                       pars_LHP_setting$pars_LHP[1,1],
                                       sizes[1],
                                       sizes[n_p],
                                       binclass,
                                       pars_LHP_setting$pars_LHP[1,3] )
  }
  
}else{

# -- Space  
source("2_Max_spatial_projection/LHP_functions/pars_LHP_map.R")  
  
  growth_par <- array(0, c(n_s,pars_LHP_setting$n_grpar))
  for (i in 1:(pars_LHP_setting$n_grpar - 1)) {
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
    
   growth_par[, i] =  par$LHP_pars_map

  }
  
  growth_par[, 3] =  pars_LHP_setting$pars_LHP[1, 3]
  
  for (k in 1:n_s) {
    growth_matrix[k, , ] = growth_trans(sizes,
                                        growth_par[k, 1],
                                        growth_par[k, 2],
                                        sizes[1],
                                        sizes[n_p],
                                        binclass,
                                        growth_par[k, 3])
  }
}

# match lat/lon with cell number
Growth_temp <- melt(growth_matrix)
colnames(Growth_temp) <- c("cell","Size_class_from","Size_class_to", "prob")
Growth_temp <- data.frame(Growth_temp)

Growth_temp$Size_class_from <- as.factor(Growth_temp$Size_class_from)
Growth_temp$Size_class_to <- as.factor(Growth_temp$Size_class_to)

Cell <- data.frame(cell = loc_x[,1], lat=loc_x[,3],lon=loc_x[,2])
Growth <- left_join(Growth_temp,Cell,by = "cell")

# Transform to appropriate format for OM input
lon2 <- data.frame(lon=lon,lon_nb = c(1:length(lon)))
lat2 <- data.frame(lat=lat,lat_nb = c(1:length(lat)))
colnames(lon2) <- c("lon","cell")
colnames(lat2) <- c("lat","cell")
Growth_temp2 <- left_join(Growth, lon2,by = c("cell", "lon"))
Growth_temp3 <- left_join(Growth_temp2, lat2,by = c("cell", "lat"))

if(plot==TRUE){
  
  
  
  # Plot
  p <- ggplot(Growth_temp %>% filter(cell %in% c(seq(1,1000,100))), aes(Size_class_from, Size_class_to, fill= prob)) + 
    geom_tile()+scale_fill_viridis() + theme_bw()
  
  if (!G_spatial) {p_size_trans <- p }else{ p_size_trans <- p + facet_wrap(~cell)}
  
  # ggsave((paste0(DateFile, '/Growth_size_trans.png')),plot=p_size_trans,
  #        width = 27,
  #        height = 18,
  #        units = "cm")
  
  test = Growth_temp2$lon %>% 
    filter(Size_class_from %in% 1:2 & Size_class_to %in% 3:5) %>% 
    mutate()
  
  p <- ggplot() +   geom_raster(test, mapping=aes(lon, lat, fill= prob))+
    scale_fill_viridis() +  
    theme_classic()+ 
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())+
    xlab("")+ylab("")+
    geom_sf(data = world_sf,fill="black",color=NA)+
    coord_sf(xlim=range(loc_x$lon), ylim=range((loc_x$lat)))
  
  plot(p + facet_grid(rows = vars(Size_class_from), cols =  vars(Size_class_to),switch="y"))
  
}

# Format for pop. dyn model
temp <- Growth_temp2 %>% mutate(value=prob) %>% dplyr::select( lat,lon,Size_class_from,Size_class_to,value)
Growth_OM <- acast(temp, lat~lon~Size_class_from~Size_class_to)
Growth_OM <- ifelse(is.na(Growth_OM),0,Growth_OM)

return(Growth_OM)

}


