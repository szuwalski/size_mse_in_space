
# -------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------
  library(dplyr)
  library(sp) 
  library(sf)
  library(rgdal)
  library(sf)
  library(sp)
  #library(INLA)
  library(RandomFields)
  library(maptools)
  library(ggplot2)
  library(gstat)
  library(BBmisc)
  library(viridis)
  library(patchwork)
  library(raster)
  library(spData) ## For `world`, an sf MULTIPOLYGON object
  require(virtualspecies)
  library(scales)
  library(BBmisc)
  library(reshape2)


# Output file -------------------------------------------------------------
# -------------------------------------------------------------------------
DateFile = paste0(getwd(), "/Outputs/")
dir.create(DateFile)  
  
# Define grid to predict the Field-----------------------------------------
# -------------------------------------------------------------------------
lat		   <-seq(70,51.5,length.out=40)
lon		   <-seq(-179,-155,length.out=40)
point_expand <- expand.grid(Lon=lon, Lat=lat) 

# Remove lat/lon inland  
data(wrld_simpl)

## Create a SpatialPoints object
set.seed(0)
crs_LL = CRS(proj4string(wrld_simpl))
pts <- SpatialPoints(point_expand, proj4string=crs_LL)
proj4string(wrld_simpl)<-CRS(proj4string(pts))

## Find which points fall over land
land <- !is.na(over(pts, wrld_simpl)$FIPS)
point_expand$Lat <- ifelse(land==TRUE,NA,point_expand$Lat)
point_expand$Lon <- ifelse(land==TRUE,NA,point_expand$Lon)

# Plot grid
world_sf <- st_as_sf(wrld_simpl,crs=CRS(proj4string(wrld_simpl)))
p_grid <- ggplot()+ geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(lon), ylim=range((lat)))+
  geom_point(data=point_expand,map=aes(x=Lon,y=Lat),col="orange") + theme_bw()

ggsave((paste0(DateFile, 'Spatial_grid.png')),plot=p_grid,
       width = 27,
       height = 18,
       units = "cm")

# Define spatial 
loc_x <- na.omit(point_expand)


# -------------------------------------------------------------------------
# Year/ climate scenarios/ preferntial habitat function
# -------------------------------------------------------------------------

# choose years
years <- c(2030,2050,2070)

# choose climate scenarios
clim_sc <- c("rcp45", "rcp85")


# -------------------------------------------------------------------------
# OMEGA
# -------------------------------------------------------------------------

# 1 - Simulate Field conditioned on climate scenario----------------------
# -------------------------------------------------------------------------
# Load Temperature scenarios  ---------------------------------------------
load("climate_scenarios/mn_var_all.RData")
En.cond <- mn_var_all %>% filter(simulation %in% clim_sc , year %in% years ) %>% dplyr::select(year,simulation, latitude,longitude, val)
# transform long to get long between -180 and 180
En.cond <- En.cond %>% mutate(Lon = ifelse(longitude <= 180, longitude, longitude -360),Lat=latitude)
dat_cond <-  data.frame(coords=En.cond[,c("Lon","Lat")],Year=En.cond$year, Temp=En.cond$val, Climate_sc = En.cond$simulation)
colnames(dat_cond) <- c( "Lon","Lat","Year","Temp","Climate_sc")
# transform in sf object
dat_cond_sf <- st_as_sf(dat_cond, coords=c("Lon","Lat"), crs=crs_LL)


# Simulate field conditioned by env. variations
En.cov_plot <- NULL
for (y in 1: length(years)) {
  for (s in 1:length(clim_sc)){
    temp_sc_year_sf <- NULL
    En.cov_plot_temp <- NULL
    temp_sc_year_sf <- (dat_cond_sf) %>% dplyr::filter(Year == years[y],Climate_sc == clim_sc[s])
    
    dat_cond_temp <- data.frame(coords=(st_coordinates(temp_sc_year_sf)),Temp=temp_sc_year_sf$Temp)
    
    #  RMexp() is a stationary isotropic covariance model whose corresponding covariance function
    #only depends on the distance r ??? 0 between two points and is given by C(r) = exp(-r)
    model_omega <-  RMexp()
    En.cov_temp <- RFsimulate(model = model_omega, x=loc_x[,1], y=loc_x[,2],data= dat_cond_temp)
    En.cov_plot_temp <- data.frame(Year = years[y],Climate_sc = clim_sc[s], Temp = En.cov_temp$Temp,Lon=loc_x[,1],Lat=loc_x[,2])
    
    En.cov_plot <- bind_rows(En.cov_plot, En.cov_plot_temp)
  }
}



p1 <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=En.cov_plot ,aes(Lon, Lat,color=Temp), size = 1.5) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")+
  facet_grid( Climate_sc~Year )+ ggtitle("Spatial field conditonned on Bottom Temperature ")

p2 <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=dat_cond ,aes(Lon, Lat,color=Temp), size = 1.5) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")+
  facet_grid( Climate_sc~Year ) + ggtitle("Bottom Temperature projection (level3)")


pomega_1 <- p2 +p1

ggsave((paste0(DateFile,"Omega_step1.png")),plot=pomega_1,
       width = 27,
       height = 18,
       units = "cm")



# 2-Apply Preferential Habitat function to Field---------------------------
# -------------------------------------------------------------------------
# choose Preference habitat
pref_hab <- c("logi","quadr","line")
# define parameters of the preferential habitat funcion 
pars <- matrix(0,2,length(pref_hab))
# logistic
pars[1,1] <- -1
pars[2,1] <-  1
# quad
pars[1,2] <- 4
pars[2,2] <- 2
# linea
pars[1,3] <- 0.5
pars[2,3] <-  0

X_s <- NULL
for( p in 1:length(pref_hab)){
  for (y in 1: length(years)) {
    for (s in 1:length(clim_sc)){
      
      X_s_temp <- NULL
      rast <- raster(ncol=150,nrow=150)
      loc_x_tmp <- loc_x
      colnames(loc_x_tmp) <- c("x" ,"y")
      extent(rast) <- extent(as.data.frame(loc_x_tmp))
      
      P2 = SpatialPoints(loc_x_tmp)
      data_temp <-   En.cov_plot %>% dplyr::filter(Year == years[y],Climate_sc == clim_sc[s])
      P2$data =     data_temp$Temp
      rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
      
      envir_stack <- stack(rast.temp)
      names(envir_stack) <- c('temp')
      
      # Apply the habitat preference function
      p_hab = pref_hab[p]
      if (p_hab=="logi") fun_hab <- c(fun='logisticFun', alpha=pars[1,1], beta=pars[2,1])
      if (p_hab=="quadr") fun_hab <-  c(fun="dnorm",mean=pars[1,2],sd=pars[2,2])
      if (p_hab=="line") fun_hab <-  c(fun="linearFun",a=pars[1,3],b=pars[2,3])
      # Plot the Preferential habitat function

       # Use Virtual Species to assign response curve
      parameters <- virtualspecies::formatFunctions(temp = fun_hab)
      envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)
      
      png(paste(DateFile,'PH_fun',p_hab,'.png',sep=''), height = 7, width = 8, units = 'in', res=600)
      virtualspecies::plotResponse(envirosuitability) #plot response curves
      dev.off()
      
      X_s_temp <- data.frame(Lat=loc_x[,2],Lon=loc_x[,1], Temp_pref = raster::extract(envirosuitability$suitab.raster,loc_x),Year = rep(years[y], length(loc_x[,2])),Climate_sc = rep(clim_sc[s],length(loc_x[,2])), PH_fun = rep(p_hab,length(loc_x[,2])) )   
      X_s <- bind_rows(X_s, X_s_temp)
    }
  }
}


# 3- Rescale X_s to get Omega varying within a range ----------------------
# -------------------------------------------------------------------------

# Average c(mean, sd)
theta_init <- c(3,0.1)
theta_last <- c(15,0.1)

# Weight of omega :  which percentage of the average value theta_init does the spatial variation omega represent
x_om_last <- 0.1
x_om_init <- 0.1

# Rescale
Omega =NULL
for( p in 1:length(pref_hab)){
  for (y in 1: length(years)) {
    for (s in 1:length(clim_sc)){
      Omega_temp=NULL
      X_s_temp = NULL
      X_s_temp <- X_s %>% filter(Year == years[y],Climate_sc == clim_sc[s],PH_fun== pref_hab[p])
      Omega_last <-  rescale(X_s_temp$Temp_pref, to = c(-x_om_last*theta_last[1],    x_om_last*theta_last[1]))
      Omega_init <-  rescale(X_s_temp$Temp_pref, to = c(-x_om_init*theta_init[1],    x_om_init*theta_init[1]))
      Omega_temp <- data.frame(init=Omega_init, last =Omega_last,
                               Year = rep(years[y],length(Omega_init)),
                               Lat=X_s_temp$Lat,Lon =X_s_temp$Lon,
                               Climate_sc = rep(clim_sc[s],length(Omega_init)),
                               PH_fun=rep(pref_hab[p],length(Omega_init)))
      Omega <- bind_rows(Omega,Omega_temp)
    }}}


p_omega_rcp45 <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=Omega %>% filter(Climate_sc=="rcp45") ,aes(Lon, Lat,color=init), size = 1.5) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")+
  facet_grid( PH_fun~Year ) + ggtitle("Omega (RCP45) \n parameter= Theta_init (growth increment in the first class)")

p_omega_rcp85 <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=Omega %>% filter(Climate_sc=="rcp85") ,aes(Lon, Lat,color=init), size = 1.5) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")+
  facet_grid( PH_fun~Year ) + ggtitle("Omega(RCP85) \n parameter Theta_init (growth increment in the first class)")


p_omega <- p_omega_rcp45 +p_omega_rcp85

ggsave((paste0(DateFile,"Omega.png")),plot=p_omega,
       width = 30,
       height =30 ,
       units = "cm")


# -------------------------------------------------------------------------
#  EPSILON ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Define variance and scale for RF
scale_g <- 1.5 # scale of the RF

Epsilon <-NULL
for(i in 1:length(years)){
  # define model
  # RMgauss a stationary isotropic covariance model. The corresponding covariance function only
  # depends on the distance r ??? 0 between two points and is given by
  # C(r)=e^{-r^2}
  model_epsilon_m_last       <- RMgauss(var=theta_last[2]^2, scale=scale_g)
  model_epsilon_m_init     <- RMgauss(var=theta_init[2]^2, scale=scale_g)
  
  # Simulate RF
  theta_last_epsilon = RFsimulate(model = model_epsilon_m_last, x=loc_x[,1], y=loc_x[,2])@data[,1]
  theta_init_epsilon = RFsimulate(model = model_epsilon_m_init, x=loc_x[,1], y=loc_x[,2])@data[,1]
  hist(theta_last_epsilon)

  # Rescale Epsilon 
  x_eps_last=0.05 
  x_eps_init=0.05 
  

  Epsilon_theta_last <- rescale(theta_last_epsilon, to = c(-x_eps_last*theta_last[1],    x_eps_last*theta_last[1]))
  Epsilon_theta_init <- rescale(theta_init_epsilon, to = c(-x_eps_init*theta_init[1],    x_eps_init*theta_init[1]))
  
  Epsilon_temp <- data.frame(init=Epsilon_theta_init, last =Epsilon_theta_last,
                             Lat=loc_x[,2],Lon =loc_x[,1],
                             Year = rep(years[i],length(Epsilon_theta_init)))
  
  
  Epsilon <- bind_rows(Epsilon,Epsilon_temp)
  
}

p_epsilon <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=Epsilon  ,aes(Lon, Lat,color=init), size = 4.5) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")+
  facet_grid( ~Year ) + ggtitle("Epsilon \n parameter Theta_init (growth increment in the first class)")


ggsave((paste0(DateFile,"Epsilon.png")),plot=p_epsilon,
       width = 30,
       height =30 ,
       units = "cm")


# -------------------------------------------------------------------------
# Simulate growth ---------------------------------------------------------
# -------------------------------------------------------------------------
source("growth/script_growth_inc.R")

# Arguments for growth function
# size range
binsize  <-5
sizes		 <-seq(27.5,132.5,binsize)
n_grpar = 3 # number of parameter to define growth 
# nber of location
n_s <- dim(loc_x)[1]
#nber of size bins
n_p = as.numeric(length(sizes))
# nber of years
n_t <- length(years)

# scale parameter of the gamma distribution used to generate the growth transition matrix
beta_scale= c(0.3)

Growth = NULL

for( p in 1:length(pref_hab)){
  for (s in 1:length(clim_sc)){
    growth_matrix   = array (data = 0, dim = c(n_t,n_s,n_p,n_p))
    growth_par = array (data = NA, dim = c(n_s,n_grpar,n_t)) 
    
    for (y in 1: length(years)) {
      Growth_temp <- NULL
      
      Omega_temp <- Omega %>% dplyr::filter(Year == years[y],Climate_sc == clim_sc[s],PH_fun== pref_hab[p])
      Epsilon_temp <- Epsilon %>% dplyr::filter(Year == years[y])
      
      growth_par[,1,y] =  theta_init[1] + Epsilon_temp$init +  Omega_temp$init
      growth_par[,2,y] = theta_last[1] + Epsilon_temp$last + Omega_temp$last
      growth_par[,3,y] = beta_scale[1]
      ### define a function with dplyr group_by
      
      for (k in 1:n_s) {
        growth_matrix[y,k,,] = growth_trans(sizes,
                                                 growth_par[k, , y ][2],
                                                 growth_par[k, , y ][1],
                                                 sizes[1],
                                                 sizes[n_p],
                                                 binsize,
                                                 beta_scale)
      }
    }
    Growth_temp <- melt(growth_matrix)
    colnames(Growth_temp) <- c("Year","cell","Size_class_from","Size_class_to", "prob") 
    Growth_temp <- data.frame(Growth_temp, Climate_sc = rep(clim_sc[s],dim(Growth_temp)[1]),PH_fun=rep(pref_hab[p],dim(Growth_temp)[1]))
    
    Growth <- bind_rows( Growth,Growth_temp)
    
  }}



# add locations to growth
Cell <- data.frame(cell = unique(Growth$cell), Lat=loc_x[,2],Lon=loc_x[,1])
Growth_df <- left_join(Growth,Cell)


# Plot some examples
# -------------------------------------------------------------------------


# Prob 1 -> 2 Plot 2 scenarios 3 years
Temp  <- Growth_df %>% filter(Size_class_from==1,Size_class_to==2, Climate_sc =="rcp85" ) 
Temp2 <- data.frame(Year = c(1:3), Years= years) # transform years from c(1,2,3) to c(2030,2050,2070)
Temp3 <- left_join(Temp,Temp2)

p_growth_1_2 <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=Temp3,aes(Lon, Lat, color = prob), size = 2) +
  scale_color_viridis()  +
  theme_bw()+ ylab("") +xlab("")+
  facet_grid( Years ~ PH_fun ) + ggtitle(" Probability of growing from size class 1 to 2")


ggsave((paste0(DateFile,"Growth_PH_funvsYear.png")),plot=p_growth_1_2,
       width = 27,
       height = 18,
       units = "cm")

