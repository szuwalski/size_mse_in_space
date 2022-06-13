

growth_map <- function(G_spatial , lon, lat,clim_sc, pref_hab, years,pars,x_om_init,x_om_last,x_eps_last,x_eps_init,beta_scale,binsize,sizes, m_last, m_init,scale_g){

  library(dplyr)
  library(sp) 
  library(sf)
  library(rgdal)
  library(sf)
  library(sp)
  library(INLA)
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
 # RFoptions(seed=0)   
   
 # lon #longitutde
  #lat #latitude
  #clim_sc=c("rcp45") # climate scenario
  #pref_hab= c("quadr") # peferential habitat function (logi, quad, lin)
  #pars = c(4,2)# parameters to define the preferential habitat function
  #years = 2050
  #x_om_init = 0.1 # if additive effect, % of the mean of the initial growth increment represented by omega (field driven by climate sc.)
  #x_om_last=0.1 # if additive effect, % of the mean of the last growth increment represented by omega (field driven by climate sc.)
  #x_eps_last=0.05 # if additive effect, % of the mean of the last growth increment represented by epsilon (random field)
  #x_eps_init=0.05 # if additive effect, % of the mean of the intial growth increment represented by epsilon (random field.)
  #beta_scale= c(0.1)# scale parameter of the gamma distribution used to generate the growth transition matric 
  #binsize # bin size
  #sizes # size classes
  
  #m_last = c(30,0.1) # last size increment : meam, sd
  #m_init = c(8,0.1) # initial size increment  : meam, sd
  
  #scale_g <- 0.1 # scale of the RF



# Output file -------------------------------------------------------------
# -------------------------------------------------------------------------
DateFile = paste0( "Outputs/")
dir.create(DateFile)  
  
# Define grid to predict the Field-----------------------------------------
# -------------------------------------------------------------------------
point_expand <- expand.grid(Lon=lon, Lat=lat) 

# Remove lat/lon inland  
point_expand$Lat <- ifelse(land==TRUE,NA,point_expand$Lat)
point_expand$Lon <- ifelse(land==TRUE,NA,point_expand$Lon)

# Plot grid
world_sf <- st_as_sf(wrld_simpl,crs=CRS(proj4string(wrld_simpl)))
p_grid <- ggplot()+ geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(lon), ylim=range((lat)))+
  geom_point(data=point_expand,map=aes(x=Lon,y=Lat),col="orange") + theme_bw()

ggsave((paste0(DateFile, '/Spatial_grid.png')),plot=p_grid,
       width = 27,
       height = 18,
       units = "cm")


# Define spatial 
loc_x <- na.omit(point_expand)
crs_LL = CRS(proj4string(wrld_simpl))

# Arguments for growth fuction
n_p = as.numeric(length(sizes))


# -------------------------------------------------------------------------
# OMEGA
# -------------------------------------------------------------------------

# 1 - Simulate Field conditionned on climate scenario----------------------
# -------------------------------------------------------------------------

# Load Temperature scenarios  ---------------------------------------------
load("climate_scenarios/mn_var_all.RData")
En.cond <- mn_var_all %>% filter(simulation %in% clim_sc , year %in% years ) %>% dplyr::select(year,simulation, latitude,longitude, val)
En.cond <- En.cond %>% mutate(Lon = ifelse(longitude <= 180, longitude, longitude -360),Lat=latitude)
dat_cond <-  data.frame(coords=En.cond[,c("Lon","Lat")],Year=En.cond$year, Temp=En.cond$val, Climate_sc = En.cond$simulation)
colnames(dat_cond) <- c( "Lon","Lat","Year","Temp","Climate_sc")
dat_cond_sf <- st_as_sf(dat_cond, coords=c("Lon","Lat"), crs=crs_LL)


# Extract year and scenarios ----------------------------------------------
temp_sc_year_sf <- NULL
En.cov_plot <- NULL
temp_sc_year_sf <- (dat_cond_sf) %>% dplyr::filter(Year == years,Climate_sc == clim_sc)
dat_cond_temp <- data.frame(coords=(st_coordinates(temp_sc_year_sf)),Temp=temp_sc_year_sf$Temp)

# Simulate Field------------------------------------------------------------
spde_cov <-  RMexp() #RMgauss()#var=0.1*0.1,scale=scale_g
En.cov_temp <- RFsimulate(model = spde_cov, x=loc_x[,1], y=loc_x[,2],data= dat_cond_temp)
En.cov_plot <- data.frame( Temp = En.cov_temp$Temp,Lon=loc_x[,1],Lat=loc_x[,2])

p_RFcond <- ggplot() + geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(loc_x$Lon), ylim=range((loc_x$Lat)))+
  geom_point(data=En.cov_plot,aes(Lon, Lat,color=Temp), size = 2) +
  scale_color_viridis()  +
  theme_bw() + ylab("") +xlab("")
p_RFcond

# 2-Apply Preferential Habitat function to Field---------------------------
# -------------------------------------------------------------------------
X_s <- NULL

X_s_temp <- NULL
rast <- raster(ncol=150,nrow=150)
loc_x_tmp <- loc_x
colnames(loc_x_tmp) <- c("x" ,"y")
extent(rast) <- extent(as.data.frame(loc_x_tmp))
  
P2 = SpatialPoints(loc_x_tmp)
data_temp <-   En.cov_plot 
P2$data =     data_temp$Temp
rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    
envir_stack <- stack(rast.temp)
names(envir_stack) <- c('temp')
  
# Apply the habitat preference function
if (pref_hab=="logi") fun_hab <- c(fun='logisticFun', alpha=pars[1], beta=pars[2])
if (pref_hab=="quadr") fun_hab <-  c(fun="dnorm",mean=pars[1],sd=pars[2])
if (pref_hab=="line") fun_hab <-  c(fun="linearFun",a=pars[1],b=pars[2])
    
parameters <- virtualspecies::formatFunctions(temp = fun_hab)
envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)

X_s <- data.frame(Lat=loc_x[,2],Lon=loc_x[,1], Temp_pref =   raster::extract(envirosuitability$suitab.raster,loc_x) )   

# Plot suitability function    
png(paste('PH_fun',pref_hab,'.png',sep=''), height = 7, width = 8, units = 'in', res=600)
virtualspecies::plotResponse(envirosuitability) #plot response curves
dev.off()


# 3- Rescale X_s to get Omega varying within a range ----------------------
# -------------------------------------------------------------------------

Omega_last <-   rescale(X_s$Temp_pref, to = c(-x_om_last*m_last[1],    x_om_last*m_last[1]))
Omega_init <-   rescale(X_s$Temp_pref, to = c(-x_om_init*m_init[1],    x_om_init*m_init[1]))
Omega <- data.frame(init=Omega_init, last =Omega_last)

# -------------------------------------------------------------------------
#  EPSILON ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Define variance and scale for RF

# define model
model_epsilon_m_last       <- RMgauss(var=m_last[2]^2, scale=scale_g)
model_epsilon_m_init     <- RMgauss(var=m_init[2]^2, scale=scale_g)

# Simulate RF
m_last_epsilon = RFsimulate(model = model_epsilon_m_last, x=loc_x[,1], y=loc_x[,2])@data[,1]
m_init_epsilon = RFsimulate(model = model_epsilon_m_init, x=loc_x[,1], y=loc_x[,2])@data[,1]
  
# Rescale Epsilon 

Epsilon_m_last <- rescale(m_last_epsilon, to = c(-x_eps_last*m_last[1],    x_eps_last*m_last[1]))
Epsilon_m_init <- rescale(m_init_epsilon, to = c(-x_eps_init*m_init[1],    x_eps_init*m_init[1]))

Epsilon <- data.frame(init=Epsilon_m_init, last =Epsilon_m_last)


# Simulate growth ---------------------------------------------------------
# -------------------------------------------------------------------------
source("function_LHP/script_growth_inc.R")

n_grpar = 3 # number of parameter to define growth 
n_s <- dim(loc_x)[1]
growth_par = array (data = NA, dim = c(n_s,n_grpar))
  


# Map of growth matrix
growth_matrix   = array (data = 0, dim = c(n_s,n_p,n_p))
growth_par = array (data = NA, dim = c(n_s,n_grpar))
      
growth_par[,1] =  m_init[1] + Epsilon$init +  Omega$init
growth_par[,2] = m_last[1] + Epsilon$last + Omega$last
growth_par[,3] = beta_scale


if(!G_spatial){
  for (k in 1:n_s) {
    growth_matrix[k,,] = growth_trans(sizes,m_last[1],
                                      m_init[1],
                                      sizes[1],
                                      sizes[n_p],
                                      binsize,
                                      beta_scale)
  }
  }else{
  for (k in 1:n_s) {
    growth_matrix[k,,] = growth_trans(sizes,growth_par[k,2],
                                           growth_par[k,1],
                                           sizes[1],
                                           sizes[n_p],
                                           binsize,
                                           growth_par[k,3])
    }
}

# match lat/lon with cell number
Growth_temp <- melt(growth_matrix)
colnames(Growth_temp) <- c("cell","Size_class_from","Size_class_to", "prob") 
Growth_temp <- data.frame(Growth_temp)

Growth_temp$Size_class_from <- as.factor(Growth_temp$Size_class_from)
Growth_temp$Size_class_to <- as.factor(Growth_temp$Size_class_to)

# Plot
p <- ggplot(Growth_temp %>% filter(cell %in% c(seq(1,1000,1000))), aes(Size_class_from, Size_class_to, fill= prob)) + 
  geom_tile()+scale_fill_viridis() +  theme_bw()

p + facet_wrap(~cell)

Cell <- data.frame(cell = unique(Growth_temp$cell), Lat=loc_x[,2],Lon=loc_x[,1])
Growth <- left_join(Growth_temp,Cell)

# Transform to appropriate format for OM input
lon2 <- data.frame(Lon=lon,lon_nb = c(1:length(lon)))
lat2 <- data.frame(Lat=lat,lat_nb = c(1:length(lat)))
Growth_temp <- left_join(Growth, lon2)
Growth_temp2 <- left_join(Growth_temp, lat2)
Growth_temp3 <- melt(Growth_temp2[,c("Size_class_from","Size_class_to","lon_nb","lat_nb","prob")])

temp <- Growth_temp2 %>% mutate(value=prob) %>% dplyr::select( lat_nb,lon_nb,Size_class_from,Size_class_to,value)
Growth_OM <- acast(temp, lat_nb~lon_nb~Size_class_from~Size_class_to)
Growth_OM <- ifelse(is.na(Growth_OM),0,Growth_OM)

return(Growth_OM)

}




