## Spatial extent of the study area
#----------------------------------
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
plot(wrld_simpl,xlim = c(min(lon), max(lon)), ylim = c(min(lat),max(lat)))
points(pts, col=1+land, pch=16)
