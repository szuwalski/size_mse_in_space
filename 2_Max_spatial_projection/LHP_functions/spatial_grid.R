spatial_grid <- function (lon,lat){
  

#-- Create a SpatialPoints object
set.seed(0)
data(wrld_simpl)
crs_LL = CRS(proj4string(wrld_simpl))
point_expand <- expand.grid(lon=lon, lat=lat) 

pts <- SpatialPoints(point_expand, proj4string=crs_LL)
proj4string(wrld_simpl)<-CRS(proj4string(pts))

#-- Find which points fall over land
land <- !is.na(over(pts, wrld_simpl)$FIPS)
point_expand$lat <- ifelse(land==TRUE,NA,point_expand$lat)
point_expand$lon <- ifelse(land==TRUE,NA,point_expand$lon)

#-- Highlight critters where there is land (this will be used for movement as well)
land_mask<-matrix(1,ncol=length(lon),nrow=length(lat),byrow=T)
for(x in 1:length(lat))
  for(y in 1:length(lon))
  {
    if(land[intersect(which(!is.na(match(pts$lon,(lon)[y]))) ,which(!is.na(match( pts$lat,(lat)[x]))))])
      land_mask[x,y]<-0
  }      

#-- Plot grid
world_sf <- st_as_sf(wrld_simpl,crs=CRS(proj4string(wrld_simpl)))
p_grid <- ggplot()+ geom_sf(data = world_sf,fill="black",color=NA)+ coord_sf(xlim=range(lon), ylim=range((lat)))+
  geom_point(data=point_expand,map=aes(x=lon,y=lat),col="orange") + theme_bw()

ggsave((paste0(DateFile, 'Spatial_grid.png')),plot=p_grid,
       width = 27,
      height = 18,
     units = "cm")

#-- Define spatial 
loc_x <- na.omit(point_expand)
n_s <- dim(loc_x)[1]
cells <- c(1:n_s)
loc_x <- cbind(cells,loc_x)

return(list("n_s"=n_s,"loc_x"=loc_x,"land_mask"=land_mask, "land"=land))
}