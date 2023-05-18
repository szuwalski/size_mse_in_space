## Demographic settings
#----------------------

##------------------------------------ Fishery selectivity -------------------------------------------------

# Parameterize with GMACS
selec_curve = repfile$selectivity$Start_Y %>% 
  filter(fleet == 1) %>% 
  dplyr::select_at(vars(starts_with("SizeC_")))

fish_sel = selec_curve

fish_sel_50_f<-NA
fish_sel_95_f<-NA
fish_sel_50_m<-50
fish_sel_95_m<-101

fish_sel<-rbind(1/(1+exp(-log(19)*(sizes-fish_sel_50_f)/(fish_sel_95_f-fish_sel_50_f))),
                1/(1+exp(-log(19)*(sizes-fish_sel_50_m)/(fish_sel_95_m-fish_sel_50_m))))
fish_sel[is.na(fish_sel)]<-0


##------------------------------- Calculate costs to a patch from harbour --------------------------------

# Dutch harbor
port_lat<-  54
port_lon<- -166.54

cost<-raster(nrow=length(lat), ncol=length(lon), 
             xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs="+proj=utm")
cost[]<-1
for(x in 1:nrow(cost))
  for(y in 1:ncol(cost))
    if(land_mask[x,y]==0) cost[x,y]<-100000

trCost <- transition(1/cost, min, directions=8)
trCostC <- geoCorrection(trCost, type="c")
trCostR <- geoCorrection(trCost, type="r")

## Create three points (representing three points in time series)
pnts <- cbind(x=c(port_lon, -155,-160), y=c(port_lat, 52,58))
costDistance(trCostC,pnts)

# ## Display results for one set of points
# plot(cost)
# plot(SpatialPoints(pnts), add=TRUE, pch=20, col="red")
# plot(shortestPath(trCostC, pnts[1,], pnts[2,], output="SpatialLines"), add=TRUE)
# plot(shortestPath(trCostC, pnts[1,], pnts[3,], output="SpatialLines"), add=TRUE)
# plot(shortestPath(trCostC, pnts[2,], pnts[3,], output="SpatialLines"), add=TRUE)

#==find the distance from harbor to all points
distance_map<-matrix(ncol=length(lon),nrow=length(lat))
colnames(distance_map)<-lon
rownames(distance_map)<-(lat)
for(x in 1:length(lon))
  for(y in 1:length(lat))
  {
    if(land_mask[y,x]!=0)
    {
      pts <- cbind(x=c(port_lon, as.numeric(colnames(distance_map)[x])), y=c(port_lat,  as.numeric(rownames(distance_map)[y])))
      distance_map[y,x]<-costDistance(trCostC, pts[1,],pts[2,])
      if(distance_map[y,x]>10000)distance_map[y,x]<-NA
    }
  }
# filled.contour(x=lon,y=rev(lat),g(distance_map*land_mask),plot.axes=c(map(add=TRUE,fill=T,col='grey'),
#                                                                       points(y=port_lat,x=port_lon,pch=16,col='red')))
# write.csv(distance_map,'dist.csv')


##------------------------------- Spatialized cost of fishing --------------------------------

cost_fish <- 0
cost_travel <- 0
cost_patch <- cost_travel * distance_map + cost_fish
price <- 1.5

fishers<-30
quota<-rep(1e7/fishers,fishers)

fishing_process="stochastic"
# fishing_process="max_benefit_min_dist"

