## Abundance at size
#-------------------

sexN = 2


## Load abundance at size from stock assessment
init_year = 37 # Time series is 1982-2022 then time step 37 is 2019
Ab_males = unlist(Snow_Out$Abundance$N_males[init_year,])
Ab_males_mat = unlist(Snow_Out$Abundance$N_males_mature[init_year,])
Ab_females = unlist(Snow_Out$Abundance$N_females[init_year,])
Ab_females_mat = unlist(Snow_Out$Abundance$N_females_mature[init_year,])

## Load and shape survey data from 2019 for ventilating abundance in space
ebs_2019 <- read_csv("4_full_MSE/data/ebs_2019.csv")
nbs_2019 <- read_csv("4_full_MSE/data/nbs_2019.csv")

ebs_2019_2 = ebs_2019 %>% 
  mutate(area = "ebs") %>% 
  mutate(stage = ifelse(WIDTH < 45,"Juv",'Mat')) %>% 
  dplyr::group_by(HAUL,MID_LATITUDE,MID_LONGITUDE,stage,area) %>% 
  tally()

nbs_2019_2 = nbs_2019 %>% 
  mutate(area = "nbs") %>% 
  mutate(stage = ifelse(WIDTH < 45,"Juv",'Mat')) %>% 
  dplyr::group_by(HAUL,MID_LATITUDE,MID_LONGITUDE,stage,area) %>% 
  tally()

count_df = rbind(ebs_2019_2,nbs_2019_2)

point_sf = st_as_sf(point_expand,coords = c("Var1","Var2"))
raster_dom = st_rasterize(point_sf %>% dplyr::select(key, geometry))
grid_sf = st_as_sf(raster_dom)

count_sf = st_as_sf(count_df,coords = c("MID_LONGITUDE","MID_LATITUDE"))
count_grid_sf = st_intersection(grid_sf,count_sf) %>%
  as.data.frame %>% 
  group_by(key,stage) %>%
  dplyr::summarise(n = mean(n)) %>% 
  full_join(grid_sf) %>%
  filter(!is.na(n)) %>% 
  filter(!is.na(stage)) %>% 
  st_as_sf

# ## Plot survey data
# lat_range = range(count_df$MID_LATITUDE)
# lon_range = range(count_df$MID_LONGITUDE)
# 
# Juv_plot = ggplot(count_df[which(count_df$stage == "Juv"),])+
#   geom_point(aes(x=MID_LONGITUDE,y=MID_LATITUDE,col=log(n)),size = 2)+
#   scale_color_distiller(palette = "Spectral")+
#   ggtitle("Juveniles")+
#   xlim(lon_range)+ylim(lat_range) # + facet_wrap(.~area)
# 
# Mat_plot = ggplot(count_df[which(count_df$stage == "Mat"),])+
#   geom_point(aes(x=MID_LONGITUDE,y=MID_LATITUDE,col=log(n)),size = 2)+
#   scale_color_distiller(palette = "Spectral")+
#   ggtitle("Mature")+
#   xlim(lon_range)+ylim(lat_range) # + facet_wrap(.~area)
# 
# count_grid_plot = ggplot()+
#   geom_sf(data=count_grid_sf,aes(fill=n))+
#   scale_fill_distiller(palette="Spectral")+
#   facet_wrap(.~stage)

## Compute spatial matrices of relative biomass distribution for juveniles and matures
init_adult<-array(0,dim=c(length(lat),length(lon)))
init_juv<-array(0,dim=c(length(lat),length(lon)))

for(x in 1:length(lat))
  for(y in 1:length(lon))
  {
    
    key = point_expand$key[which(point_expand$Var1 == lon[y] & point_expand$Var2 == lat[x])]
    count_grid_mat = count_grid_sf$n[count_grid_sf$key == key & count_grid_sf$stage == "Mat"]
    count_grid_immat = count_grid_sf$n[count_grid_sf$key == key & count_grid_sf$stage == "Juv"]
    if(length(count_grid_mat) > 0) init_adult[x,y] = count_grid_sf$n[count_grid_sf$key == key & count_grid_sf$stage == "Mat"]
    if(length(count_grid_immat) > 0) init_juv[x,y] = count_grid_sf$n[count_grid_sf$key == key & count_grid_sf$stage == "Juv"]
    
  }

init_adult = init_adult / sum(init_adult)
init_juv = init_juv / sum(init_juv)

# ## Plot matrices
# cowplot::plot_grid(Juv_plot,Mat_plot)
# par(mfrow=c(1,2))
# plot(init_juv * land_mask_na, main = "", asp = 1)
# plot(init_adult * land_mask_na, main = "", asp = 1)

## Ventilate total abundance of stock assessment in space
# Make matrices
imm_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
mat_N_at_Len<-array(dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))

# Fill abundance matrices
imm_N_at_Len[,,1,1,1]<-init_juv * (Ab_females[1] - Ab_females_mat[1])
imm_N_at_Len[,,2,1,1]<-init_juv * (Ab_males[1] - Ab_males_mat[1])
mat_N_at_Len[,,1,1,1]<-init_juv * Ab_females_mat[1]
mat_N_at_Len[,,2,1,1]<-init_juv * Ab_males_mat[1]

for(x in 2:length(sizes))
{
  
  if(sizes[x] < ad_size){
    imm_N_at_Len[,,1,x,1]<-init_juv * (Ab_females[x] - Ab_females_mat[x])
    imm_N_at_Len[,,2,x,1]<-init_juv * (Ab_males[x] - Ab_males_mat[x])
    mat_N_at_Len[,,1,x,1]<-init_juv * Ab_females_mat[x]
    mat_N_at_Len[,,2,x,1]<-init_juv * Ab_males_mat[x]
  }else if(sizes[x] > 45){
    imm_N_at_Len[,,1,x,1]<-init_adult * (Ab_females[x] - Ab_females_mat[x])
    imm_N_at_Len[,,2,x,1]<-init_adult * (Ab_males[x] - Ab_males_mat[x])
    mat_N_at_Len[,,1,x,1]<-init_adult * Ab_females_mat[x]
    mat_N_at_Len[,,2,x,1]<-init_adult * Ab_males_mat[x]
  }
  
}

imm_N_at_Len[imm_N_at_Len=="NaN"]<-0
mat_N_at_Len[mat_N_at_Len=="NaN"]<-0

imm_N_at_Len[imm_N_at_Len<0]<-0
mat_N_at_Len[mat_N_at_Len<0]<-0

#==ensures no critters on land
for(x in 1:2)
  for(y in 1:length(sizes))
  {
    imm_N_at_Len[,,x,y,1]<- imm_N_at_Len[,,x,y,1]*land_mask
    mat_N_at_Len[,,x,y,1]<- imm_N_at_Len[,,x,y,1]*land_mask
  }

#==SHOULD THERE ALSO BE A 'DEPTH MASK'???
#==MAYBE A 'FISHERY MASK'??  Depth should limit the fishery.
# filled.contour(x=lon,y=rev(lat),g(imm_N_at_Len[,,1,1,1]))

#==create a file with bottom temperature for all time periods
avg_bot_tmp<-rep(c(1,2,3,3,3,2,1,0,0,0,0,1),year_n)
for(x in 1:length(proj_period))
{
  temp<-land_mask*rnorm(length(land_mask),avg_bot_tmp[x],.5)
  #write.csv(temp,paste("temp_data/bot_temp_",x,".csv",sep=""),row.names=FALSE)
}

if(size_class_settings == "rough") area_mask = "EBS_only"

