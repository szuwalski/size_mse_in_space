## Shape commercial data
#-----------------------
if(print_messages) print("Make commercial data")

com_it = 0
for(x in 1:length(lat))
{
  for(y in 1:length(lon))
  {
    for(i in 1:length(sizes))
    {
      
      include_loc = T
      if(sp_domain == "EBS" & EBS_mask_matrix[x,y] == 1) include_loc = F
      
      if(land_mask[x,y]!=0 & include_loc)
      {
        
        #------------------------------------------------------------------------------------------------------------------------
        # Need to implement observation pdf
        # print(paste0("x:",x,"|y:",y,"|i:",i))
        if(com_it==0){
          
          catch_mat = matrix(data = c(Years_climsc[t], # Year
                                      i, # size_bin
                                      sum(total_spatial_catch[x,y,2,i,which(Years_climsc == Years_climsc[t])]), # Catches_N: only male
                                      lon[y], # Lon
                                      lat[x]), # Lat
                             nrow = 1)
          com_it = com_it + 1
          
        }else{
          
          line_vec = matrix(data = c(Years_climsc[t], # Year
                                     i, # size_bin
                                     sum(total_spatial_catch[x,y,2,i,which(Years_climsc == Years_climsc[t])]), # Catches_N: only male
                                     lon[y], # Lon
                                     lat[x]), # Lat
                            nrow = 1)
          
          catch_mat = rbind(catch_mat,line_vec)
          
        }
        
      }
    }
  }
}

catch_df = as.data.frame(catch_mat)
colnames(catch_df) = colnames(catch_N)
catch_df$Year = as.character(catch_df$Year)
catch_df$size_bin = as.character(catch_df$size_bin)
catch_N = rbind(catch_N,catch_df)

catch_N %>% 
  filter(Year > 2018) %>% 
  group_by(Year) %>% 
  dplyr::summarise(Catches_N = sum(Catches_N))
