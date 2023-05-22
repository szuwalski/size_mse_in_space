## Sample scientific data
#------------------------
if(print_messages) print("Simulate scientific data")
## Implemented through matrix because much faster
sci_it = 0
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
        
        # Need to implement observation pdf
        # print(paste0("x:",x,"|y:",y,"|i:",i))
        if(sci_it == 0){
          
          sci_mat = matrix(data = c(i, # size_class
                                    Years_climsc[t], # year
                                    imm_N_at_Len[x,y,2,i,t] + mat_N_at_Len[x,y,2,i,t], # Catch_N: only male
                                    cell_area[x,y], # AreaSwept_km2
                                    0, # Vessel
                                    lat[x], # Lat
                                    lon[y]), # Lon
                           nrow = 1)
          
          sci_it = sci_it + 1
          
        }else{
          
          line_vec = matrix(data = c(i, # size_class
                                     Years_climsc[t], # year
                                     imm_N_at_Len[x,y,2,i,t] + mat_N_at_Len[x,y,2,i,t], # Catch_N: only male
                                     cell_area[x,y], # AreaSwept_km2
                                     0, # Vessel
                                     lat[x], # Lat
                                     lon[y]), # Lon
                            nrow = 1)
          sci_mat = rbind(sci_mat,line_vec)
          
        }
      }
    }
  }
}

sci_df = as.data.frame(sci_mat)
colnames(sci_df) = colnames(Data_Geostat)
sci_df$size_class = as.factor(sci_df$size_class)
Data_Geostat = rbind(Data_Geostat,sci_df)