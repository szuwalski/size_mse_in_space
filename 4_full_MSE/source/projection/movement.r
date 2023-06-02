## Movement
#----------
# From Jim's Equation
# N.b. this does dot model ontongnic movement
# only preferential taxis habitat movement

for(sex in 1:2)
  for(x in 1:length(sizes))
  {
    
    if(sizes[x] < ad_size){
      
      mov_imm_temp = temp_imm_N[,,sex,x]
      mov_mat_temp = temp_mat_N[,,sex,x]
      
      mov_imm_v = mfraction_gg_juv %*% as.vector(t(mov_imm_temp))
      mov_mat_v = mfraction_gg_juv %*% as.vector(t(mov_imm_temp))
      
      mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
      mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
      
    }else{
      
      mov_imm_temp = temp_imm_N[,,sex,x]
      mov_mat_temp = temp_mat_N[,,sex,x]
      
      mov_imm_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp))
      mov_mat_v = mfraction_gg_ad %*% as.vector(t(mov_mat_temp))
      
      mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
      mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
      
    }
    
    temp_imm_N[,,sex,x] = mov_imm_temp_2
    temp_mat_N[,,sex,x] = mov_mat_temp_2
    
    
    # sum(mov_mat_temp) / sum(mov_mat_temp_2)
    
    # ## Check that movement happens
    # x11()
    # par(mfrow = c(4,4))
    # 
    # mov_imm_temp = temp_imm_N[,,sex,x]
    # mov_mat_temp = temp_mat_N[,,sex,x]
    # 
    # plot(mov_imm_temp,main=paste0("t = 0"),asp = 1)
    # 
    # for(t in 1:15){
    # 
    #   mov_imm_v = mfraction_gg_ad %*% as.vector(t(mov_imm_temp)) # mrate_gg %*%
    #   mov_mat_v = mfraction_gg_ad %*% as.vector(t(mov_mat_temp))  # mrate_gg %*%
    # 
    #   mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
    #   mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
    # 
    #   # print(which((mov_imm_temp != mov_mat_temp_2)))
    # 
    #   mov_imm_temp = mov_imm_temp_2
    #   mov_mat_temp = mov_mat_temp_2
    # 
    #   plot(mov_mat_temp,main=paste0("t = ",t),asp = 1)
    # 
    # }
    # 
    # # and compare with
    # test = matrix(stationary_g_ad,nrow=40,ncol=40,byrow = T)
    # plot(test)
    # 
    # # 
    # ## Or
    # land_mask_na = land_mask
    # land_mask_na[which(land_mask_na == 0)] = NA
    # for(t in 1:(12*3)){
    #   if(t %in% c(1,13,25)){
    #     x11()
    #     par(mfrow = c(4,3))
    #   }
    #   
    #   plot(mat_N_at_Len[,,2,5,t] * land_mask_na,main=paste0(t),asp = 1)
    # }
    
  }

