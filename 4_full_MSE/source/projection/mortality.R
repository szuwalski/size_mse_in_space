## Natural mortality
#-------------------

if(morta_model == "cody_model"){
  
  print("Cody morta")
  imm_N_at_Len[,,1,,t+1] <-  temp_imm_N[,,1,]*exp(-imm_fem_M*1/year_step)
  imm_N_at_Len[,,2,,t+1] <-  temp_imm_N[,,2,]*exp(-imm_male_M*1/year_step)
  mat_N_at_Len[,,1,,t+1] <-  temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
  mat_N_at_Len[,,2,,t+1] <-  temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
  
}else if(morta_model == "max_model"){
  
  print("Max morta")
  imm_N_at_Len[,,1,,t+1] <- temp_imm_N[,,1,]*exp(-imm_fem_M_size*1/year_step)
  imm_N_at_Len[,,2,,t+1] <- temp_imm_N[,,2,]*exp(-imm_male_M_size*1/year_step)
  mat_N_at_Len[,,1,,t+1] <- temp_mat_N[,,1,]*exp(-mat_fem_M_size*1/year_step)
  mat_N_at_Len[,,2,,t+1] <- temp_mat_N[,,2,]*exp(-mat_male_size*1/year_step)
  
}else if(morta_model == "coldpool"){
  
  df = data.frame(coldpool = proj_cold_pool_df_3$fracbelow2_scaled)
  
  imm_fem_M_size_cp = exp(predict(mod,df))
  imm_male_M_size_cp = exp(predict(mod,df))
  
  imm_N_at_Len[,,1,,t+1] <- temp_imm_N[,,1,] * exp(-as.numeric(imm_fem_M_size_cp)*1/year_step)
  imm_N_at_Len[,,2,,t+1] <- temp_imm_N[,,2,]*exp(-as.numeric(imm_male_M_size_cp)*1/year_step)
  mat_N_at_Len[,,1,,t+1] <- temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
  mat_N_at_Len[,,2,,t+1] <- temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
  
}
