## Harvest Control Rule
#----------------------

## Implement HCR for next year quota
if(HCR == "cody"){
  
  f_mort = rep(0,length(sizes))
  temp_imm = apply(imm_N_at_Len[,,2,,t],c(3),sum,na.rm=T)
  temp_mat = apply(mat_N_at_Len[,,2,,t],c(3),sum,na.rm=T)
  fish_sel_hcr = fish_sel[2,]
  weight_at_size = wt_at_len[2,]
  functional_mat = 0
  mmb_def = 'morphometric'
  
  source("4_full_MSE/source/projection/hcr_cody.R")
  
}
