##################################
## Function to be used for the HCR
##################################

proj_to_fish<-function(temp_imm,
                       temp_mat,
                       fish_sel,
                       weight_at_size,
                       in_f,
                       functional_mat,
                       mmb_def)
{
  
  temp_imm2 <- temp_imm * exp(-(in_f)*fish_sel)
  temp_mat2 <- temp_mat * exp(-(in_f)*fish_sel)
  
  if(mmb_def == 'morphometric')
    proj_mmb <- sum(temp_mat2 * weight_at_size)
  if(mmb_def == 'functional')
    proj_mmb <- sum(temp_mat2*weight_at_size*functional_mat) + sum(temp_imm2*weight_at_size*functional_mat)
  return(proj_mmb)
  
}