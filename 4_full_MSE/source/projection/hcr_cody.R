############
## Apply HCR
############

#==project foward under fishing at fmsyproxy
in_f<-fmsy_proxy
proj_mmb<-proj_to_fish(temp_imm,
                       temp_mat,
                       fish_sel_hcr,
                       weight_at_size,
                       in_f,
                       functional_mat,
                       mmb_def)

#==if MMB > bmsy_proxy you're good otherwise jump into the loop to find F on HCR
if(proj_mmb < bmsy_proxy)
{
  in_f<-0
  proj_mmb<-proj_to_fish(temp_imm,
                         temp_mat,
                         fish_sel_hcr,
                         weight_at_size,
                         in_f,
                         functional_mat,
                         mmb_def)
  
  if(proj_mmb > b_hcr*bmsy_proxy)
  {
    proj_mmb<-proj_to_fish(temp_imm,
                           temp_mat,
                           fish_sel_hcr,
                           weight_at_size,
                           in_f=fmsy_proxy*0.5,
                           functional_mat,
                           mmb_def)
    
    for(z in 1:15)
    {
      in_f<-fmsy_proxy * (((proj_mmb/bmsy_proxy)-a_hcr)/(1-a_hcr))
      proj_mmb<-proj_to_fish(temp_imm,
                             temp_mat,
                             fish_sel_hcr,
                             weight_at_size,
                             in_f,
                             functional_mat,
                             mmb_def)
    }
  }
}

f_mort<-in_f

## Compute quota
quota_nb = f_mort / (f_mort + m) * (1 - exp(-(f_mort + m))) * n_at_size * fish_sel_hcr
quota_biom = sum(weight_at_size * quota_nb)
quota<-rep(quota_biom/fishers,fishers)

in_f_vec = c(in_f_vec,in_f)
quota_vec = c(quota_vec,quota_biom)
