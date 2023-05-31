## Morta at time step t
#----------------------
pars_Morta_setting_m_imm$year_LHP <- pars_Morta_setting_m_imm$Years_climsc[t]
#pars_Growth_setting_m_imm$clim_sc <- pars_Growth_setting_m_imm$clim_sc_test[c]

pars_Morta_setting_m_mat$year_LHP <- pars_Morta_setting_m_mat$Years_climsc[t]
#pars_Growth_setting_m_mat$clim_sc <- pars_Growth_setting_m_mat$clim_sc_test[c]

pars_Morta_setting_f_imm$year_LHP <- pars_Morta_setting_f_imm$Years_climsc[t]
#pars_Growth_setting_f_imm$clim_sc <- pars_Growth_setting_f_imm$clim_sc_test[c]

pars_Morta_setting_f_mat$year_LHP <- pars_Morta_setting_f_mat$Years_climsc[t]
#pars_Growth_setting_f_mat$clim_sc <- pars_Growth_setting_f_mat$clim_sc_test[c]

imm_male_M_size = array(0,dim = c(length(lat),length(lon),length(sizes)))
mat_male_size = array(0,dim = c(length(lat),length(lon),length(sizes)))
imm_fem_M_size = array(0,dim = c(length(lat),length(lon),length(sizes)))
mat_fem_M_size = array(0,dim = c(length(lat),length(lon),length(sizes)))

imm_male_M_mat <- mortality(E_morta_par = imm_male_M, # <-- could be size specific (but very long)
                            pars_Morta_setting_m_imm,
                            n_s,
                            M_spatial=F,
                            plot=FALSE)

mat_male_M_mat <- mortality(E_morta_par = mat_male_M, # <-- could be size specific (but very long)
                            pars_Morta_setting_m_mat,
                            n_s,
                            M_spatial=F,
                            plot=FALSE)

  
imm_fem_M_mat <- mortality(E_morta_par = imm_fem_M, # <-- could be size specific (but very long)
                           pars_Morta_setting_f_imm,
                           n_s,
                           M_spatial=F,
                           plot=FALSE)

mat_fem_M_mat <- mortality(E_morta_par = mat_fem_M, # <-- could be size specific (but very long)
                           pars_Morta_setting_f_mat,
                           n_s,
                           M_spatial=F,
                           plot=FALSE)

for(s in 1:n_p){
  
  # print(s)
    
  imm_male_M_size[,,s] = imm_male_M_mat
  mat_male_size[,,s] = mat_male_M_mat
  imm_fem_M_size[,,s] = imm_fem_M_mat
  imm_male_M_size[,,s] = mat_fem_M_mat
  
}
