## Growth at time step t
#-----------------------
pars_Growth_setting_m_imm$year_LHP <- pars_Growth_setting_m_imm$Years_climsc[t]
#pars_Growth_setting_m_imm$clim_sc <- pars_Growth_setting_m_imm$clim_sc_test[c]

pars_Growth_setting_m_mat$year_LHP <- pars_Growth_setting_m_mat$Years_climsc[t]
#pars_Growth_setting_m_mat$clim_sc <- pars_Growth_setting_m_mat$clim_sc_test[c]

pars_Growth_setting_f_imm$year_LHP <- pars_Growth_setting_f_imm$Years_climsc[t]
#pars_Growth_setting_f_imm$clim_sc <- pars_Growth_setting_f_imm$clim_sc_test[c]

pars_Growth_setting_f_mat$year_LHP <- pars_Growth_setting_f_mat$Years_climsc[t]
#pars_Growth_setting_f_mat$clim_sc <- pars_Growth_setting_f_mat$clim_sc_test[c]

growth_m_imm[t,,,,] <- growth(sizes,
                              binclass,
                              pars_Growth_setting_m_imm,
                              n_s,
                              n_p,
                              G_spatial,
                              plot = F)

growth_m_mat[t,,,,] <- growth(sizes,
                              binclass,
                              pars_Growth_setting_m_mat,
                              n_s,
                              n_p,
                              G_spatial,
                              plot = F)

growth_f_imm[t,,,,] <- growth(sizes,
                              binclass,
                              pars_Growth_setting_f_imm ,
                              n_s,
                              n_p,
                              G_spatial,
                              plot = F)

growth_f_mat[t,,,,] <- growth(sizes,
                              binclass,
                              pars_Growth_setting_f_mat,
                              n_s,
                              n_p,
                              G_spatial,
                              plot = F)

size_transition_mat_f_imm  <- growth_f_imm
size_transition_mat_f_mat  <- growth_f_mat
size_transition_mat_m_imm  <- growth_m_imm
size_transition_mat_m_mat  <- growth_m_mat
