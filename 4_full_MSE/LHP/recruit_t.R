## Recruit at time step t
#-------------------------
pars_Recruit_setting_m$year_LHP <- pars_Recruit_setting_m$Years_climsc[t]
#pars_Growth_setting_m_imm$clim_sc <- pars_Growth_setting_m_imm$clim_sc_test[c]

pars_Recruit_setting_f$year_LHP <- pars_Recruit_setting_f$Years_climsc[t]
#pars_Growth_setting_f_imm$clim_sc <- pars_Growth_setting_f_imm$clim_sc_test[c]

# Simulate potential recruitment on the full area
logrecruit_male_mat <- recruitment(meanlog_recruit_male,
                                   pars_Recruit_setting_m,
                                   n_s,
                                   R_spatial)

logrecruit_fem_mat <- recruitment(meanlog_recruit_female,
                                  pars_Recruit_setting_f,
                                  n_s,
                                  R_spatial)

# Pass log recruitment to standard recruitment and allocate 
# proportionally the recruitment in Juveniles habitats 
recruit_male_mat = exp(logrecruit_male_mat) * init_juv
recruit_fem_mat = exp(logrecruit_fem_mat) * init_juv
# --> Imply that temperature affect recruitment at the cell level

# Aggregate to model stochasticity
meanlog_recruit_male_2 = sum(recruit_male_mat)
meanlog_recruit_female_2 = sum(recruit_fem_mat)
# Modeling recruitment at aggregated level is equivalent to assume 
# that this random noise accounts for all the processes unmodelled 
# that impact the recruitment from spawning to settlement and being recruited