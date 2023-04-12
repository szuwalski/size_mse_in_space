# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Growth settings ------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
LHP = "growth"

# Spatial configuration  --------------------------------------------------
# -------------------------------------------------------------------------

# - Is growth spatially defined ? -----------------------------------------
Growth_spatial = TRUE

# - Parameters for spatial field ------------------------------------------
scale_g <- 1.5 # scale of the RF
x_omega = 0.1 # how much of the total variance represent omega (in %)
x_epsilon = 0.05 # how much of the total variance represent epsilon (in %)


# Parameters --------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Growth : mean, sd and scale parameter for growth transition matrix
growth_par_beta <- array(0,c(2,3,4))
n_grpar_growth = dim(growth_par_beta)[2] # number of parameter to define growth 

# --- Growth is define by 2 parameters : initial and last size increment 
# n.b. I divide the size by 10 to convert mm into cm --> required to make the transition between GMACS model and Max model
# ---- Female immature
growth_par_beta[1,1,1] <- f_postmolt_imm[1] / 10 # initial size increment  
growth_par_beta[1,2,1] <- f_postmolt_imm[length(f_postmolt_imm)] / 10 # last size increment : mean
growth_par_beta[1,3,1] <- scale_f_imm # scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,1] <- growth_sd_imm[1] / 10 # initial size increment  : sd
growth_par_beta[2,2,1] <- growth_sd_imm[1] / 10 # last size increment : sd
growth_par_beta[2,3,1] <- NA

# ---- Female mature
growth_par_beta[1,1,2] <- f_postmolt_mat[1] / 10 # initial size increment  
growth_par_beta[1,2,2] <- f_postmolt_mat[length(f_postmolt_mat)] / 10 # last size increment : mean
growth_par_beta[1,3,2] <- scale_f_mat # scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,2] <- growth_sd_mat[1] / 10 # initial size increment  : sd
growth_par_beta[2,2,2] <- growth_sd_mat[1] / 10 # last size increment : sd
growth_par_beta[2,3,2] <-  NA 

# ---- male immature
growth_par_beta[1,1,3] <- m_postmolt_imm[1] / 10 # initial size increment  
growth_par_beta[1,2,3] <- m_postmolt_imm[length(m_postmolt_imm)] / 10 # last size increment : mean
growth_par_beta[1,3,3] <- scale_m_imm # scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,3] <- growth_sd_imm[2] / 10 # initial size increment  : sd
growth_par_beta[2,2,3] <- growth_sd_imm[2] / 10 # last size increment : sd
growth_par_beta[2,3,3] <-  NA 

# ---- male mature
growth_par_beta[1,1,4] <- m_postmolt_mat[1] / 10 # initial size increment  
growth_par_beta[1,2,4] <- m_postmolt_mat[length(m_postmolt_mat)] / 10 #last size increment : mean
growth_par_beta[1,3,4] <- scale_m_mat # scale parameter of the gamma distribution used to generate the growth transition matrix 

growth_par_beta[2,1,4] <- growth_sd_mat[2] / 10 # initial size increment: sd
growth_par_beta[2,2,4] <- growth_sd_mat[2] / 10 # last size increment: sd
growth_par_beta[2,3,4] <-  NA 


# Uncertainty/sensitivity -------------------------------------------------
# -------------------------------------------------------------------------
# -- Multiplicative or additive effect
ad_eff_growth <- TRUE # additive or multiplicative effect

# -- Which preferential habitat function 
pref_hab_growth= c("quadr") # peferential habitat function (logi, quad, lin)
pars_pref_hab_growth = c(4,2)



# Growth Settings ---------------------------------------------------------
# -------------------------------------------------------------------------
source("2_Max_spatial_projection/LHP_functions/pars_LHP_setting.R")
pars_Growth_setting_m_imm <- pars_LHP_f(pref_hab_growth,
                                        pars_pref_hab_growth,
                                        x_omega,
                                        x_epsilon,
                                        scale_g,
                                        LHP,
                                        ad_eff_growth ,
                                        Growth_spatial,
                                        growth_par_beta[,,1],
                                        n_grpar_growth,
                                        Years_climsc)

pars_Growth_setting_m_mat <- pars_LHP_f(pref_hab_growth,
                                        pars_pref_hab_growth,
                                        x_omega,
                                        x_epsilon,
                                        scale_g,
                                        LHP,
                                        ad_eff_growth ,
                                        Growth_spatial,
                                        growth_par_beta[,,2],
                                        n_grpar_growth,
                                        Years_climsc)

pars_Growth_setting_f_imm <- pars_LHP_f(pref_hab_growth,
                                        pars_pref_hab_growth,
                                        x_omega,
                                        x_epsilon,
                                        scale_g,
                                        LHP,
                                        ad_eff_growth ,
                                        Growth_spatial,
                                        growth_par_beta[,,3],
                                        n_grpar_growth,
                                        Years_climsc)

pars_Growth_setting_f_mat <- pars_LHP_f(pref_hab_growth,
                                        pars_pref_hab_growth,
                                        x_omega,
                                        x_epsilon,
                                        scale_g,
                                        LHP,
                                        ad_eff_growth,
                                        Growth_spatial,
                                        growth_par_beta[,,4],
                                        n_grpar_growth,
                                        Years_climsc)
