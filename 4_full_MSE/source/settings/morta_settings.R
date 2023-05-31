# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Morta settings -------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
LHP = "morta"

# Spatial configuration  --------------------------------------------------
# -------------------------------------------------------------------------

# - Is growth spatially defined ? -----------------------------------------
M_spatial = F

# - Parameters for spatial field ------------------------------------------
scale_g <- 1.5 # scale of the RF
x_omega = 0.1 # how much of the total variance represent omega (in %)
x_epsilon = 0.05 # how much of the total variance represent epsilon (in %)

# Parameters --------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Morta : mean, sd parameter for morta parameter
morta_par <- array(0,c(2,4))
n_grpar_morta = 1 # number of parameter to define growth 

# --- Morta is define by one single parameter 
# ---- Female immature
morta_par[1,1] <- imm_fem_M
morta_par[2,1] <- morta_sd_imm

# ---- Female mature
morta_par[1,2] <- mat_fem_M
morta_par[2,2] <- morta_sd_mat

# ---- male immature
morta_par[1,3] <- imm_male_M
morta_par[2,3] <- morta_sd_imm

# ---- male mature
morta_par[1,4] <- mat_male_M
morta_par[2,4] <- morta_sd_mat


# Uncertainty/sensitivity -------------------------------------------------
# -------------------------------------------------------------------------
# -- Multiplicative or additive effect
ad_eff_morta <- TRUE # additive or multiplicative effect

# -- Which preferential habitat function 
pref_hab_morta= c("quadr") # peferential habitat function (logi, quad, lin)
pars_pref_hab_morta = c(4,2)



# Growth Settings ---------------------------------------------------------
# -------------------------------------------------------------------------
source("2_Max_spatial_projection/LHP_functions/pars_LHP_setting.R")
pars_Morta_setting_m_imm <- pars_LHP_f(pref_hab_morta,
                                       pars_pref_hab_morta,
                                       x_omega,
                                       x_epsilon,
                                       scale_g,
                                       LHP,
                                       ad_eff_morta,
                                       Morta_spatial,
                                       morta_par[,1],
                                       n_grpar_morta,
                                       Years_climsc)

pars_Morta_setting_m_mat <- pars_LHP_f(pref_hab_morta,
                                       pars_pref_hab_morta,
                                       x_omega,
                                       x_epsilon,
                                       scale_g,
                                       LHP,
                                       ad_eff_morta,
                                       Morta_spatial,
                                       morta_par[,2],
                                       n_grpar_morta,
                                       Years_climsc)

pars_Morta_setting_f_imm <- pars_LHP_f(pref_hab_morta,
                                       pars_pref_hab_morta,
                                       x_omega,
                                       x_epsilon,
                                       scale_g,
                                       LHP,
                                       ad_eff_morta,
                                       Morta_spatial,
                                       morta_par[,3],
                                       n_grpar_morta,
                                       Years_climsc)

pars_Morta_setting_f_mat <- pars_LHP_f(pref_hab_morta,
                                       pars_pref_hab_morta,
                                       x_omega,
                                       x_epsilon,
                                       scale_g,
                                       LHP,
                                       ad_eff_morta,
                                       Morta_spatial,
                                       morta_par[,4],
                                       n_grpar_morta,
                                       Years_climsc)
