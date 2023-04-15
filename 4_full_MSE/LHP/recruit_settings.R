# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Morta settings -------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
LHP = "recruitment"

# Spatial configuration  --------------------------------------------------
# -------------------------------------------------------------------------

# - Is growth spatially defined ? -----------------------------------------
R_spatial = TRUE

# - Parameters for spatial field ------------------------------------------
scale_g <- 1.5 # scale of the RF
x_omega = 0.1 # how much of the total variance represent omega (in %)
x_epsilon = 0.05 # how much of the total variance represent epsilon (in %)


# Parameters --------------------------------------------------------------
# -------------------------------------------------------------------------
# -- Recruit : mean, sd parameter for morta parameter
recruit_par <- array(0,c(2,2))
n_grpar_recruit = 1 # number of parameter to define growth 

# --- Recruitment is define by one single parameter 
# ---- Female
recruit_par[1,1] <- meanlog_recruit_female
recruit_par[2,1] <- morta_sd_imm

# ---- Male
recruit_par[1,2] <- meanlog_recruit_male
recruit_par[2,2] <- morta_sd_mat


# Uncertainty/sensitivity -------------------------------------------------
# -------------------------------------------------------------------------
# -- Multiplicative or additive effect
ad_eff_recruit <- TRUE # additive or multiplicative effect

# -- Which preferential habitat function 
pref_hab_recruit= c("quadr") # peferential habitat function (logi, quad, lin)
pars_pref_hab_recruit = c(4,2)


# Growth Settings ---------------------------------------------------------
# -------------------------------------------------------------------------
source("2_Max_spatial_projection/LHP_functions/pars_LHP_setting.R")
pars_Recruit_setting_m <- pars_LHP_f(pref_hab_recruit,
                                     pars_pref_hab_recruit,
                                     x_omega,
                                     x_epsilon,
                                     scale_g,
                                     LHP,
                                     ad_eff_recruit,
                                     Recruit_spatial,
                                     recruit_par[,1],
                                     n_grpar_recruit,
                                     Years_climsc)

pars_Recruit_setting_f <- pars_LHP_f(pref_hab_recruit,
                                     pars_pref_hab_recruit,
                                     x_omega,
                                     x_epsilon,
                                     scale_g,
                                     LHP,
                                     ad_eff_recruit,
                                     Recruit_spatial,
                                     recruit_par[,2],
                                     n_grpar_recruit,
                                     Years_climsc)

