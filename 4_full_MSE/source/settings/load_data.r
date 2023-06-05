## Load data
#-----------

load("4_full_MSE/data/Snow_GMACS_repfile.Rdata") # Outputs provided by Mathieu from GMACS
load("4_full_MSE/data/Abundance_Recruitment_Assessment.Rdata")
load("4_full_MSE/data/NRS_vars_wide_op.Rdata") # Environmental covariates --> all variables
load("2_Max_spatial_projection/Climate_scenarios/mn_var_all.Rdata") # --> bottom temperature
load("4_full_MSE/data/m_gam_dat.RData") # gam for mortality

## Cold pool extent
scenarios_vec = c("ACLIMregion_B10K-K20P19_CMIP5_GFDL_rcp45",
                  "ACLIMregion_B10K-K20P19_CMIP5_GFDL_rcp85",
                  "ACLIMregion_B10K-K20P19_CMIP6_gfdl_ssp126",
                  "ACLIMregion_B10K-K20P19_CMIP6_gfdl_ssp585")

proj_cold_pool_df = NRS_vars_wide_op %>%
  filter(sim %in% scenarios_vec) %>% 
  dplyr::select(year,sim,fracbelow2)

proj_cold_pool_df$simulation = NA
proj_cold_pool_df$simulation[which(proj_cold_pool_df$sim == "ACLIMregion_B10K-K20P19_CMIP5_GFDL_rcp45")] = "rcp45"
proj_cold_pool_df$simulation[which(proj_cold_pool_df$sim == "ACLIMregion_B10K-K20P19_CMIP5_GFDL_rcp85")] = "rcp85"
proj_cold_pool_df$simulation[which(proj_cold_pool_df$sim == "ACLIMregion_B10K-K20P19_CMIP6_gfdl_ssp126")] = "ssp126"
proj_cold_pool_df$simulation[which(proj_cold_pool_df$sim == "ACLIMregion_B10K-K20P19_CMIP6_gfdl_ssp585")] = "ssp585"

hind_cold_pool_df = NRS_vars_wide_op %>% 
  filter(sim == "ACLIM + Operational Hindcast")

# ## Check consistency with data used to fit GAM
# m_gam_dat$year = m_gam_dat$Year
# test = inner_join(m_gam_dat,hind_cold_pool_df)
# plot(lm(test$coldpool~test$fracbelow2))


