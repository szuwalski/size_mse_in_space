##########################
## Save projection outputs
##########################

save_file_name = paste0("4_full_MSE/res/",simu_name,"/")
if(!dir.exists(save_file_name)) dir.create(save_file_name)

save(data=follow_ab_df,file = paste0(save_file_name,"follow_ab_df.RData"))

imm_N_at_Len_full[[t]] = temp_imm_N
save(data=imm_N_at_Len_full,file = paste0(save_file_name,"imm_N_at_Len_full.RData"))

mat_N_at_Len_full[[t]] = temp_mat_N
save(data=temp_mat_N,file = paste0(save_file_name,"temp_mat_N.RData"))

total_spatial_catch_full[[t]] = total_spatial_catch
save(data=total_spatial_catch,file = paste0(save_file_name,"total_spatial_catch.RData"))

catch_by_fisher_full[[t]] = catch_by_fisher
save(data=catch_by_fisher,file = paste0(save_file_name,"catch_by_fisher.RData"))
# 
# profit_by_fisher_full[[t]] = profit_by_fisher
# save(data=profit_by_fisher,file = paste0(save_file_name,"profit_by_fisher.RData"))
# 
# cost_by_fisher_full[[t]] = cost_by_fisher
# save(data=cost_by_fisher,file = paste0(save_file_name,"cost_by_fisher.RData"))
# 
# all_net_benefit_patch_full[[t]] = all_net_benefit_patch
# save(data=all_net_benefit_patch,file = paste0(save_file_name,"all_net_benefit_patch.RData"))
# 
# chosen_patch_full[[t]] = chosen_patch
# save(data=chosen_patch,file = paste0(save_file_name,"chosen_patch.RData"))
# 
# all_chosen_patch_full[[t]] = all_chosen_patch
# save(data=all_chosen_patch,file = paste0(save_file_name,"all_chosen_patch.RData"))
