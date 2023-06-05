###############################################
## Stocking and debugging param for projections
###############################################
list_it = 0

## Objects to stock variables
#----------------------------
# Abundance list
imm_N_at_Len_full = list()
mat_N_at_Len_full = list()

# Fishery list
total_spatial_catch_full = list()
catch_by_fisher_full = list()
profit_by_fisher_full = list()
cost_by_fisher_full = list()
all_net_benefit_patch_full = list()
chosen_patch_full = list()
all_chosen_patch_full = list()

# Time-steps array
total_spatial_catch<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
total_spatial_catch_nb<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period)))
catch_by_fisher<-array(0,dim=c(length(lat),length(lon),sexN,length(sizes),length(proj_period),fishers))
profit_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
cost_by_fisher<-array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
all_net_benefit_patch=array(0,dim=c(length(lat),length(lon),length(proj_period),fishers))
all_chosen_patch=array(0,dim=c(2,length(proj_period),fishers))

## Quota vector
quota_vec = c()
in_f_vec = c()


## For debugging
#---------------
## Fix some values so that the loop works fine
max_quota_it = 100 # maximum iteration for filling the quota (so that we are not blocked in a while loop)
use_harv = 1

## Following abundance throughout the loop
follow_ab = T
follow_ab_iter = 1
follow_ab_df = data.frame()

## Print distribution throughout the loop
print_distrib = F

## Print messages for identifying where are bugs
print_messages <- F

## Function for following abundance throughout the MSE
follow_ab_f = function(follow_ab_df, # data frame of abundance
                       simu_name, # Simulation name
                       phase, # phase of the MSE
                       ab_imm_matrix, # spatial matrix of abundance at size (immature)
                       ab_mat_matrix, # spatial matrix of abundance at size (mature)
                       follow_ab_iter, # line of the data frame to be filled
                       t, # Time step
                       follow_ab = F){
  
  if(follow_ab){
    
    for(size in 1:dim(ab_imm_matrix)[3]){
      
      follow_ab_df[follow_ab_iter,"simu_name"] = simu_name
      follow_ab_df[follow_ab_iter,"phase"] = phase
      follow_ab_df[follow_ab_iter,"size"] = size
      follow_ab_df[follow_ab_iter,"ab_imm"] = sum(ab_imm_matrix[,,size]) 
      follow_ab_df[follow_ab_iter,"ab_mat"] = sum(ab_mat_matrix[,,size])
      follow_ab_df[follow_ab_iter,"t"] = t
      follow_ab_df[follow_ab_iter,"iter"] = follow_ab_iter
      follow_ab_iter=follow_ab_iter+1
      
    }
    
  }
  
  res = list(follow_ab_df,follow_ab_iter)
  
  return(res)
  
}
