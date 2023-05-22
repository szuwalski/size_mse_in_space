## Fishery projection
#--------------------

quota_remaining <- quota[f]

#==THISNEEDS TO BE FIXED====
#==calculate net benefits by patch
#==this needs to be based on the amount of quota available
# if there is a patch nearby they can get their quota filled, they will go there
#==there should be some relationship between cost and biomass in a cell?

temp_catch<-array(dim=c(length(lat),length(lon),length(sizes)))
for(sex in 1:2)
  for(x in 1:length(sizes))
  {
    temp_catch[,,x]<-temp_imm_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
      temp_mat_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
  }

## Choose fishing location
catch_patch<-apply(temp_catch,c(1,2),sum)
catch_patch[catch_patch>quota[f]]<-quota[f] # this makes it so they don't travel a long way if they can get it close < # We could modify this
net_benefit_patch<-catch_patch*price-cost_patch
all_net_benefit_patch[,,t,f] = net_benefit_patch

if(fishing_process == "max_benefit_min_dist"){
  
  max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
  chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
  
}else if(fishing_process == "stochastic"){
  
  prob_net = net_benefit_patch[which(!is.na(net_benefit_patch) & net_benefit_patch > 0)]
  chosen_patch<-which(net_benefit_patch == sample(prob_net,size=1,prob=prob_net),arr.ind=T)
  chosen_patch=chosen_patch[sample(1:nrow(chosen_patch),1),]
  
}

#========================================================
#==subtract catch from locations while quota is remaining
quota_it = 1
while(quota_it < max_quota_it & quota_remaining>0.1 & net_benefit_patch[chosen_patch[1],chosen_patch[2]]>0)
{
  
  quota_it = quota_it + 1
  
  if(fishing_process == "max_benefit_min_dist"){
    
    max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
    chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
    
  }else if(fishing_process == "stochastic"){
    
    prob_net = net_benefit_patch[which(!is.na(net_benefit_patch) & net_benefit_patch > 0)]
    chosen_patch<-which(net_benefit_patch == sample(prob_net,size=1,prob=),arr.ind=T)
    chosen_patch=chosen_patch[sample(1:nrow(chosen_patch),1),]
    
  }
  
  if(is.array(chosen_patch)) all_chosen_patch[,t,f] = chosen_patch[1,]
  if(is.integer(chosen_patch) & length(chosen_patch) == 2) all_chosen_patch[,t,f] = chosen_patch
  
  #==calculate total potential catch in a patch
  potential_catch<-0
  for(sex in 1:2)
    for(x in 1:length(sizes))
    {
      if(print_messages) print(paste0("sex=",sex," | sizes=",x," | ",
                                      " | temp_imm_N=",temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x],
                                      " | temp_mat_N=",temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x],
                                      " | potential catch= ",potential_catch))
      potential_catch<-potential_catch + temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
        temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
    }
  
  #==patch has less than needed to fill quota
  if(potential_catch<=quota_remaining)
  {
    
    if(print_messages) print(paste0("potential_catch<=quota_remaining | potential_catch=",potential_catch,"| quota_remaining=",quota_remaining,"| use_harv=",use_harv))
    
    for(sex in 1:2)
      for(x in 1:length(sizes))
      {
        
        total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
          temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
        
        total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] + 
          temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x] + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
        
        temp_catch <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
        
        catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
        
      }
    cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]   <- cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f] + cost_patch[chosen_patch[1],chosen_patch[2]]        
    quota_remaining<-quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
    #==update temp array of n at len
    for(sex in 1:2)
      for(x in 1:length(sizes))
      {
        temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
        temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
      }
    
  }
  
  #==patch has more than needed to fill quota
  if(potential_catch>quota_remaining)
  {
    
    #==find harvest rate that would fill quota
    maxHarv<-1
    minHarv<-.0000001
    for(o in 1:25)
    {
      use_harv<-(maxHarv+minHarv)/2
      temp_cat<-0
      for(sex in 1:2)
        for(x in 1:length(sizes))
        {
          temp_cat<-temp_cat+
            temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
            temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
        }
      if(temp_cat<quota_remaining)
        minHarv<-use_harv
      if(temp_cat>quota_remaining)
        maxHarv<-use_harv
    }
    
    if(print_messages) print(paste0("potential_catch<=quota_remaining | potential_catch=",potential_catch,"| quota_remaining=",quota_remaining,"| use_harv=",use_harv))
    
    temp_catch<-0
    for(sex in 1:2)
      for(x in 1:length(sizes))
      {
        total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
          temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv  
        
        total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch_nb[chosen_patch[1],chosen_patch[2],sex,x,t] + 
          temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x] + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]
        
        temp_catch<- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
          temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
        
        catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
      }
    #sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
    quota_remaining <- quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
    cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]   <- cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f] + cost_patch[chosen_patch[1],chosen_patch[2]]
    
    #==update temp array of n at len
    for(sex in 1:2)
      for(x in 1:length(sizes))
      {
        temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv
        temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] <- temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x] - temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv
      }
  }
  
  #===update the net benefits in while loop
  temp_catch<-array(dim=c(length(lat),length(lon),length(sizes)))
  for(sex in 1:2)
    for(x in 1:length(sizes))
      temp_catch[,,x]<-temp_imm_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + temp_mat_N[,,sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
  
  catch_patch<-apply(catch_patch,c(1,2),sum)
  catch_patch[catch_patch>quota_remaining]<-quota_remaining
  net_benefit_patch<-catch_patch*price - cost_patch
}

profit_by_fisher[chosen_patch[1],chosen_patch[2],t,f] <- sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])*price - cost_by_fisher[chosen_patch[1],chosen_patch[2],t,f]

