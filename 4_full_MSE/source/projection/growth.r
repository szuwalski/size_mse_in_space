## Growth
#--------

# bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
for(x in 1:nrow(imm_N_at_Len[,,,,t])){
  for(y in 1:ncol(imm_N_at_Len[,,,,t]))
  {
    
    ## Cody's growth model
    #---------------------
    if(growth_model == "cody_model"){
      
      if(land_mask[x,y]!=0)
      {
        
        #======================================
        #==make size transition matrix immature
        size_transition_mat_m_imm<-matrix(ncol=length(sizes),nrow=length(sizes))
        size_transition_mat_f_imm<-matrix(ncol=length(sizes),nrow=length(sizes))
        
        for(n in 1:nrow(size_transition_mat_m_imm))
        {
          size_transition_mat_m_imm[n,]<-dnorm(x=sizes,mean=m_postmolt_imm[n],sd=growth_sd_imm[2])
          size_transition_mat_f_imm[n,]<-dnorm(x=sizes,mean=f_postmolt_imm[n],sd=growth_sd_imm[1])
        }
        
        #==ensure no crab shrinks in size after molting
        size_transition_mat_f_imm[lower.tri(size_transition_mat_f_imm)]<-0
        size_transition_mat_m_imm[lower.tri(size_transition_mat_m_imm)]<-0
        
        #==standardize to sum to 1
        for(z in 1:nrow(size_transition_mat_f_imm))
        {
          size_transition_mat_f_imm[z,]<-size_transition_mat_f_imm[z,]/sum(size_transition_mat_f_imm[z,],na.rm=T)
          size_transition_mat_m_imm[z,]<-size_transition_mat_m_imm[z,]/sum(size_transition_mat_m_imm[z,],na.rm=T)
        }
        
        #==immature crab molt, some mature, some remain immature
        if(molt_time[1,t]==1)
        {
          tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm
          temp_imm_N[x,y,1,]<-tmp_molt*term_molt_prob
          temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (1-term_molt_prob)*tmp_molt
          
        }
        if(molt_time[2,t]==1)
        {
          tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm
          temp_imm_N[x,y,2,]<-tmp_molt*term_molt_prob
          temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (1-term_molt_prob)*tmp_molt
        }
        
        if(terminal_molt==0)
        { 
          #======================================
          #==make size transition matrix mature
          size_transition_mat_m_mat<-matrix(ncol=length(sizes),nrow=length(sizes))
          size_transition_mat_f_mat<-matrix(ncol=length(sizes),nrow=length(sizes))
          
          for(n in 1:nrow(size_transition_mat_m_mat))
          {
            size_transition_mat_m_mat[n,]<-dnorm(x=sizes,mean=m_postmolt_mat[n],sd=growth_sd_mat[2])
            size_transition_mat_f_mat[n,]<-dnorm(x=sizes,mean=f_postmolt_mat[n],sd=growth_sd_mat[1])
          }
          
          #==ensure no crab shrinks in size after molting
          size_transition_mat_f_mat[lower.tri(size_transition_mat_f_mat)]<-0
          size_transition_mat_m_mat[lower.tri(size_transition_mat_m_mat)]<-0
          
          #==normalize
          for(z in 1:nrow(size_transition_mat_m_mat))
          {
            size_transition_mat_f_mat[z,]<-size_transition_mat_f_mat[z,]/sum(size_transition_mat_f_mat[z,],na.rm=T)
            size_transition_mat_m_mat[z,]<-size_transition_mat_m_mat[z,]/sum(size_transition_mat_m_mat[z,],na.rm=T)
          }
          if(!is.na(match(molt_time[1,t],t)) ) 
            temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat
          if(!is.na(match(molt_time[2,t],t)))
            temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat     
        }
        
      }
      
    }
    
    
    ## Maximes' growth model
    #-----------------------
    if(growth_model == "max_model"){
      
      if(land_mask[x,y]!=0)
      {
        
        #==immature crab molt, some mature, some remain immature
        if(molt_time[1,t]==1)
        {
          
          tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm[t,x,y,,] ##### Is the indexing right ??????
          temp_imm_N[x,y,1,]<-tmp_molt*term_molt_prob
          temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (1-term_molt_prob)*tmp_molt
          
        }
        
        if(molt_time[2,t]==1)
        {
          tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm[t,x,y,,] ##### Is the indexing right ??????
          temp_imm_N[x,y,2,]<-tmp_molt*term_molt_prob
          temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (1-term_molt_prob)*tmp_molt
        }
        
        if(terminal_molt==0)
        {
          
          if(!is.na(match(molt_time[1,t],t)))
            temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat[t,x,y,,] ##### Is the indexing right ??????
          if(!is.na(match(molt_time[2,t],t)))
            temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat[t,x,y,,] ##### Is the indexing right ??????
        }
        
      }
      
    }
    
  }
  
}

## Follow up chunk (only male at the moment)
follow_res = follow_ab_f(follow_ab_df,
                         simu_name,
                         phase = "post growth",
                         ab_imm_matrix = temp_imm_N[,,2,],
                         ab_mat_matrix = temp_mat_N[,,2,],
                         follow_ab_iter = follow_ab_iter,
                         t = t,
                         follow_ab = follow_ab)
follow_ab_df = follow_res[[1]]
follow_ab_iter = follow_res[[2]]


## Ontogenic movement
#--------------------
# that's rough, but cannot make better at the moment
# the best would be to make some kind of transition 
# movement from adult to juvenile with Jim's model 

# Take the size 
which_ad_size = which(sizes < ad_size)
which_ad_size = which_ad_size[length(which_ad_size)] + 1

# Compute last size at which we model ontogenic migration
if(length(dim(size_transition_mat_m_imm)) == 2){ # for Cody's model (size transition matrix dimensions = 2)
  last_size_mig = which(size_transition_mat_m_imm[which_ad_size,] < 0.001)[which_ad_size]
}else if(length(dim(size_transition_mat_m_imm)) == 5){ # for Max's model (size transition matrix dimensions = 5)
  last_size_mig2 = c()
  for(x in 1:length(lat)){
    for(y in 1:length(lon)){
      if(init_juv[x,y] > 0){
        last_size_mig2 = c(last_size_mig2,which(size_transition_mat_m_imm[t,x,y,which_ad_size,] < 0.001)[which_ad_size])
      }
    }
  }
  last_size_mig = max(last_size_mig2)
}

# For each size class where we apply ontogenic migration
for(s in which_ad_size:last_size_mig){
  # print(s)
  ontog_mig_N_imm_female = 0
  ontog_mig_N_imm_male = 0
  ontog_mig_N_mat_female = 0
  ontog_mig_N_mat_male = 0
  
  # For each cell of the grid in the juvenile habitats,
  # compute 
  for(x in 1:length(lon)){
    for(y in 1:length(lat)){
      if(init_juv[x,y] > 0){
        
        ontog_mig_N_imm_female = ontog_mig_N_imm_female + temp_imm_N[x,y,1,s]
        ontog_mig_N_imm_male = ontog_mig_N_imm_male + temp_imm_N[x,y,2,s]
        ontog_mig_N_mat_female = ontog_mig_N_mat_female + temp_mat_N[x,y,1,s]
        ontog_mig_N_mat_male = ontog_mig_N_mat_male + temp_mat_N[x,y,2,s]

        temp_imm_N[x,y,1,s] = 0
        temp_imm_N[x,y,2,s] = 0
        temp_mat_N[x,y,1,s] = 0
        temp_mat_N[x,y,2,s] = 0

      }
    }
  }

  Ab_ontog_mig_N_imm_female = sum(ontog_mig_N_imm_female)
  Ab_ontog_mig_N_imm_male = sum(ontog_mig_N_imm_male)
  Ab_ontog_mig_N_mat_female = sum(ontog_mig_N_mat_female)
  Ab_ontog_mig_N_mat_male = sum(ontog_mig_N_mat_male)

  Map_ontog_mig_N_imm_female = Ab_ontog_mig_N_imm_female * init_adult
  Map_ontog_mig_N_imm_male = Ab_ontog_mig_N_imm_male * init_adult
  Map_ontog_mig_N_mat_female = Ab_ontog_mig_N_mat_female * init_adult
  Map_ontog_mig_N_mat_male = Ab_ontog_mig_N_mat_male * init_adult

  temp_imm_N[,,1,s] = temp_imm_N[,,1,s] + Map_ontog_mig_N_imm_female
  temp_imm_N[,,2,s] = temp_imm_N[,,2,s] + Map_ontog_mig_N_imm_male
  temp_mat_N[,,1,s] = temp_mat_N[,,1,s] + Map_ontog_mig_N_mat_female
  temp_mat_N[,,2,s] = temp_mat_N[,,2,s] + Map_ontog_mig_N_mat_male

}

## Follow up chunk (only male at the moment)
follow_res = follow_ab_f(follow_ab_df,
                         simu_name,
                         phase = "post ontogenic migration",
                         ab_imm_matrix = temp_imm_N[,,2,],
                         ab_mat_matrix = temp_mat_N[,,2,],
                         follow_ab_iter = follow_ab_iter,
                         t = t,
                         follow_ab = follow_ab)
follow_ab_df = follow_res[[1]]
follow_ab_iter = follow_res[[2]]
