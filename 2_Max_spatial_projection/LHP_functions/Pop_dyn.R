
Function <- DynPop(lat,lon,
                   n_lat,n_lon,n_n,n_p,n_t,rec_sizes,
                   imm_N_at_Len,mat_N_at_Len,
                   quota,fishers,cost_patch,price,
                   size_transition_mat_f_imm,size_transition_mat_f_mat,size_transition_mat_m_imm,size_transition_mat_m_mat,
                   survey_time,  
                   fish_time,		
                   recruit_time,	
                   move_time,		  
                   molt_time_month_m,
                   molt_time_month_f,
                   molt_time_m,
                   molt_time_f,
                   molt_time,
                   mate_time,
                   term_molt_prob,
                   terminal_molt,
                   fish_50)

# -------------------------------------------------------------------------
# 0- Some good questions to have in mind ----------------------------------
# -------------------------------------------------------------------------

#==what is the first size class to be modeled?
#==this depend on how size at maturity changes
#==do crab mature after a set number of molts?
#==do they molt no matter what, just different increments?
#==OR do they molt the same size increment, but fewer times?
#==size dependent molting?
#==or start the model at the point that they are already only molting once a year
#==and the size they enter the model change based on the temperature during the time period

# -------------------------------------------------------------------------
# 1- Projection------------------------------------------------------------
# -------------------------------------------------------------------------

total_spatial_catch<-array(0,dim=c(n_lat,n_lon,n_n,n_p,n_t))
catch_by_fisher<-array(0,dim=c(n_lat,n_lon,n_n,n_p,n_t,fishers))
profit_by_fisher<-array(0,dim=c(n_lat,n_lon,n_t,fishers))
cost_by_fisher <-array(0,dim=c(n_lat,n_lon,n_t,fishers))

#==indices: lat,lon,sex,size,time
for(t in 1:(n_t-1))
  #for(t in 1:320)
{
  #==create a 'working' array for a given time step of both mature and immature critters
  
  temp_imm_N<-imm_N_at_Len[,,,,t]
  temp_mat_N<-mat_N_at_Len[,,,,t]
  #filled.contour(x=lon,y=rev(lat),g(temp_mat_N[,,1,5]),plot.axes=map(add=TRUE,fill=T,col='grey') )
  # if(survey_time[t]==1)
  #   collect_survey_data()
  
  #==========================
  #==FISHERY OCCURS
  #==========================
  #==indices: lat,lon,sex,size,time
  
  if(fish_time[t]==1)
  {
    for(f in 1:fishers)
    {
      quota_remaining <- quota[f,t]
      
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
      
      catch_patch<-apply(temp_catch,c(1,2),sum)
      catch_patch[catch_patch>quota[f,t]]<-quota[f,t] # this makes it so they don't travel a long way if they can get it close
      
      #filled.contour(x=lon,y=rev(lat),g(catch_patch),plot.axes=map(add=TRUE,fill=T,col='grey') )
      net_benefit_patch < -catch_patch*price-cost_patch
      
      #filled.contour(x=lon,y=rev(lat),g(net_benefit_patch),plot.axes=map(add=TRUE,fill=T,col='grey'),zlim=c(0,max(net_benefit_patch,na.rm=T)) )
      max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
      chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
      # filled.contour(x=lon,y=rev(lat),g(net_benefit_patch),
      #               plot.axes=c(map(add=TRUE,fill=T,col='grey'),
      #                          points(x=lon[chosen_patch[2]],y=lat[chosen_patch[1]],col=2,pch=16)),
      #                                zlim=c(0,max(net_benefit_patch,na.rm=T)) )
      # 
      #========================================================
      #==subtract catch from locations while quota is remaining
      while(quota_remaining>0.1 & net_benefit_patch[chosen_patch[1],chosen_patch[2]]>0)
      {
        #==find closest, highest value, fishable patch
        max_net_benefit<-which(net_benefit_patch==max(net_benefit_patch,na.rm=T),arr.ind=T)
        #==have to do this if there are two patches with identical net benefits
        chosen_patch<-max_net_benefit[which(distance_map[max_net_benefit]==min(distance_map[max_net_benefit])),]
        
        #==calculate total potential catch in a patch
        potential_catch<-0
        for(sex in 1:2)
          for(x in 1:length(sizes))
          {
            potential_catch<-potential_catch + temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
              temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]
          }
        
        #==patch has less than needed to fill quota
        if(potential_catch<=quota_remaining)
        { 
          for(sex in 1:2)
            for(x in 1:length(sizes))
            {
              total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x] + 
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]       
              
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
          
          temp_catch<-0
          for(sex in 1:2)
            for(x in 1:length(sizes))
            {
              total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] <- total_spatial_catch[chosen_patch[1],chosen_patch[2],sex,x,t] + 
                temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv + 
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*use_harv  
              
              
              temp_catch<- temp_imm_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv + 
                temp_mat_N[chosen_patch[1],chosen_patch[2],sex,x]*fish_sel[sex,x]*wt_at_len[sex,x]*use_harv
              
              catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f]  <- catch_by_fisher[chosen_patch[1],chosen_patch[2],sex,x,t,f] + temp_catch
            }
          #sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
          quota_remaining<-quota_remaining - sum(catch_by_fisher[chosen_patch[1],chosen_patch[2],,,t,f])
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
    }
    
  }
  
  #==========================
  #==MOVEMENT OCCURS
  #==========================  
  #==movement input as a .csv?
  #==movement constant?
  #==movement follows gradient?
  
  if(move_time[t]==1)
  {
    #==two options: gaussian and temperature mediated
    #==gaussian
    #==create disperal kernel for each space 
    
  }
  
  #==========================
  #==GROWTH OCCURS
  #==========================
  if(molt_time[1,t]==1 | molt_time[2,t]==1 )
  {
    # bot_temp_dat<-read.csv(paste("temp_data/bot_temp_",time,".csv",sep=""),header=T)
    for(x in 1:nrow(imm_N_at_Len[,,,,t]))
      for(y in 1:ncol(imm_N_at_Len[,,,,t]))
      {
        if(land_mask[x,y]!=0)
        {
          #==immature crab molt, some mature, some remain immature
          if(molt_time[1,t]==1)
          {
            tmp_molt          <-temp_imm_N[x,y,1,]%*%size_transition_mat_f_imm[t,x,y,,]
            temp_imm_N[x,y,1,]<-tmp_molt*term_molt_prob
            temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,] + (1-term_molt_prob)*tmp_molt
            
          }
          if(molt_time[2,t]==1)
          {
            tmp_molt          <-temp_imm_N[x,y,2,]%*%size_transition_mat_m_imm[t,x,y,,]
            temp_imm_N[x,y,2,]<-tmp_molt*term_molt_prob
            temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,] + (1-term_molt_prob)*tmp_molt
          }
          
          if(terminal_molt==0)
          { 
            if(!is.na(match(molt_time[1,t],t)) ) 
              temp_mat_N[x,y,1,]<-temp_mat_N[x,y,1,]%*%size_transition_mat_f_mat[t,x,y,,]
            if(!is.na(match(molt_time[2,t],t)))
              temp_mat_N[x,y,2,]<-temp_mat_N[x,y,2,]%*%size_transition_mat_m_mat[t,x,y,,]    
          }
        }
      }  
  } 
  
  #==========================
  #==SPAWNING OCCURS
  #==========================      
  #==this makes a map of spawning biomass to be used with transition matrices for recruitment
  
  if(mate_time[t]==1 )
  {
    #==aggregate spawnign biomass
    #==just count female biomass?
    #==include sperm reserves?
    #==include biennial spawning?
    
    #==THIS IS JUST NUMBERS RIGHT NOW...
    spbiom<-apply(temp_mat_N[,,1,],c(1,2),sum,na.rm=T)
    plot_spb<-spbiom
    plot_spb[plot_spb==0]<-NA
    # filled.contour(x=lon,y=rev(lat),g(plot_spb),plot.axes=map(add=TRUE,fill=T,col='grey') )
  }
  
  
  #==========================
  #==RECRUITMENT OCCURS
  #========================== 
  #==how do we determine which bins they drop into?
  #==will temperature determine the size they reach in the time before they settle?
  if(recruit_time[t]==1)
  {
    #==this is a dumb temporary fix
    #==this implants the original recruitment with some error
    #==ultimately we need an algorithm to determine location and intensity of recruitment
    #==teleconnections postdoc will hopefully help with this
    tmp_rec_1<- matrix(rnorm(length(imm_N_at_Len[,,1,1,1]),1,1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,1,1,1]))
    tmp_rec_1[tmp_rec_1<0]<-0
    tmp_rec_2<- matrix(rnorm(length(imm_N_at_Len[,,2,1,1]),1,1),ncol=ncol(imm_N_at_Len[,,1,1,1]),nrow=nrow(imm_N_at_Len[,,2,1,1]))
    tmp_rec_2[tmp_rec_2<0]<-0  
    
    for(r in 1:rec_sizes)
    {
      temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + imm_N_at_Len[,,1,r,1]*tmp_rec_1
      temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + imm_N_at_Len[,,2,r,1]*tmp_rec_2
    }
  }
  
  #==update dynamics
  imm_N_at_Len[,,1,,t+1] <-  temp_imm_N[,,1,]*exp(-imm_fem_M*1/year_step)
  imm_N_at_Len[,,2,,t+1] <-  temp_imm_N[,,2,]*exp(-imm_male_M*1/year_step)
  mat_N_at_Len[,,1,,t+1] <-  temp_mat_N[,,1,]*exp(-mat_fem_M*1/year_step)
  mat_N_at_Len[,,2,,t+1] <-  temp_mat_N[,,2,]*exp(-mat_male_M*1/year_step)
  
  
}

tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len,c(5),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len,c(5),sum,na.rm=T)

par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,.1,1,1))
plot(tot_imm[100:length(tot_imm)],type='l',las=1,xaxt='n',ylim=c(0,9000000000))
lines(tot_mat[100:length(tot_imm)],lty=2)
legend('topright',bty='n',lty=c(1,2),legend=c("Immature N","Mature N"))
# plot(tot_catch,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total catch"))
# plot(tot_cost,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total cost"))
# plot(tot_profit,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total profits"))
plot(tot_catch[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16,ylim=c(0,60000000))
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit[(tot_profit>0)],las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total profits"))
