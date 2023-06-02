## Plot projection
#-----------------

## Check plot
vec_fact = c("init","post fishery","post growth", "post ontogenic migration","post monthly movement","post recruitment","post natural mortality")

follow_ab_df$phase = factor(follow_ab_df$phase,levels = vec_fact)

test = follow_ab_df %>% 
  filter(size > 14) %>% 
  group_by(phase,t) %>% 
  dplyr::summarise(ab_male_mat  = ab_male_mat,
                   ab_male_imm = ab_male_imm) %>% 
  filter(t < 13)


ab_male_mat_plot = ggplot(test,aes(x=t,y=ab_male_mat,fill=phase))+
  geom_bar(stat="identity", position=position_dodge())
ab_male_imm_plot = ggplot(test,aes(x=t,y=ab_male_imm,fill=phase))+
  geom_bar(stat="identity", position=position_dodge())
plot_grid(ab_male_mat_plot,ab_male_imm_plot,ncol = 2)


## Key variable plot (catch, cost, profit, abundance)
tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len[,,2,,],c(4),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len[,,2,,],c(4),sum,na.rm=T)

plot(tot_imm[1:11])
plot(tot_mat[1:36])

x11()
par(mfrow=c(4,1),mar=c(.1,.1,.1,.1),oma=c(4,.1,1,1))
plot(tot_imm,type='l',las=1,xaxt='n',ylim=c(0,max(tot_imm,tot_mat)))
lines(tot_mat,lty=2)
legend('topright',bty='n',lty=c(1,2),legend=c("Immature N","Mature N"))
# plot(tot_catch,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total catch"))
# plot(tot_cost,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total cost"))
# plot(tot_profit,xaxt='n',las=1)
# legend('right',bty='n',legend=c("Total profits"))
plot(tot_catch[(tot_catch>0)],xaxt='n',las=1,type='b',pch=16,ylim=c(0,max(tot_catch)))
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost[(tot_catch>0)],xaxt='n',las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit[(tot_catch>0)],las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total profits"))

# ## Plot size classes
# par(mfrow = c(4,3))
# for(x in 1:12){
#   
#   plot(imm_N_at_Len[,,2,x,600] * land_mask_na,
#        main = paste0("Size: [",sizes[x],", ",sizes[x+1],"] mm"),
#        asp=1)
#   
# }
