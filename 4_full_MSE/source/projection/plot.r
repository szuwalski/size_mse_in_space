## Plot projection
#-----------------

tot_catch<-apply(catch_by_fisher,c(5),sum,na.rm=T)
tot_cost<-apply(cost_by_fisher,c(3),sum,na.rm=T)
tot_profit<-apply(profit_by_fisher,c(3),sum,na.rm=T)

tot_imm<-apply(imm_N_at_Len,c(5),sum,na.rm=T)
tot_mat<-apply(mat_N_at_Len,c(5),sum,na.rm=T)

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
plot(tot_catch[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16,ylim=c(0,max(tot_catch)))
legend('right',bty='n',legend=c("Total catch"))
plot(tot_cost[(tot_profit>0)],xaxt='n',las=1,type='b',pch=16)
legend('right',bty='n',legend=c("Total cost"))
plot(tot_profit[(tot_profit>0)],las=1,type='b',pch=16)
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
