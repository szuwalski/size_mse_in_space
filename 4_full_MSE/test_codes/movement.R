#########################
## Simulation of movement
#########################

# Compute adjacency matrix
A_gg = dnearneigh(point_expand[,c("Var1","Var2")],d1=0.45,d2=0.65)
A_gg_mat = nb2mat(A_gg)
A_gg_mat[which(A_gg_mat > 0)] = 1
# plot(A_gg_mat[1:100,1:100])

## Diffusion
#-----------
# D * DeltaT / A
# with D: diffusion coefficient, here btwn 0.1 and 1.1 km(^2?) per day (Cf. Olmos et al. SM)
# DeltaT: time interval btwn time steps, 
# A: area of grid cells
# par(mfrow = c(3,2))
D = 2
DeltaT = (1/30) # convert day in month
A = mean(cell_area)
diffusion_coefficient = D * DeltaT / A
diffusion_gg = A_gg_mat * diffusion_coefficient
diag(diffusion_gg) = -1 * colSums(diffusion_gg)

# plot(diffusion_gg[1:100,1:100], breaks=20)
# hist(diffusion_gg)

## Taxis for adults
#------------------
taxis_coef_ad = 10^3 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations

mov_mat_temp = 4
preference_g_ad = log(init_adult + 0.0001)
preference_g_ad[which(is.na(preference_g_ad))] = min(preference_g_ad)
preference_g_ad = as.vector(t(preference_g_ad))
# # check
# test = matrix(preference_g,nrow = 40,ncol = 40,byrow = T)
# plot(test)

taxis_gg_ad = A_gg_mat *  taxis_coef_ad * DeltaT / A * exp(outer(preference_g_ad, preference_g_ad, "-") * DeltaT / sqrt(A))

diag(taxis_gg_ad) = -1 * colSums(taxis_gg_ad)

# Total
mrate_gg_ad = taxis_gg_ad # + diffusion_gg
if( any(mrate_gg_ad-diag(diag(mrate_gg_ad))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")

mfraction_gg_ad = Matrix::expm(mrate_gg_ad)
# test = matrix(mfraction_gg@x,nrow=mfraction_gg@Dim[1],ncol=mfraction_gg@Dim[2])
# plot(test[1:100,1:100])

stationary_g_ad = eigen(mfraction_gg_ad)$vectors[,1]
stationary_g_ad = stationary_g_ad / sum(stationary_g_ad)
# test = matrix(stationary_g_ad,nrow=40,ncol=40,byrow = T)
# x11();plot(test)


## Transition pattern
#--------------------
taxis_coef_trans = 10^3 # This value is set so that movement happen rapidly enough --> should be refined by some ecological considerations

trans_distrib = init_adult + init_juv
trans_distrib[which(trans_distrib > 0)] = 1
preference_g_trans = trans_distrib
preference_g_trans = as.vector(t(preference_g_trans))
# # check
# test = matrix(preference_g,nrow = 40,ncol = 40,byrow = T)
# plot(test)

taxis_gg_trans = A_gg_mat *  taxis_coef_trans * DeltaT / A * exp(outer(preference_g_trans, preference_g_trans, "-") * DeltaT / sqrt(A))

diag(taxis_gg_trans) = -1 * colSums(taxis_gg_trans)

# Total
mrate_gg_trans = taxis_gg_trans # + diffusion_gg
if( any(mrate_gg_trans-diag(diag(mrate_gg_trans))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")

mfraction_gg_trans = Matrix::expm(mrate_gg_trans)
# test = matrix(mfraction_gg@x,nrow=mfraction_gg@Dim[1],ncol=mfraction_gg@Dim[2])
# plot(test[1:100,1:100])

stationary_g_trans = eigen(mfraction_gg_trans)$vectors[,1]
stationary_g_trans = stationary_g_trans / sum(stationary_g_trans)
# # Check
test = matrix(stationary_g_trans,nrow=40,ncol=40,byrow = T)
x11();plot(test)


## Movement
##---------
sex = 2
x = 1
t = 1
mov_imm_temp<-imm_N_at_Len[,,sex,x,t]
mov_mat_temp<-mat_N_at_Len[,,sex,x,t]

x11()
par(mfrow = c(2,2))
plot(mov_imm_temp,main=paste0("t = 0"),asp = 1)
for(t in 1:1000){
  
  mov_imm_v = mfraction_gg_trans %*% as.vector(t(mov_imm_temp)) # mrate_gg %*%
  mov_mat_v = mfraction_gg_trans %*% as.vector(t(mov_mat_temp))  # mrate_gg %*%
  
  mov_imm_v = mfraction_gg_ad %*% mov_imm_v # mrate_gg %*%
  mov_mat_v = mfraction_gg_ad %*% mov_mat_v  # mrate_gg %*%
  
  mov_imm_temp_2 = matrix(mov_imm_v, nrow = length(lat), ncol = length(lon),byrow = T)
  mov_mat_temp_2 = matrix(mov_mat_v, nrow = length(lat), ncol = length(lon),byrow = T)
  
  # print(which((mov_imm_temp != mov_mat_temp_2)))
  
  mov_imm_temp = mov_imm_temp_2
  mov_mat_temp = mov_mat_temp_2

  
}

plot(mov_mat_temp,main=paste0("t = ",t),asp = 1)

