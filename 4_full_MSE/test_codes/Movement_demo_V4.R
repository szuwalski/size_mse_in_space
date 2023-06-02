
#library(Matrix)
n_g = 50

# Domain characteristics
lat_g = seq(55, 65, length=n_g)
Temp_g = seq(15, 5, length=n_g)

# Parameters
# if testing diffusion via comparison with variance of displacement
  diffusion_coefficient = 2 ^ 2
  preference_g = -0 * (Temp_g - 10)^2
# if testing dynamics for taxis
  #diffusion_coefficient = 1 ^ 2
  #preference_g = -0.1 * (Temp_g - 10)^2

# Movement operator
A_gg = ifelse( round(abs(outer(lat_g,lat_g,"-")),2) == round(mean(diff(lat_g)),2), 1, 0 )
# Diffusion
diffusion_gg = A_gg * diffusion_coefficient   
diag(diffusion_gg) = -1 * colSums(diffusion_gg)
# Taxis
taxis_gg = A_gg * outer(preference_g, preference_g, "-")
diag(taxis_gg) = -1 * colSums(taxis_gg)
# Total
mrate_gg = diffusion_gg + taxis_gg
if( any(mrate_gg-diag(diag(mrate_gg))<0) ) stop("`mrate_gg` must be a Metzler matrix. Consider changing parameterization")
# Annualized
mfraction_gg = Matrix::expm(mrate_gg)

# plot
matplot( scale(cbind(lat_g,Temp_g,preference_g,mfraction_gg[,ceiling(0.25*n_g)],mfraction_gg[,ceiling(0.75*n_g)])), type="l", lty="solid", lwd=2 )

#
Matrix::image(mfraction_gg[1:150,1:150], xlab="From", ylab="To",breaks = )
test = matrix(mfraction_gg@x,nrow = mfraction_gg@Dim[1],ncol = mfraction_gg@Dim[2])
plot(test)

# Stationary distribution
stationary_g = eigen(mfraction_gg)$vectors[,1]
stationary_g = stationary_g / sum(stationary_g)

#
matplot( x=lat_g, y=cbind(preference_g-min(preference_g),stationary_g), type="l", col=c("black","blue") )

# n(t+1) = Mrate_gg * n(t)

# Diffusion should equal variance of displacement
  # Using central cell to minimize boundary issues
# sum( mfraction_gg[,6] * 1:n_g )
mean1 = weighted.mean(1:n_g, w=mfraction_gg[,ceiling(n_g/2)])
# sum( mfraction_gg[,6] * (1:n_g-mean1)^2 )
(var1 = weighted.mean( (1:n_g-mean1)^2, w=mfraction_gg[,ceiling(n_g/2)]) )
diffusion_coefficient * 2  # 2*D because variance of displacement is 2*D in one-dimension

# ----------------------------------------------------------------------------------------
