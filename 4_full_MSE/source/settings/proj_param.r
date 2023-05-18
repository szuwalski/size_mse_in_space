## Projection setting
#-------------------
year_n	 <- 50
year_step <- 12
proj_period	<- seq(1,year_step*year_n)
Years_climsc <- rep(2022:(2022+year_n), each=year_step)
n_t <- length(proj_period)
clim_sc = c("rcp45") # climate scenario