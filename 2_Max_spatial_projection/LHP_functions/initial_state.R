initial_state <- function(fake_dist_data, n_s, n_n, n_p, n_t,plot=FALSE) {
  
  #==need to track numbers at size by sex by maturity by location
  imm_N_at_Len <- array(dim = c(n_lat, n_lon, n_n, n_p, n_t))
  dim(imm_N_at_Len)
  mat_N_at_Len <- array(dim = c(n_lat, n_lon, n_n, n_p, n_t))
  
  load(file = "2_Max_spatial_projection/Input_data/initial_states/smooth_mat_N_at_Len_2017.RData") # smooth_mat
  load(file = "2_Max_spatial_projection/Input_data/initial_states/smooth_imm_N_at_Len_2017.RData") # smooth_imm
  
  
  #==FIX THIS SO THAT THESE DISTRIBUTIONS ARE MADE TO TRUE DISTRIBUTIONS AND MULTIPLIED BY NUMBERS AT LENGTH FROM ASSESSMENT
  imm_N_at_Len[, , , , 1] <- exp(smooth_imm)
  mat_N_at_Len[, , , , 1] <- exp(smooth_mat)
  
  imm_N_at_Len[imm_N_at_Len == "NaN"] <- 0
  mat_N_at_Len[mat_N_at_Len == "NaN"] <- 0
  
  imm_N_at_Len[imm_N_at_Len < 0] <- 0
  mat_N_at_Len[mat_N_at_Len < 0] <- 0
  
  
  #==This makes up a random distribution for the population
  #==only used for testing
  fake_dist_data <- 0
  if (fake_dist_data == 1)
  {
    #==set dummy initial distribution--take this from the survey for real
    #==this is only for getting mechanics down
    imm_N_at_Len[, , 1, 1, 1] <-
      rnorm(dim(imm_N_at_Len[, , 1, 1, 1])[1] * dim(imm_N_at_Len[, , 1, 1, 1])[2], 1000, 100)
    imm_N_at_Len[, , 2, 1, 1] <-
      rnorm(dim(imm_N_at_Len[, , 1, 1, 1])[1] * dim(imm_N_at_Len[, , 1, 1, 1])[2], 1000, 100)
    mat_N_at_Len[, , 1, 1, 1] <-
      rnorm(dim(imm_N_at_Len[, , 1, 1, 1])[1] * dim(imm_N_at_Len[, , 1, 1, 1])[2], 1000, 100)
    mat_N_at_Len[, , 2, 1, 1] <-
      rnorm(dim(imm_N_at_Len[, , 1, 1, 1])[1] * dim(imm_N_at_Len[, , 1, 1, 1])[2], 1000, 100)
    
    for (x in 2:length(sizes))
    {
      imm_N_at_Len[, , 1, x, 1] <- imm_N_at_Len[, , 1, x - 1, 1] * exp(-M)
      imm_N_at_Len[, , 2, x, 1] <- imm_N_at_Len[, , 2, x - 1, 1] * exp(-M)
      mat_N_at_Len[, , 1, x, 1] <- mat_N_at_Len[, , 1, x - 1, 1] * exp(-M)
      mat_N_at_Len[, , 2, x, 1] <- mat_N_at_Len[, , 2, x - 1, 1] * exp(-M)
    }
    
  }
  if(plot == TRUE) plot(mat_N_at_Len)
  #==ugh. this is dumb, but how to automate?
  g <- function(m)
    t(m)[, nrow(m):1]
  
  land_mask[25:40, 20:40]
  land_mask[31, 32] <- 0
  land_mask[32, 29] <- 0
  filled.contour(g(land_mask))
  filled.contour(x = lon, y = rev(lat), g(imm_N_at_Len[, , 1, 1, 1]))
  #write.csv(land_mask,'landmask.csv')
  
  #==ensures no critters on land
  for (x in 1:2)
    for (y in 1:length(sizes))
    {
      imm_N_at_Len[, , x, y, 1] <- imm_N_at_Len[, , x, y, 1] * land_mask
      mat_N_at_Len[, , x, y, 1] <- imm_N_at_Len[, , x, y, 1] * land_mask
    }
  dim(mat_N_at_Len)
  #==SHOULD THERE ALSO BE A 'DEPTH MASK'???
  #==MAYBE A 'FISHERY MASK'??  Depth should limit the fishery.
  filled.contour(x = lon, y = rev(lat), g(imm_N_at_Len[, , 1, 1, 1]))
  
  
  #==create a file with bottom temperature for all time periods
  avg_bot_tmp <- rep(c(1, 2, 3, 3, 3, 2, 1, 0, 0, 0, 0, 1), year_n)
  for (x in 1:length(proj_period))
  {
    temp <- land_mask * rnorm(length(land_mask), avg_bot_tmp[x], .5)
    #write.csv(temp,paste("temp_data/bot_temp_",x,".csv",sep=""),row.names=FALSE)
  }
  
  return(list("imm_N_at_Len"=imm_N_at_Len,"mat_N_at_Len"=mat_N_at_Len))
  
}
