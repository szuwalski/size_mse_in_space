#####################
## Time series setting
#####################
temp_resolution = "year" # "month" or "year"
year_n	 <- 30 # Number of years
if(temp_resolution == "year") year_step <- 1 # Number of time steps
if(temp_resolution == "month") year_step <- 12
proj_period	<- seq(1,year_step*year_n) # Projection period
Start_Y = 2019 # Starting year
Years_climsc <- rep(Start_Y:(Start_Y+year_n), each=year_step) # Vector of years on the full time series
n_t <- length(proj_period) # Number of time steps
clim_sc = c("rcp45") # Climate scenario

## Seasonality vectors
#---------------------
#==Binary vectors related to period that determine when life events happen
#==July,Aug,Sept,Oct,Nov,Dec,Jan,Feb,Mar,Apr,May,Jun

## Yearly seasonality
# Easier to manage for GMACS ----
if(temp_resolution == "month"){
  
  survey_time_Yr   <-c(1,0,0,0,0,0,0,0,0,0,0,0)
  fish_time_Yr		 <-c(0,0,0,0,0,0,0,1,1,1,0,0)
  recruit_time_Yr	 <-c(0,0,0,0,0,0,0,0,0,1,0,0)
  move_time_Yr		 <-c(1,1,1,1,1,1,1,1,1,1,1,1)
  Fem_molt_time_Yr <- c(0,0,0,0,0,0,0,0,0,1,0,0)
  Mal_molt_time_Yr <- c(0,0,0,0,0,0,0,0,0,1,0,0)
  mate_time_Yr 	   <-c(0,0,0,0,0,0,0,1,0,1,0,0)
  SA_time_Yr       <-c(0,0,0,0,0,0,0,0,0,1,0,0)
  january_month_Yr <-c(1,0,0,0,0,0,0,0,0,0,0,0)
  
}else if(temp_resolution == "year"){
  
  survey_time_Yr   <-c(1)
  fish_time_Yr		 <-c(1)
  recruit_time_Yr	 <-c(1)
  move_time_Yr		 <-c(1)
  Fem_molt_time_Yr <- c(1)
  Mal_molt_time_Yr <- c(1)
  mate_time_Yr 	   <-c(1)
  SA_time_Yr       <-c(1)
  january_month_Yr <-c(1)
  
}

## Seasonality on the whole projection period
survey_time   <-rep(survey_time_Yr,year_n)
fish_time		  <-rep(fish_time_Yr,year_n)
recruit_time	<-rep(recruit_time_Yr,year_n)
move_time		  <-rep(move_time_Yr,year_n)
molt_time  	  <-rbind(rep(Fem_molt_time_Yr,year_n),rep(Mal_molt_time_Yr,year_n)) # females first, males second
mate_time 	  <-rep(mate_time_Yr,year_n)
SA_time       <-rep(SA_time_Yr,year_n)
january_month <-rep(january_month_Yr,year_n)

