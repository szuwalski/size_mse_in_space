## Recruitment
#-------------
#==this is a dumb temporary fix
#==this implants the original recruitment with some error
#==ultimately we need an algorithm to determine location and intensity of recruitment
#==teleconnections postdoc will hopefully help with this

if(recruit_model == "cody_model"){
  
  meanlog_recruit_male_2 = meanlog_recruit_male
  meanlog_recruit_female_2 = meanlog_recruit_female
  
}else if(recruit_model == "max_model"){
  
  source("4_full_MSE/LHP/recruit_t.R")
  
}

if(R_stochastic == "yes"){
  
  aggreg_rec_1 = rlnorm(n = 1, mean = meanlog_recruit_female_2, sdlog = sdlog_female)
  aggreg_rec_2 = rlnorm(n = 1, mean = meanlog_recruit_male_2, sdlog = sdlog_male)
  
}else if(R_stochastic == "no"){
  
  aggreg_rec_1 = exp(meanlog_recruit_female_2)
  aggreg_rec_2 = exp(meanlog_recruit_male_2)
  
}else if(R_stochastic == "pre-simulated"){
  
  aggreg_rec_1 = vec_rec_1[Years_climsc[t]-Start_Y+1]
  aggreg_rec_2 = vec_rec_2[Years_climsc[t]-Start_Y+1]
  
}

tmp_rec_1 <- init_juv * aggreg_rec_1
tmp_rec_2 <- init_juv * aggreg_rec_2

for(r in 1:rec_sizes)
{
  # print(sum(temp_imm_N[,,2,]))
  temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + tmp_rec_1 * prop_rec[r]
  temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + tmp_rec_2*prop_rec[r]
  # print(sum(temp_imm_N[,,2,]))
  
}

