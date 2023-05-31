## Recruitment
#-------------
#==this is a dumb temporary fix
#==this implants the original recruitment with some error
#==ultimately we need an algorithm to determine location and intensity of recruitment
#==teleconnections postdoc will hopefully help with this

print('Recruit')

if(recruit_model == "cody_model"){
  
  meanlog_recruit_male_2 = meanlog_recruit_male
  meanlog_recruit_female_2 = meanlog_recruit_female
  
}else if(recruit_model == "max_model"){
  
  source("4_full_MSE/LHP/recruit_t.R")
  
}

if(R_stochastic){
  
  aggreg_rec_1 = rlnorm(n = 1, mean = meanlog_recruit_female_2, sdlog = sdlog_female)
  aggreg_rec_2 = rlnorm(n = 1, mean = meanlog_recruit_male_2, sdlog = sdlog_male)
  
}else{
  
  aggreg_rec_1 = exp(meanlog_recruit_female_2)
  aggreg_rec_2 = exp(meanlog_recruit_male_2)
}


tmp_rec_1 <- init_juv * aggreg_rec_1
tmp_rec_2 <- init_juv * aggreg_rec_2

for(r in 1:rec_sizes)
{
  temp_imm_N[,,1,r] <- temp_imm_N[,,1,r] + imm_N_at_Len[,,1,r,1]*tmp_rec_1*prop_rec[r]
  temp_imm_N[,,2,r] <- temp_imm_N[,,2,r] + imm_N_at_Len[,,2,r,1]*tmp_rec_2*prop_rec[r]
}
