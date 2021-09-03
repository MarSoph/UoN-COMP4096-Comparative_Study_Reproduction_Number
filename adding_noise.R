adding_noise_to_J=function(data,parms){
  'simulates data for reported cases 
  (both symptomatic and asymptomatic)
  with poisson random noise'
  u=parms[1]
  uv=parms[2]
  g=parms[3]
  uA=parms[4]
  lambda=numeric(151)
  for (i in 1:151){
    lambda[i+1]=(u*g*data$I[i]+uA*g*data$A[i]+
                   uv*g*data$Iv[i])*888000 #calculate the number of reported cases each day
  }
  noise=(rpois(length(lambda),lambda)) #simulate random poisson noise with lambda equal 
  #to the number of reported cases each day
  return(noise/888000) #returns the proportion of the population not the number of cases
}

parms=c(0.9,0.9,1/3,0.75)
J0_013=adding_noise_to_J(data=out0_013, parms) 
J0_005=adding_noise_to_J(data=out0_005, parms) 
J0_025=adding_noise_to_J(data=out0_025, parms)
#Each simulated dataset was saved into an .csv file for later use

