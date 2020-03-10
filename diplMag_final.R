rm(list = ls())
setwd(getwd()) # Set the working directory

options(scipen = 999)
source('shifts.R')



temperature = ThermalConduct()
#fill 0 layer
SHIFTSU_RESULT = list(shiftU)
SHIFTSV_RESULT = list(shiftU)


checkResidual = function(U, it)
{
  MIN_RESIDUAL = 10000
  max_residuals=numeric()
  for (t in it:number_of_iterations) {
    residuals = numeric()
    
    for(i in 1:(rows-1))
    {
      for(j in 1:(cols-1))
      {
        if(abs(U[[t]][i,j]-U[[t-1]][i,j]) == 0)
          next()
        #print(abs(U[[t]][i,j]-U[[t-1]][i,j]))
        residuals <- c(residuals, abs(U[[t]][i,j]-U[[t-1]][i,j]) )
      }
    }
    
    max_residuals <-c(max_residuals,max(residuals))
    
  }
  print(max_residuals)
  for(i in 1:length(max_residuals))
  {
    maximum = max_residuals[i]
    if(i!= 1)
    {
      if(max_residuals[i] > max_residuals[i-1])
        return(FALSE)
      if(maximum < MIN_RESIDUAL)
        assign("MIN_RESIDUAL", maximum, envir = .GlobalEnv)
    }
    
  }
  
  if(tail(max_residuals, n=1) < epsilon)
  {
    #print(tail(max_residuals, n=1))
    print("Found minimum residual")
    return(TRUE)
  }
  
  return(FALSE)
  
}



currentTemperature = numeric()

for(it in 2:nrow(temperature))
{
  res = calculate_Shifts(it)
  SHIFTSU_RESULT[[it]] = res$shift.u
  SHIFTSV_RESULT[[it]] = res$shift.v
  
}



RESULTS_DIRECTORY <- paste(Sys.getenv("HOME"), "\\GitHub\\grad_work\\GradResults\\", sep="")

#file.path = paste(RESULTS_DIRECTORY, "result_", "matrix", ".csv", sep = "")
#Uacc = formatC(SHIFTS_U[[number_of_iterations]], digits = 6, format = "f")
#write.table(Uacc, file = file.path, sep = ";")
