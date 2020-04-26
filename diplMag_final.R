rm(list = ls())
setwd(getwd()) # Set the working directory

options(scipen = 999)
source('shifts.R')

#temperature = ThermalConduct()
#fill 0 layer
SHIFTSU_RESULT = list()
SHIFTSV_RESULT = list()
ksi_result   = list()
sigma_result = list()
U_res=matrix(data = NA, nrow = rows, ncol = cols)

#currentTemperature = numeric()
Test_U = numeric(n)

RESULTS_DIRECTORY <- paste(Sys.getenv("HOME"), "\\GitHub\\grad_work\\GradResults\\", sep="")


for(it in 3:n)
{
  res = calculate_Shifts(it)
  SHIFTSU_RESULT[[it]] = res$shift.u
  SHIFTSV_RESULT[[it]] = res$shift.v
  ksi_result[[it]] = res$ksi
  sigma_result[[it]] = res$sigma
  Test_U[it] = res$test.u
  #file.path = paste(RESULTS_DIRECTORY, "result_", it, ".csv", sep = "")
  #Uacc = formatC(SHIFTSU_RESULT[[it]], digits = 6, format = "f")
  #write.table(Uacc, file = file.path, sep = ";")
}

#print(SHIFTSU_RESULT[[5]][2,2])
 
#mizes calculation
sig = sigma_result[[4]]
mizes_s = matrix(data = 0, nrow = rows, ncol = cols)
mizes_res = list(mizes_s)
for( k in 1:cols)
{
  for(m in 1:rows)
  {
      mizes_s[k,m] = sqrt((sig[[1]][k,m] - sig[[4]][k,m])^2 + (sig[[4]][k,m] - sig[[5]][k,m])^2 + (sig[[5]][k,m] - sig[[1]][k,m])^2 + 6 * sig[[2]][k,m]^2)
  }
}



#write(Test_U, file = "tempU.txt")
write.table(Test_U, "Test_U_3.csv", col.names = "Test_U values", row.names = FALSE, sep = ";")

