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
#U_res=matrix(data = NA, nrow = rows, ncol = cols)
#colnames(U_res) = x
#row.names(U_res) = y
#currentTemperature = numeric()
Test_U = numeric(n)

RESULTS_DIRECTORY <- paste(Sys.getenv("HOME"), "\\GitHub\\grad_work\\GradResults\\", sep="")


for(it in 3:n)
{
  res = calculate_Shifts(it)
  if(is.null(res))
  {
    break
  }
  SHIFTSU_RESULT[[it]] = res$shift.u
  SHIFTSV_RESULT[[it]] = res$shift.v
  ksi_result[[it]] = res$ksi
  sigma_result[[it]] = res$sigma
  Test_U[it] = res$test.u
  #file.path = paste(RESULTS_DIRECTORY, "result_", it, ".csv", sep = "")
  #Uacc = formatC(SHIFTSU_RESULT[[it]], digits = 6, format = "f")
  #write.table(Uacc, file = file.path, sep = ";")
}

#plot(SHIFTSU_RESULT[[3]][3,], xaxt="n")
#axis(1, at = c(1:(cols)), labels = x)
#lines(SHIFTSU_RESULT[[3]][3,],col="red")
last = length(SHIFTSU_RESULT)
#mizes calculation
sig = sigma_result[[last]]
mizes_s = matrix(data = 0, nrow = rows, ncol = cols) #критерий прочности
#mizes_res = list(mizes_s)
for( k in 1:cols)
{
  for(m in 1:rows)
  {
      mizes_s[k,m] = sqrt( (sig[[1]][k,m] - sig[[4]][k,m])^2 + (sig[[4]][k,m] - sig[[5]][k,m])^2 + (sig[[5]][k,m] - sig[[1]][k,m])^2 + 6 * sig[[2]][k,m]^2)
  }
}


#shifts calculation
shiftsV=SHIFTSU_RESULT[[last]]
shiftsH=SHIFTSV_RESULT[[last]]

#old_matrix = new_matrix = list()
list_p_x = list_p_y = list()
list_new_x = list_new_y = list()
#for(t in 1:n)
#{
  for(i in 1:nrow(shiftsV))
  {
    vecp_x = vecp_new_x = c()
    vecp_y = vecp_new_y = c()
    
    
    for(j in 1:ncol(shiftsV))
    {
      vecp_x[j] <- x[i]
      vecp_y[j] <- y[j]
      vecp_new_x[j] <- vecp_x[j] + shiftsV[i,j]
      vecp_new_y[j] <- vecp_y[j] + shiftsH[i,j]
      #list_p[[j]] <- c(x[i], y[j])
      #v_point <- c(shiftsV[i,j], shiftsH[i,j])
      #list_newp[[j]] <- list_p[[j]] + v_point
    }
    list_p_x[[i]] = vecp_x
    list_p_y[[i]] = vecp_y
    
    list_new_x[[i]] = vecp_new_x
    list_new_y[[i]] = vecp_new_y
    #old_matrix[[i]] = list_p
    #new_matrix[[i]] = list_newp
  }
#}

#axis(1, at = c(1:(cols)), labels = x)
#lines(new_matrix[[1]],col="red")

library(animation)
oopt = ani.options(interval = 0.3)

ani.record(reset = TRUE)
for(j in 1:nrow(shiftsV))
{
  file.path = paste(RESULTS_DIRECTORY,"\\graph\\" , "result_", j, ".bmp", sep = "")
  #bmp(file.path)
  
  plot(list_p_y[[j]], list_p_x[[j]], xaxt="n",  xlab="value of y" , ylab="value of x")
  points(list_new_y[[j]], list_new_x[[j]], type = "p")
  axis(1, at = list_p_y[[j]])
  #points(old_point_x, point_x,type = "p")
  lines(list_p_y[[j]], list_p_x[[j]], col="red")
  lines(list_new_y[[j]], list_new_x[[j]], col="blue")
  
  #dev.off()

  ani.pause()
  ani.record()
  
}

while(TRUE)
{ 
  ani.replay()
  ani.options(oopt, loop = 100)
}
#write(Test_U, file = "tempU.txt")
#write.table(Test_U, "Test_U_3.csv", col.names = "Test_U values", row.names = FALSE, sep = ";")
#write.table(mizes_s, "mizes.csv", col.names = x, row.names = y, sep = ";")
