source('temperature.R')

temperature = ThermalConduct()

rows = length(x)
cols = length(x)
y = x

shiftU = shiftV = matrix(data = 0, nrow = rows, ncol = cols)
colnames(shiftU)= colnames(shiftV) = x
row.names(shiftU) = row.names(shiftV) = y

ksi11= ksi12= ksi21=ksi22=matrix(data = 0, nrow = rows, ncol = cols)
sigma11=sigma21=sigma12=sigma22=sigma33=matrix(data = 0, nrow = rows, ncol = cols)


SHIFTS_U = list(shiftU)
SHIFTS_V = list(shiftV)

checkBorders = function(iter)
{
  x_truncated = numeric()
  #y_truncated = numeric()
  for(i in 1:length(x))
  if(temperature[iter, i] < 0)
    x_truncated[i] = temperature[iter, i]

  return(x_truncated)
}

calculate_Shifts = function(iteration = 0)
{
  x_tr = checkBorders(iteration)
  for (t in 2:number_of_iterations) {
    for (i in 2:(length(x_tr) - 1))
    {
      for (j in 2:(cols - 1))
      {
        shiftU[i, j] = step_t * ((lambda + 2 * nu) * (SHIFTS_U[[t - 1]][i + 1, j] - 2 * SHIFTS_U[[t - 1]][i, j] + SHIFTS_U[[t - 1]][i - 1, j]) / step_x ^ 2 + (lambda + nu) * (SHIFTS_V[[t - 1]][i + 1, j + 1] - SHIFTS_V[[t - 1]][i - 1, j + 1] - SHIFTS_V[[t - 1]][i + 1, j - 1] + SHIFTS_V[[t - 1]][i - 1, j - 1]) / (4 * step_x * step_y) + nu * (SHIFTS_U[[t - 1]][i, j - 1] - 2 * SHIFTS_U[[t - 1]][i, j] + SHIFTS_U[[t - 1]][i, j + 1]) / step_y ^ 2 - alpha * (temperature[iteration, i + 1] - temperature[iteration, i - 1]) / (2 * step_x) * (3 * lambda + 2 * nu)) + SHIFTS_U[[t - 1]][i, j]
        shiftV[i, j] = step_t * (nu * (SHIFTS_V[[t - 1]][i + 1, j] - 2 * SHIFTS_V[[t - 1]][i, j] + SHIFTS_V[[t - 1]][i - 1, j]) / step_x ^ 2 + (lambda + nu) * (SHIFTS_U[[t - 1]][i + 1, j + 1] - SHIFTS_U[[t - 1]][i - 1, j + 1] - SHIFTS_U[[t - 1]][i + 1, j - 1] + SHIFTS_U[[t - 1]][i - 1, j - 1]) / (4 * step_x * step_y) + (lambda + 2 * nu) * (SHIFTS_V[[t - 1]][i, j - 1] - 2 * SHIFTS_V[[t - 1]][i, j] + SHIFTS_V[[t - 1]][i, j + 1]) / step_y ^ 2) + SHIFTS_V[[t - 1]][i, j]
       }
    }
    SHIFTS_U[[t]] = shiftU
    SHIFTS_V[[t]] = shiftV
    
    ##fit first and last layers
    for (k in 2:(cols - 1))
    {
      SHIFTS_U[[t]][1, k] = lambda * step_x / (2 * step_y) * (SHIFTS_V[[t - 1]][1, k + 1] - SHIFTS_V[[t - 1]][1, k - 1]) / (lambda + 2 * nu) - alpha * step_x * temperature[iteration, 1] / (lambda + 2 * nu) * (3 * lambda + 2 * nu) + SHIFTS_U[[t]][2, k]
      SHIFTS_V[[t]][1, k] = step_x / (2 * step_y) * (SHIFTS_U[[t - 1]][1, k + 1] - SHIFTS_U[[t - 1]][1, k - 1]) + SHIFTS_V[[t]][2, k]
      SHIFTS_U[[t]][length(x_tr), k] = SHIFTS_U[[t]][length(x_tr) - 1, k] - lambda * step_x /
(2 * step_y) * (SHIFTS_V[[t - 1]][length(x_tr), k + 1] - SHIFTS_V[[t - 1]][length(x_tr), k - 1]) / (lambda + 2 * nu) + alpha * step_x * temperature[iteration, length(x_tr)] / (lambda + 2 * nu) * (3 * lambda + 2 * nu)
      SHIFTS_V[[t]][length(x_tr) , k] = SHIFTS_V[[t]][length(x_tr) - 1, k] - step_x / (2 * step_y) * (SHIFTS_U[[t - 1]][length(x_tr), k + 1] - SHIFTS_U[[t - 1]][length(x_tr), k - 1])
    }
}
  
  shiftU = SHIFTS_U[[number_of_iterations]]
  shiftV = SHIFTS_V[[number_of_iterations]]
  
  for (t in 2:number_of_iterations) {
    for(i in 2:(rows-1))
    {
      for(j in 2:(cols-1))
      {
        ksi11[i,j] = (SHIFTS_U[[t]][i+1, j] - SHIFTS_U[[t]][i-1,j]) / (2 * step_x)
        ksi12[i,j] = ksi21[i,j] = ((SHIFTS_U[[t]][i,j+1] - SHIFTS_U[[t]][i,j-1]) / (2 * step_y) + (SHIFTS_V[[t]][i+1,j] - SHIFTS_V[[t]][i-1,j]) / (2*step_x)) / 2
        ksi22[i,j] = (SHIFTS_V[[t]][i,j+1] - SHIFTS_V[[t]][i,j-1]) / (2 * step_y)
        sigma11[i,j] = lambda * (ksi11[i,j] + ksi22[i,j] - 3 * alpha * temperature[iteration,i]) + 2 * nu * (ksi11[i,j] - alpha * temperature[iteration,i])
        sigma12[i,j] = sigma21[i,j] = 2 * nu * ksi12[i,j]
        sigma22[i,j] = lambda * (ksi11[i,j] + ksi22[i,j] - 3 * alpha * temperature[iteration,i]) + 2 * nu * (ksi22[i,j] - alpha * temperature[iteration,i])
        sigma33[i,j]= lambda * (ksi11[i,j] + ksi22[i,j] - 3 * alpha * temperature[iteration,i]) - 2 * nu * alpha * temperature[iteration,i]
      }
    }
  }
  
  ksi = list(ksi11,ksi12,ksi21,ksi22)
  sigma = list(sigma11,sigma12,sigma21,sigma22,sigma33)
  
  delta_T = abs(temperature[iteration,1]-temperature[iteration,length(x_tr)])
  #for(i in 1:rows)
  #  U[i] = y[rows]*alpha*delta_T/(8*x[rows])
  Test_U = y[rows]*alpha*delta_T/(8*x[length(x_tr)])
    
  result = list(shift.u = shiftU, shift.v = shiftV, ksi = ksi, sigma = sigma, test.u = Test_U)
  return(result)
}