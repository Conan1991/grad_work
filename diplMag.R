rm(list = ls())
setwd(getwd()) # Set the working directory

options(scipen = 999)
source('constants.R')


ThermalConduct = function()
{
  k1=k2=0
  
  
  tj = numeric(n) #tau j
  nj = numeric(n) #ветор моментов времени
  
  
  fij = numeric(N+1) #сила
  
  Ai = numeric(N) #коэффициенты
  Bi = numeric(N)
  Ci = numeric(N)
  Fi = numeric(N-1)
  
  alpha = numeric(N)
  beta = numeric(N)
  
  for(i in 1:N+1)
  {
    fij[i]=0
  }
  
  U = matrix(data=NA,nrow=n,ncol=N+1)
  colnames(U)=x
  row.names(U)= c(0:(n-1))
  
  
  kU = function(U_i)
  {
    if( U_i < 0)
      return(k1_U)
    else return(k2_U)
  }
  
  CU = function(U_i)
  {
    if( U_i < 0)
      return(c1_U)
    else return(c2_U)
  }
  
  
  A = function(t)
  {
    #0
    #1/(1+t)
    -10
  }
  B = function(t)
  {
    #0
    #1/(1+t)
    10
  }
  
  for(j in 1:n)
  {
    tj[j]=(j-1)*tau
    nj[j]=j-1
  }
  
  Ux0=function(x)
  {
    #A(tj[1]) + (B(tj[1])-A(tj[1]))*(x/L)^2;
    #1
    10
  }
  
  
  for(i in 0:N+1) #Считаем 0-й слой
  {
    U[1,i] = Ux0(x[i])
  }
  
  CoeffF = function(j) #Поправка на коэфф
  {
    alpha[1] = k1
    beta[1] = A(tj[j])
    Ai[1] = - kU(U[j,1]) / h^2
    Bi[1] = - kU(U[j,1]) / h^2
    Ci[1] = CU(U[j,1]) / tau + 2 * kU(U[j,1]) / h^2
    
    for(i in 2:(N-1)) # Считаем очередные Ai, Bi , Ci
    {
      Ai[i] = - kU(U[j,i-1]) / h^2
      Bi[i] = - kU(U[j,i+1]) / h^2
      Ci[i] = CU(U[j,i]) / tau + (kU(U[j,i+1]) + kU(U[j,i-1]))  / h^2
    }
    
    for(i in 1:(N-1)) #Считаем Fi
      Fi[i] = CU(U[j,i]) * U[j-1,i+1] / tau
    
    for(i in 2:N) #Считаем альфа и бета коэффициенты
    {
      alpha[i]=-Bi[i-1]/(Ci[i-1]+Ai[i-1]*alpha[i-1])
      beta[i]= (Fi[i-1]-Ai[i-1]*beta[i-1])/(Ci[i-1]+Ai[i-1]*alpha[i-1])
    }
    
    U[j,N+1]=(B(tj[j])+k2*beta[N])/(1-k2*alpha[N])
    for(i in N:1) #Считаем Uj
      U[j,i]= alpha[i] * U[j,i+1] + beta[i]
    
    return(U[j,])
  }
  
  IterF = function(j) #Считаем остальные слои
  {
    
    alpha[1] = k1
    beta[1] = A(tj[j])
    
    for(i in 1:(N-1)) # Считаем очередные Ai, Bi , Ci
    {
      Ai[i] = - kU(U[j-1,i]) / h^2
      Bi[i] = - kU(U[j-1,i]) / h^2
      Ci[i] = CU(U[j-1,i]) / tau + 2 * kU(U[j-1,i]) / h^2
    }
    
    for(i in 1:(N-1)) #Считаем Fi
      Fi[i] = CU(U[j-1,i])*U[j-1,i+1]/tau
    
    for(i in 2:N) #Считаем альфа и бета коэффициенты
    {
      alpha[i]=-Bi[i-1]/(Ci[i-1]+Ai[i-1]*alpha[i-1])
      beta[i]= (Fi[i-1]-Ai[i-1]*beta[i-1])/(Ci[i-1]+Ai[i-1]*alpha[i-1])
    }
    
    U[j,N+1]=(B(tj[j])+k2*beta[N])/(1-k2*alpha[N])
    
    for(i in N:1) #Считаем Uj
      U[j,i]= alpha[i] * U[j,i+1] + beta[i]
    #U[j,1]=k1*U[j,2]-h*A(tj[j])
    
    return(U[j,])
  }
  
  for(j in 2:n)
  {
    U[j,] = IterF(j)
    U[j,] = CoeffF(j)
  }
  
  return(U)
}


temperature = ThermalConduct()
#fill 0 layer

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


for (t in 2:number_of_iterations) {
  for(i in 2:(rows-1))
  {
    for(j in 2:(cols-1))
    {
      shiftU[i,j] = step_t*((lambda+2*nu)*(SHIFTS_U[[t-1]][i+1,j]-2*SHIFTS_U[[t-1]][i,j]+SHIFTS_U[[t-1]][i-1,j])/step_x^2+(lambda+nu)*(SHIFTS_V[[t-1]][i+1,j+1]-SHIFTS_V[[t-1]][i-1,j+1]-SHIFTS_V[[t-1]][i+1,j-1]+SHIFTS_V[[t-1]][i-1,j-1])/(4*step_x*step_y)+nu*(SHIFTS_U[[t-1]][i,j-1]-2*SHIFTS_U[[t-1]][i,j]+SHIFTS_U[[t-1]][i,j+1])/step_y^2-3*lambda*alpha*(temperature[2,i+1]-temperature[2,i-1])/(2*step_x))+SHIFTS_U[[t-1]][i,j]
      shiftV[i,j] = step_t*(nu*(SHIFTS_V[[t-1]][i+1,j]-2*SHIFTS_V[[t-1]][i,j]+SHIFTS_V[[t-1]][i-1,j])/step_x^2+(lambda+nu)*(SHIFTS_U[[t-1]][i+1,j+1]-SHIFTS_U[[t-1]][i-1,j+1]-SHIFTS_U[[t-1]][i+1,j-1]+SHIFTS_U[[t-1]][i-1,j-1])/(4*step_x*step_y)+(lambda+2*nu)*(SHIFTS_V[[t-1]][i,j-1]-2*SHIFTS_V[[t-1]][i,j]+SHIFTS_V[[t-1]][i,j+1])/step_y^2)+SHIFTS_V[[t-1]][i,j]
    }
    
  }
  SHIFTS_U[[t]] = shiftU
  SHIFTS_V[[t]] = shiftV
  ##fit first layers
 for(k in 2:(cols-1))
 {
   SHIFTS_U[[t]][1,k]=lambda*step_x/(2*step_y)*(SHIFTS_V[[t-1]][1,k+1]-SHIFTS_V[[t-1]][1,k-1])/(lambda+2*nu) - 3*lambda*alpha*step_x*temperature[2,1]/(lambda+2*nu) + SHIFTS_U[[t]][2,k]
   SHIFTS_V[[t]][1,k] = step_x/(2*step_y)*(SHIFTS_U[[t-1]][1,k+1]-SHIFTS_U[[t-1]][1,k-1])+SHIFTS_V[[t]][2,k]
   SHIFTS_U[[t]][rows, k]= SHIFTS_U[[t]][rows-1, k] - lambda*step_x/(2*step_y)*(SHIFTS_V[[t-1]][rows, k+1]-SHIFTS_V[[t-1]][rows,k-1])/(lambda+2*nu)+3*lambda*alpha*step_x*temperature[1,rows]/(lambda+2*nu)
   SHIFTS_V[[t]][rows , k]= SHIFTS_V[[t]][rows-1, k] - step_x/(2*step_y)*(SHIFTS_U[[t-1]][rows,k+1]-SHIFTS_U[[t-1]][rows,k-1])
 }
}

RESULTS_DIRECTORY <- paste(Sys.getenv("HOME"), "\\GitHub\\grad_work\\GradResults\\", sep="");

#file.path = paste(RESULTS_DIRECTORY, "result_", "matrix", ".csv", sep = "")
#Uacc = formatC(SHIFTS_U[[number_of_iterations]], digits = 6, format = "f")
#write.table(Uacc, file = file.path, sep = ";")

#fit first and last layers
#for (t in 2:number_of_iterations)
  

#print(checkResidual(SHIFTS_U, number_of_iterations - 500))
#print(checkResidual(SHIFTS_U, number_of_iterations - 500))

ksi11= ksi12= ksi21=ksi22=matrix(data = 0, nrow = rows, ncol = cols)
sigma11=sigma21=sigma12=sigma22=matrix(data = 0, nrow = rows, ncol = cols)

for (t in 2:number_of_iterations) {
  for(i in 2:(rows-1))
  {
    for(j in 2:(cols-1))
    {
      ksi11[i,j] = (SHIFTS_U[[t]][i+1, j] - SHIFTS_U[[t]][i-1,j])/(2*step_x)
      ksi12[i,j]=ksi21[i,j] = ((SHIFTS_U[[t]][i,j+1] - SHIFTS_U[[t]][i,j-1])/(2*step_y) + (SHIFTS_V[[t]][i+1,j] - SHIFTS_V[[t]][i-1,j])/(2*step_x))/2
      ksi22[i,j] = (SHIFTS_V[[t]][i,j+1] - SHIFTS_V[[t]][i,j-1])/(2*step_y)
      sigma11[i,j] = lambda*(ksi11[i,j]+ksi22[i,j]-3*alpha*temperature[2,i]) + 2*nu*ksi11[i,j]
      sigma12[i,j] = sigma21[i,j] = 2*nu*ksi12[i,j]
      sigma22[i,j] = lambda*(ksi11[i,j]+ksi22[i,j]-3*alpha*temperature[2,i]) + 2*nu*ksi22[i,j]
    }
  }
}

ksi = list(ksi11,ksi12,ksi21,ksi22)
sigma = list(sigma11,sigma12,sigma21,sigma22)



plot(SHIFTS_U[[number_of_iterations]][10,], xaxt="n")
axis(1, at = c(1:(cols)), labels = x)
lines(SHIFTS_U[[number_of_iterations]][10,],col="red")
