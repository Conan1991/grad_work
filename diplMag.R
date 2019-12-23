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

lambda =  0.1
alpha = 0.33
nu = 0.25

rows = length(x)
cols = ncol(temperature)
y = x

step_x = 0.5
step_y = 0.5
step_t = 0.5

shiftU = shiftV = matrix(data = 0, nrow = rows, ncol = cols)
colnames(shiftU)= colnames(shiftV) = x
row.names(shiftU) = row.names(shiftV) = y

SHIFTS_U = list(shiftU)
SHIFTS_V = list(shiftV)


#fill 0 layer



for (t in 2:15) {
  for(i in 2:(rows-1))
  {
    for(j in 2:(cols-1))
    {
      shiftU[i,j] = step_t*((lambda+2*nu)*(SHIFTS_U[[t-1]][i+1,j]-2*SHIFTS_U[[t-1]][i,j]+SHIFTS_U[[t-1]][i-1,j])/step_x^2+(lambda+nu)*(SHIFTS_V[[t-1]][i+1,j+1]-SHIFTS_V[[t-1]][i-1,j+1]-SHIFTS_V[[t-1]][i+1,j-1]+SHIFTS_V[[t-1]][i-1,j-1])/(4*step_x*step_y)+nu*(SHIFTS_U[[t-1]][i,j-1]-2*SHIFTS_U[[t-1]][i,j]+SHIFTS_U[[t-1]][i,j+1])/step_y^2-3*lambda*alpha*(temperature[i+1,j]-temperature[i-1,j])/(2*step_x))+SHIFTS_U[[t-1]][i,j]
      shiftV[i,j] = step_t*(nu*(SHIFTS_V[[t-1]][i+1,j]-2*SHIFTS_V[[t-1]][i,j]+SHIFTS_V[[t-1]][i-1,j])/step_x^2+(lambda+nu)*(SHIFTS_U[[t-1]][i+1,j+1]-SHIFTS_U[[t-1]][i-1,j+1]-SHIFTS_U[[t-1]][i+1,j-1]+SHIFTS_U[[t-1]][i-1,j-1])/(4*step_x*step_y)+(lambda+2*nu)*(SHIFTS_U[[t-1]][i,j-1]-2*SHIFTS_U[[t-1]][i,j]+SHIFTS_U[[t-1]][i,j+1])/step_y^2-3*lambda*alpha*(temperature[i+1,j]-temperature[i-1,j])/(2*step_y))+SHIFTS_V[[t-1]][i,j]
    }
  }
  SHIFTS_U[[t]] = shiftU
  SHIFTS_V[[t]] = shiftV
}

#fit first and last layers
for (t in 2:4)
  for(i in 2:(rows-1))
  {
    SHIFTS_U[[t]][i,1]=step_y/2/step_x*(SHIFTS_V[[t]][i+1,2]-SHIFTS_V[[t]][i-1,2]) + SHIFTS_U[[t]][i,2]
    SHIFTS_V[[t]][i,1] = lambda*step_y/2/step_x*(SHIFTS_U[[t]][i+1,2]-SHIFTS_U[[t]][i-1,2])/(lambda+2*nu)-3*lambda*alpha*step_y*temperature[i,1]/(lambda+2*nu)+SHIFTS_V[[t]][i,2]
    SHIFTS_U[[t]][i,cols]= SHIFTS_U[[t]][i, cols-1] - step_y/2/step_x*(SHIFTS_V[[t]][i+1,cols]-SHIFTS_V[[t]][i-1,cols])
    SHIFTS_V[[t]][i,cols]= SHIFTS_V[[t]][i, cols-1] - lambda*step_y/2/step_x*(SHIFTS_U[[t]][i+1,cols]-SHIFTS_U[[t]][i-1,cols])/(lambda+2*nu)+3*lambda*alpha*step_y*temperature[i,cols]/(lambda+2*nu)
#print(SHIFTS_U[[t]][i,1])
  }


plot(SHIFTS_U[[2]][2,], xaxt="n")
axis(1, at = c(1:(cols)), labels = x)
