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
    10
  }
  B = function(t)
  {
    #0
    #1/(1+t)
    -10
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
    -10
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