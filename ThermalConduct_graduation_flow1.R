rm(list = ls())
setwd(getwd()) # Set the working directory

#k1=k2=1
##
k1=k2=0

k1_U = 2.3 #Теплопроводность льда
#k2_U = k1_U
k2_U = 0.58 #Теплопроводность воды

density_1  = 918.7 #Плотность льда
density_2 = 999.7 #Плотность воды
#mean_density= mean(c(density1,density2))
с1 = 2000 #Теплоёмкость льда
с2 = 4195 #Теплоёмкость воды

c1_U = с1 * density_1 
c2_U = с2 * density_2 
#c2_U = c1_U

mean_K = mean(c(k1_U,k2_U))
mean_C = mean(c(c1_U,c2_U))

a=(mean_K/mean_C)^(1/2)
h = 0.00045 #шаг по x
tau = h^2/a^2

L = 0.008 #длина
N = round(L/h) #Число шагов
n = 100 #Число шагов по времени
tj = numeric(n) #tau j
nj = numeric(n) #ветор моментов времени

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

x = numeric(N+1)
fij = numeric(N+1) #сила

Ai = numeric(N) #коэффициенты
Bi = numeric(N)
Ci = numeric(N)
Fi = numeric(N-1)


alpha = numeric(N)
beta = numeric(N)

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

for(i in 1:N+1)
{
  x[i]= (i-1)*h
  fij[i]=0
}

U = matrix(data=NA,nrow=n,ncol=N+1)
colnames(U)=x
row.names(U)= c(0:(n-1))


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

options(scipen = 999) # Disable exponential notation (e.g. 1.81e+09)
#print("Численное решение")
#print(U)
#print("Точное решение")
#print(Uacc)
RESULTS_DIRECTORY <- "../Results/"

file.path = paste(RESULTS_DIRECTORY, "result_", "matrixx", ".csv", sep = "")
#write.matrix(U , file = file.path, sep = ";")
#options(digits = 4)
#Uacc = format(U, digits = 4)
#write.table(Uacc , file = file.path, sep = ";")

library(animation)
oopt = ani.options(interval = 0.3)
ani.record(reset = TRUE)

for(j in 2:n)
{
  plot(U[j,], xaxt="n", xlab = 'X values', ylab = 'U[x,t] values')
  lines(U[j,],col="red")
  if(j == 1 || j %% 25 == 0 || j == 2 )
  {
    file.path = paste(RESULTS_DIRECTORY, "result_", j, ".bmp", sep = "")
    bmp(file.path)
    plot(U[j,], xaxt="n", xlab = 'X values', ylab = 'U[x,t] values')
    lines(U[j,],col="red")
    axis(1, at = c(1:(N+1)), labels = x)
    dev.off()
  }
  #lines(Uacc[j,],col="blue")
  axis(1, at = c(1:(N+1)), labels = x)
  ani.pause()
  ani.record()
  
}


# while(TRUE)
# { 
# ani.replay()
# ani.options(oopt, loop = 100)
# }