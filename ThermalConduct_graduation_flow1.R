rm(list = ls())

#k1=k2=1
k1=k2=0


k1_U = 2.3
#k2_U = k1_U
k2_U = 0.58
density1  = 918.7
density2 = 999.7
mean_density= mean(c(density1,density2))


c1_U = 2000 * mean_density
#c2_U = 4195 * mean_density
c2_U = c1_U
a=(k1_U/c1_U)^(1/2)
h = 0.000350 #шаг по x
#print(h)
tau = h^2/a^2



L = 0.008 #длина
N = round(L/h) #Число шагов
n = 150 #Число шагов по времени
tj = numeric(n) #tau j
nj = numeric(n) #ветор моментов времени

# k2_U= 1
# k1_U= 1
# 
# c2_U = 1
# c1_U = 1

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
Ai[1] = - kU(U[j-1,1]) / h^2
Bi[1] = - kU(U[j-1,1]) / h^2
Ci[1] = CU(U[j-1,1]) / tau + 2 * kU(U[j-1,1]) / h^2

for(i in 2:(N-1)) # Считаем очередные Ai, Bi , Ci
{
  Ai[i] = - kU(U[j-1,i-1]) / h^2
  Bi[i] = - kU(U[j-1,i+1]) / h^2
  Ci[i] = CU(U[j-1,i]) / tau + (kU(U[j-1,i+1]) + kU(U[j-1,i-1]))  / h^2
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
  #U[j,] = CoeffF(j)
}

# Uacc = matrix(data=NA,nrow=n,ncol=N+1)
# colnames(Uacc)=x
# row.names(Uacc)= c(0:(n-1))
# 
# summm=function(x,t)
# {
# summ=0
#      for(k in 1:10000)
#      {
#         #summ= summ + 1/(2*k-1)^3*sin((2*k-1)*pi*x/L) * exp(-(2*k-1)^2*pi^2*a^2*t/L^2)
#         #if(k==10) print(summ)
#        summ= summ + (-1)^(k+1)/k^2*exp(-k^2*pi^2*a^2*t/L^2)*cos(k*pi*x/L)
#      }
# return(summ)
# }
# 
# 
# for(ix in 1:(N+1))
#   for(it in 1:n) {
#     #Uacc[it,ix]=A(tj[it]) + (B(tj[it]) - A(tj[it])) * x[ix]/L - 8*(B(tj[it])-A(tj[it]))/pi^3 * summm(x[ix], tj[it])
#     Uacc[it,ix]=B(tj[it])*(a^2*tj[it]/L+(3*x[ix]^2-L^2)/(6*L))+2*L/pi^2*summm(x[ix], tj[it])
# }


options(scipen = 999) # Disable exponential notation (e.g. 1.81e+09)
print("Численное решение")
print(U)
#print("Точное решение")
#print(Uacc)

library(animation)
oopt = ani.options(interval = 0.3)

ani.record(reset = TRUE)
for(j in 1:n)
{
  plot(U[j,], xaxt="n", xlab = 'X values', ylab = 'U[x,t] values')
  lines(U[j,],col="red")
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