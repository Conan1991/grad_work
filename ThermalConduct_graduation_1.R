rm(list = ls())


a=1
tau = 0.01
h = 0.025 #шаг по x
L = 1 #длина
N = L/h #Число шагов

k1=k2=1
A = function(t)
{
  #0
  #2/(1+t)
  2
}
B = function(t)
{
  #0
  #2/(1+t)
  2
}

n = 15 #Число шагов по времени
tj = numeric(n) #tau j
nj = numeric(n) #ветор моментов времени


x = numeric(N+1)
fij = numeric(N+1) #сила

Ai = numeric(N) #коэффициенты
Bi = numeric(N)
Ci = numeric(N)
Fi = numeric(N-1)


alpha = numeric(N)
beta = numeric(N)

k_U = numeric(N)
C_U = numeric(N)

# k2_U= 0.5786
# k1_U= 2.3
# 
# c2_U = 4195
# c1_U = 2000

k2_U= 2.3
k1_U= 2.3

c2_U = 2000
c1_U = 2000


c2=10
c1=-10


#j=1 0-й слой
kU = function(U)
{
  if( U < 0)
    return(k1_U)
  return(k2_U)
}

CU = function(U)
{
  if( U < 0)
    return(c1_U)
  return(c2_U)
}



for(j in 1:n)
{
  tj[j]=(j-1)*tau
  nj[j]=j-1
}

Ux0=function(x)
{
  #A(tj[1]) + (B(tj[1])-A(tj[1]))*(x/L)^2;
  #cos(pi*x)
  1
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

U0= U[1,]
CoeffF = function(j) #Поправка на коэфф
{
  
alpha[1] = k1
beta[1]= - h * A(tj[j])

for(i in 1:N) # Считаем очередные Ai, Bi , Ci
{
  Ai[i] = - kU(U[j,i]) / h^2
  Bi[i] = - kU(U[j,i]) / h^2
  Ci[i] = CU(U[j,i]) / tau + 2 * kU(U[j,i]) / h^2
}

for(i in 1:(N-1)) #Считаем Fi
Fi[i] = U[j-1,i+1] / tau
  
for(i in 2:N) #Считаем альфа и бета коэффициенты
{
alpha[i]=-Bi[i-1]/(Ci[i-1]+Ai[i-1]*alpha[i-1])
beta[i]= (Fi[i-1]-Ai[i-1]*beta[i-1])/(Ci[i-1]+Ai[i-1]*alpha[i-1])
}

U[j,N+1]=(h*B(tj[j])+k2*beta[N])/(1-k2*alpha[N])
for(i in N:1) #Считаем Uj
  U[j,i]= alpha[i] * U[j,i+1] + beta[i]

return(U[j,])
}

IterF = function(j) #Считаем остальные слои
{
  
  alpha[1] = k1
  beta[1] = - h * A(tj[j])
  
  for(i in 1:N) # Считаем очередные Ai, Bi , Ci
  {
    Ai[i] = -kU(U[j-1,i])/h^2
    Bi[i] = -kU(U[j-1,i])/h^2
    Ci[i] = CU(U[j-1,i])/tau + 2*kU(U[j-1,i])/h^2
  }
  
  for(i in 1:(N-1)) #Считаем Fi
    Fi[i]=U[j-1,i+1]/tau
  
  for(i in 2:N) #Считаем альфа и бета коэффициенты
  {
    alpha[i]=-Bi[i-1]/(Ci[i-1]+Ai[i-1]*alpha[i-1])
    beta[i]= (Fi[i-1]-Ai[i]*beta[i-1])/(Ci[i-1]+Ai[i-1]*alpha[i-1])
  }
  
  U[j,N+1]=(h*B(tj[j])+k2*beta[N])/(1-k2*alpha[N])
  for(i in N:1) #Считаем Uj
    U[j,i]= alpha[i] * U[j,i+1] + beta[i]
  
  return(U[j,])
}

for(j in 2:n)
{
  U[j,] = IterF(j)
  #temp = CoeffF(j)
  #U[j,] = temp
}



options(scipen=999)

# Uacc = matrix(data=NA,nrow=n,ncol=N+1)
# colnames(Uacc)=x
# row.names(Uacc)= c(0:(n-1))
# 
# summm=function(x,t)
# {
# summ=0
#      for(k in 1:10000)
#      {
#         summ= summ + 1/(2*k-1)^3*sin((2*k-1)*pi*x/L) * exp(-(2*k-1)^2*pi^2*a^2*t/L^2)
#         #if(k==10) print(summ)
#      }
# return(summ)
# }


# for(ix in 1:(N+1))
#   for(it in 1:n) 
#     Uacc[it,ix]=A(tj[it]) + (B(tj[it]) - A(tj[it])) * x[ix]/L - 8*(B(tj[it])-A(tj[it]))/pi^3 * summm(x[ix], tj[it])


library(animation)
oopt = ani.options(interval = 0.3)

ani.record(reset = TRUE)
for(j in 1:n)
{
  plot(U[j,], xaxt="n", xlab = 'X values', ylab = 'U[x,t] values')
  lines(U[j,],col="red")
  lines(U0,col="blue")
  axis(1, at = c(1:(N+1)), labels = x)
  ani.pause()
  ani.record()
  
}



#print("Численное решение")
print(U)
for(j in 1:n)
{
plot(U[j,], xaxt="n", xlab = 'X values', ylab = 'U[x,t] values')
lines(U[j,],col="red")
lines(U0,col="blue")
}


while(TRUE)
{ 
ani.replay()
ani.options(oopt, loop = 100)
}
