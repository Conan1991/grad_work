k1_U = 2.3 #Теплопроводность льда
#k2_U = k1_U
k2_U = 0.58 #Теплопроводность воды

density_1  = 918.7 #Плотность льда
density_2 = 999.7 #Плотность воды
#mean_density= mean(c(density1,density2))
heat_capacity1 = 2000 #Теплоёмкость льда
heat_capacity2 = 4195 #Теплоёмкость воды

c1_U = heat_capacity1 * density_1 
c2_U = heat_capacity2 * density_2 
#c2_U = c1_U

mean_K = mean(c(k1_U,k2_U))
mean_C = mean(c(c1_U,c2_U))

a=(mean_K/mean_C)^(1/2)
h = 0.00045#/4#шаг по x
tau = h^2/a^2 #/ 16

L = 0.008 #длина
N = round(L/h) #Число шагов
n = 100 #*16 #Число шагов по времени

x = numeric(N+1)


for(i in 1:N+1)
  x[i]= (i-1)*h