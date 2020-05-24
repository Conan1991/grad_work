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
h = 0.00045 #/ 2#шаг по x
tau = h^2/a^2 #/ 4

L = 0.008 #длина
N = round(L/h) #Число шагов
n = 100 #* 4 #Число шагов по времени

x = numeric(N+1)


for(i in 1:N+1)
  x[i]= (i-1)*h

#модуль объёмного раширения постоянная Ламе
lambda =  0.1
#коэффициент линейного теплового расширения
alpha = 0.0033# для металлов ~10^-4

nu = 0.25 #Модуль сдвига постоянная Ламе

epsilon = 0.01
number_of_iterations = 2000


step_x = 0.00045#/2
step_y = 0.00045#/2
step_t = (step_x^2/lambda)/100#*4

