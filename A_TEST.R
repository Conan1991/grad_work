#source('shifts.R')

l = 20
k = round(20/2)

source('temperature.R')

temperature = ThermalConduct()

#v_point <- temperature[2,]
#list_p = list(v_point)
#list_p[[2]] <- c(1,0)
#old_matrix = list()

#old_matrix[[1]] = list_p
shift = function(vecl , counter)
{
  copy = c()
  for(i in 1:(length(vecl) - counter))
  {
    copy[i] = vecl[i+counter]
  }
  return(copy)
}

checkMelting = function(iter)
{
  x_truncated = temperature[iter,]
  #copy = x_truncated
  counter = 0
  for(i in 1 : length(x_truncated))
    {
      if(x_truncated[i] > 0)
        counter=counter+1
  }
  copy = shift(x_truncated, counter)
  
  return (copy)
}
res= checkMelting(2)


checkColor = function(percent)
{
  red = 0
  green = 0
  blue = 0
  const = 255 / 50
  
  if(percent >= 50)
  {
    red = 255
    green = 255 - round((percent-50)*const)
    if(green <= 0)
      green = 0
  }
  else
  {
    green = 255
    red = 255 - round((50 - percent)*const)
    if(red <= 0)
      red = 0
  }
  vec = c(red, green, blue)
  
  return(vec)
}

color = checkColor(15)


