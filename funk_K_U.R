N= 100
k1 = 1
k2 = 3
kj = numeric(N)

for(i in 1:N)
{
  kj[i] = k1+ i/N*(k2-k1)
}
plot(kj)