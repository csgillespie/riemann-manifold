##Prob(x=n(t)|n0=n(0)) for immigration death model
pn = function (k1, k2, t, n0, x) {
  k3 = k1/k2*(1-exp(-k2*t))
  value = 0
  for(i in 0:x) {
      value = value + choose(n0, i)*
          exp(-k2*t*i)*
          (1-exp(-k2*t))^(n0-i)*
          dpois(x-i, k3)
  }
  return(value)
}


