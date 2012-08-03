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

log_likelihood = function(k1, k2, t, x) {
    ll = 0
    #k1 = theta[1];k2 = theta[2]
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        ll = ll + log(pn(k1, k2, tau, n0, x[i+1]))
    }
    return(ll)
}

dd = read.csv("data/immigration_death_k1_1_k2_1.csv")
log_likelihood(1, 0.1, dd[,1], dd[,2])

k1 = seq(log10(1.1),log10(100),  length.out=500)
k2 = seq(log10(1.1),log10(100), length.out=500)
d = outer(k1, k2, function(k1, k2) log_likelihood(k1, k2, dd[,1], dd[,2]))
image(log10(k1), log10(k2), exp(d))
points(1, 1,pch=18)
