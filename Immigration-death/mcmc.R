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

dd = read.csv("data/immigration_death.csv")
log_likelihood(k1, k2, dd[,1], dd[,2])

k1 = seq(0,3, 0.005)
k2 = seq(0,0.4, 0.005)
d = outer(k1, k2, function(k1, k2) log_likelihood(k1, k2, dd[,1], dd[,2]))

image(k1, k2, exp(d))
# 
# prop_theta = function(theta_cur){
#     theta_cur
# }
# 
# 
# 
# 
# 
# theta_cur = c(0.1, 0.1)
# ll_cur = -Inf
# for(i in 1:iters) {
#     theta_prop  = prop_theta(theta_cur)
#     ll_prop = log_likelihood(theta_prop)
#     
#     
#     
#     
#     
#     
#     
# }
# 
# 
