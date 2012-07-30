##Moments
##k1 = immigration rate, k2 =death
E = function(k1, k2, n0, t) {
    (k1-k1*exp(-k2*t)+n0*exp(-k2*t)*k2)/k2
}

V = function(k1, k2, n0, t) {
    (-n0*exp(-2*k2*t)*k2+k1-k1*exp(-k2*t)+n0*exp(-k2*t)*k2)/k2
}

##First-order 
dE_dk1 = function(k1, k2,n0, t) {
    (1-exp(-k2*t))/k2    
}

dE_dk2 = function(k1, k2,n0, t) {
    -k1*(1-exp(-k2*t))/k2^2+k1*t*exp(-k2*t)/k2-n0*t*exp(-k2*t)   
}

##Derivative of Variance wrt to theta
dV_dk1 = function(k1, k2,n0, t) {
    (1-exp(-mu*t))/mu
}
dV_dk2 = function(k1, k2,n0, t) {
    (2*n0*t*exp(-2*mu*t)*mu-n0*exp(-2*mu*t)+alpha*t*exp(-mu*t)-n0*t*exp(-mu*t)*mu+n0*exp(-mu*t))/mu-(-n0*exp(-2*mu*t)*mu+alpha-alpha*exp(-mu*t)+n0*exp(-mu*t)*mu)/mu^2
}

Firstorder = function(k1, k2, t, x) {
    value = 0
    x = dd$n
    t = dd$times
    i=1;k1=1;k2=0.1
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        
        m = E(k1, k2, n0, tau)
        ##For k1
        part1 = dE_dk1(k1,k2,n0,tau)*1/V(k1,k2,n0,tau)*(m-x[i+1]) 
        part1 = part1 + (m-x[i+1])+d1V_dk1(k1,k2,n0,tau) *(m-x[i+1])
        
        ##For k2
        part2 = dE_dk2(k1,k2,n0,tau)*1/V(k1,k2,n0,tau)*(m-x[i+1])
        part2 = part2 + (m-x[i+1])+d1V_dk2(k1,k2,n0,tau) *(m-x[i+1])
        value = value + part1 + part2
    }
    return(value)
}

dd = read.csv("data/immigration_death.csv")
Firstorder(1, 0.1, dd$times, dd$n)
