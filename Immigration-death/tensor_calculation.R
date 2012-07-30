##Moments
##k1 = birth, k2 =death
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

##Second-order
##Derivative of 1/Variance wrt to theta
d1V_dk1 = function(k1, k2,n0, t) {
    -mu*(1-exp(-mu*t))/(-n0*exp(-2*mu*t)*mu+alpha-alpha*exp(-mu*t)+n0*exp(-mu*t)*mu)^2
}
d1V_dk1 = function(k1, k2,n0, t) {
    -mu*(2*n0*t*exp(-2*mu*t)*mu-n0*exp(-2*mu*t)+alpha*t*exp(-mu*t)-n0*t*exp(-mu*t)*mu+n0*exp(-mu*t))/(-n0*exp(-2*mu*t)*mu+alpha-alpha*exp(-mu*t)+n0*exp(-mu*t)*mu)^2+1/(-n0*exp(-2*mu*t)*mu+alpha-alpha*exp(-mu*t)+n0*exp(-mu*t)*mu)
}

Firstorder = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        m = E(k1, k2, n0, tau)
        tau = t[i+1] - t[i]
        n0 = x[i]
        
        ##For k1
        value = value + dE_dk1(k1,k2,n0,tau)*1/V(k1,k2,n0,tau)*(m-x[i+1])
        value = value + (m-x[i+1])+d1V_dk1 *(m-x[i+1])
        
        ##For k2
        value = value + dE_dk2(k1,k2,n0,tau)*1/V(k1,k2,n0,tau)*(m-x[i+1])
        value = value + (m-x[i+1])+d1V_dk2 *(m-x[i+1])
    }
    return(value)
}

dd = read.csv("sim.csv", header=FALSE)

