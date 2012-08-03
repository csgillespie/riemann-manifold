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
    (1-exp(-k2*t))/k2
}
dV_dk2 = function(k1, k2,n0, t) {
    (2*n0*t*exp(-2*k2*t)*k2-n0*exp(-2*k2*t)+k1*t*exp(-k2*t)-n0*t*exp(-k2*t)*k2+n0*exp(-k2*t))/k2-(-n0*exp(-2*k2*t)*k2+k1-k1*exp(-k2*t)+n0*exp(-k2*t)*k2)/k2^2
}

##Derivative of 1/Variance wrt to theta
d1V_dk1 = function(k1, k2, n0, t){
    v = V(k1,k2,n0,tau)
    -1/v*dV_dk1(k1,k2,n0,tau)/v
}

d1V_dk2 = function(k1, k2, n0, t){
    v = V(k1, k2, n0, tau)
    -1/v*dV_dk2(k1,k2,n0,tau)/v
}
#########################################################


##k1 = immigration rate, k2 =death
E = function(k1, k2, n0, t) {
    (k1 - k2*n0)*t
}

V = function(k1, k2, n0, t) {
    (k1 + k2*n0)*t
}

##First-order 
dE_dk1 = function(k1, k2,n0, t) {
    t
}

dE_dk2 = function(k1, k2,n0, t) {
    n0*t
}

##Derivative of Variance wrt to theta
dV_dk1 = function(k1, k2,n0, t) {
    t
}
dV_dk2 = function(k1, k2,n0, t) {
    n0*t
}

##Derivative of 1/Variance wrt to theta
d1V_dk1 = function(k1, k2, n0, t){
    v = V(k1,k2,n0,tau)
    -1/v*dV_dk1(k1,k2,n0,tau)/v
}

d1V_dk2 = function(k1, k2, n0, t){
    v = V(k1, k2, n0, tau)
    -1/v*dV_dk2(k1,k2,n0,tau)/v
}








##Work out Ll
L = function(k1, k2, t, x) {
    n0 = x[1:10]
    n1 = x[2:11] 
    means = E(k1, k2, n0, 5)
    vars = V(k1,k2,n0, 5)
    
    value = 0 ;i =1
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        value = value + log(dnorm(x[i+1], means[i], sqrt(vars[i])))
    }
    return(value)
}

##First-order derivative
L_deriv_k1 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        ##For k1
        part = dE_dk1(k1,k2,n0,tau)*1/v*(x[i+1]-m) 
        part = part - 0.5*(x[i+1]-m)*d1V_dk1(k1,k2,n0,tau) *(x[i+1]-m)
        part = part - 0.5*1/v*dV_dk1(k1,k2,n0,tau)
        value = value + part 
    }
    return(value)
}

L_deriv_k2 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        ##For k2
        part = dE_dk2(k1,k2,n0,tau)*1/v*(x[i+1]-m)
        part = part - 0.5*(x[i+1]-m)*d1V_dk2(k1,k2,n0,tau) *(x[i+1]-m)
        part = part - 0.5*1/v*dV_dk2(k1,k2,n0,tau)
        value = value + part
    }
    return(value)
}

L_deriv = function(k1, k2, t, x) {
    c(L_deriv_k1(k1,k2, t, x), L_deriv_k2(k1,k2, t, x))
}
    
#Numerical check - very handy!
# epsilon = 1e-6
# (L(k1+epsilon, 0.1, dd$times, dd$n)- L(k1, k2, dd$times, dd$n))/epsilon
# (L(k1, k2 + epsilon, dd$times, dd$n)- L(k1, k2, dd$times, dd$n))/epsilon
# L_deriv_k1(1, 0.1, dd$times, dd$n)
# L_deriv_k2(1, 0.1, dd$times, dd$n)


##Second-order derivative
L_deriv2_k1k1 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        #For k1
        part = dE_dk1(k1,k2,n0,tau)*1/v*dE_dk1(k1,k2,n0,tau)
        part = part + 1/2*1/v*dV_dk1(k1,k2,n0,tau)*1/v*dV_dk1(k1,k2,n0,tau)
        value = value + part
    }
    return(value)    
}

L_deriv2_k1k2 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        
        part = dE_dk1(k1,k2,n0,tau)*1/v*dE_dk2(k1,k2,n0,tau)
        part = part + 1/2*1/v*dV_dk1(k1,k2,n0,tau)*1/v*dV_dk2(k1,k2,n0,tau)
        value = value + part
    }
    return(value)    
}

L_deriv2_k2k1 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        
        part = dE_dk2(k1,k2,n0,tau)*1/v*dE_dk1(k1,k2,n0,tau)
        part = part + 1/2*1/v*dV_dk2(k1,k2,n0,tau)*1/v*dV_dk1(k1,k2,n0,tau)
        value = value + part
    }
    return(value)    
}

L_deriv2_k2k2 = function(k1, k2, t, x) {
    value = 0
    for(i in 1:(length(x)-1)) {
        tau = t[i+1] - t[i]
        n0 = x[i]
        m = E(k1, k2, n0, tau)
        v = V(k1,k2,n0,tau)
        
        part = dE_dk2(k1,k2,n0,tau)*1/v*dE_dk2(k1,k2,n0,tau)
        part = part + 1/2*1/v*dV_dk2(k1,k2,n0,tau)*1/v*dV_dk2(k1,k2,n0,tau)
        value = value + part
    }
    return(value)    
}

L_deriv2 = function(k1, k2, t, x) {
    c(L_deriv2_k1k1(k1, k2, t, x), 
      L_deriv2_k1k2(k1, k2, t, x), 
      L_deriv2_k2k1(k1, k2, t, x),
      L_deriv2_k2k2(k1, k2, t, x))
}
    

dd = read.csv("data/immigration_death1.csv")
k1 = 1; k2 = 0.1
L(k1, k2, dd$times, dd$n)
L_deriv(k1, k2, dd$times, dd$n)
L_deriv2(k1, k2, dd$times, dd$n)

L_deriv2_k1k1()

1/100
(1.3-0.07)/20

#k1 = seq(0.5,1.5, length.out=100)
k1 = seq(0.5, 1.5, 0.05)
k2 = seq(0.07, 1.3, 0.0615)
d11 = outer(k1, k2, function(k1, k2) L_deriv2_k1k1(k1, k2, dd[,1], dd[,2]))
d12 = outer(k1, k2, function(k1, k2) L_deriv2_k1k2(k1, k2, dd[,1], dd[,2]))
d21 = outer(k1, k2, function(k1, k2) L_deriv2_k2k1(k1, k2, dd[,1], dd[,2]))
d22 = outer(k1, k2, function(k1, k2) L_deriv2_k2k2(k1, k2, dd[,1], dd[,2]))

setnicepar(mfrow=c(1,2))
image(d11)


setnicepar(mfrow=c(2,2))
image(k1, k2, d11)
image(k1, k2, d12)
image(k1, k2, d21)
image(k1, k2, d22)
#
exact = read.csv("data/tmp.csv", header=FALSE)[1:441,]
k1 = rep(seq(0.5, 1.5, 0.05), each=21)
k2 = rep(seq(0.07, 1.3, 0.0615), 21)
image(seq(0.5, 1.5, 0.05), seq(0.07, 1.3, 0.0615), matrix(exact, nrow=21))


