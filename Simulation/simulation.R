##Discretise a fully observed stochastic simulation
discretise = function(dd, tstep) {
    time = tstep
    dd_dis = dd[1:(max(dd$time)/tstep+1),]
    j = 2
    for(i in 1:(nrow(dd)-1)) {
        while(dd$time[i] < time && dd$time[i+1] > time) {
            dd_dis[j,] = c(time, dd[i, 2:ncol(dd)])
            time = time + tstep
            j = j + 1
        }
    }
    return(dd_dis)
}

gillespie = function(IC, pars, maxtime=10) {
    time = 0; tstep=0; i=1 
    X1 = IC[1]; X2 = IC[2]
    
    i = 1
    pop1 = numeric(1000); pop2 = numeric(1000); times = numeric(1000)
    while(time < maxtime) {
        
        #Update state vectors
        pop1[i] = X1; pop2[i] = X2; times[i] = time
        i = i + 1
        #Calculate the overall rate
        total_rate = pars[1] + pars[2]*X1
        time = time + rexp(1, total_rate)

        #Don't go by maxtime
        if(time >= maxtime)
            break;
        
        #Choose which reaction happens
        u = runif(1)
        if(u<pars[1]/total_rate) {
            X1 = X1 + 1
        }else {
            X1 = X1 - 1
            X2 = X2 + 1
        }
    }
    #Update state vectors
    pop1[i] = X1; pop2[i] = X2; times[i] = time
    dd = data.frame(times, n=pop1, c=pop2)
    return(dd[1:i,])
}

set.seed(1)
dd = gillespie(c(0, 0), c(1, 0.1), 50)
dd_dis = discretise(dd, 5)

##Plot the data
par(mfrow=c(1,2))
plot(dd$times, dd$n, type='s')
points(dd_dis$times, dd_dis$n, col=2)    
    
plot(dd$times, dd$c, type='b')
points(dd_dis$times, dd_dis$c, col=2)

##Output data as a csv file
write.table(dd, file="sim.csv", sep=",",row.names=FALSE, col.names=FALSE)
