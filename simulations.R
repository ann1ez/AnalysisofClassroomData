set.seed(2021)
library(tidyverse)

#defining the function 
#for a simulation of the teacher process
sim_aw <- function(time,mu) {
  
  path <- matrix(0, nrow = 1, ncol = 2)
  
  jumps_number <- rpois(1, lambd = mu * time)
  jumps_time <- runif(n = jumps_number, min = 0, max = time) %>% sort()
  
  jumps_time
}

#defining the function 
#for a simulation of the student process
two.pois.sim2 <- function(time, mu, alpha0, gamma){
  
  tset <- sim_aw(time, mu)
  M <- length(tset)
  w<-rep(0,M+2)
  if(M>0){
    u.diff = c(tset, time)- c(0, tset)
    w[2:(M+2)] = exp(alpha0)*cumsum(exp(gamma*0:M)*u.diff)
    Ftinv<- function(u){
      k = findInterval(u, w/w[M+2])+1
      (w[M+2]*u- w[k-1])*exp(-alpha0-gamma*(k-2))+c(0, tset)[k-1]
    }
  } else{
    tset = NA
    Ftinv = function(u){time *u}
    w[2] = time*exp(alpha0)
  }
  
  N=rpois(1, w[M+2])
  if(N>0){
    X0=rep(NA,N)
    for(i in 1:N){
      X0[i]=Ftinv(runif(1))
    }
  }else{
    X0= NA
  }
  list(t.times = tset, s.times = if(is.numeric(X0)){sort(X0)}else{NA})
}


library(rootSolve)
#defining a function
#for calculating the estimated gamma value
mle.sim.gamma<- function(gamma, n.sim =50,time= 50, mu=0.5, alpha0=0){
  sim.res =  matrix(0, nrow = n.sim, ncol = 6)
  colnames(sim.res)<- c("time", "mu", "alpha0", "gamma","alpha0hat","gammahat")
  param <- c(time, mu, alpha0, gamma)
  u.diff = list(0)
  r = rep(0, 20)
  N= rep(0, 20)
  M = rep(0, 20)
  for(i in 1:n.sim){for( p in seq(from=1, to=20, by=1)){
    two.pois.sim2(time, mu, alpha0, gamma)-> res
    if(is.numeric(res$t.times)&is.numeric(res$s.times) ){
      r[p]= sum(findInterval(res$s.times, c(0, res$t.times, 50))-1)
      N[p]= length(res$s.times)
      M[p] = length(res$t.times)
      u.diff[[p]] = c(res$t.times, time)- c(0, res$t.times) 
    }else{
      r[p] = NA  
    } 
  }
    if(sum(is.na(r))<1){     
      model=function(par){
        l.prime.adp = 0 #alpha partial derivative
        l.prime.gdp = 0 #gamma partial derivative
        for (p in 1:20){   
          l.prime.adp=(-exp(par[1])*sum(exp(par[2]*0:M[p])*u.diff[[p]]) +N[p]) + l.prime.adp
          l.prime.gdp=(-exp(par[1])*sum((0:M[p])*exp(par[2]*0:M[p])*u.diff[[p]]) +r[p]) + l.prime.gdp }
        c(l.prime.adp=l.prime.adp,l.prime.gdp=l.prime.gdp) }
      p3<- multiroot(model,c(-2,2))$root[1]
      p5<-multiroot(model,c(-2,2))$root[2]
      
      sim.res[i,] = c(param, c(p3,p5)) } else {
        sim.res[i,] = c(param, rep(NA,2))   
      }
  }
  return(sim.res)
}

#running the simulations
#500 simulations for each gamma, aplha pair
sim.multi <- mle.sim.gamma(gamma=0, n.sim=500, time=50, mu=0.5, alpha0=0)
for( i in seq(from=0.06, to=.1, by=0.01)){ for( j in seq(from=-2.0, to=-1.3, by=0.1)){
  sim.multi= rbind(sim.multi,mle.sim.gamma(i,500,50,0.5,j))
  print(i)
  }
}

#organizing the data into a metadata
summultimeta<-as.data.frame(sim.multi)
my.mean = function(x){mean(x, na.rm = TRUE)}
meta6_15<-as.data.frame(sim.multi)%>%mutate(g.diff = gammahat - gamma, a.diff = alpha0hat-alpha0)%>%
  group_by(gamma,alpha0)%>% 
  summarise(cnt = sum(!is.na(gammahat)), g.bias = my.mean(g.diff),   a.bias = my.mean(a.diff), g.mse = my.mean(g.diff^2) , a.mse= my.mean(a.diff^2))

write.table(meta6_15,file = "metadatafor6_15.txt", sep = "\t")
write.table(sim.multi,file = "simmultifor06_15.txt", sep = "\t")