# Jeff Rouder, June, 2025
# Models in this file

# facManMarg. Marginal implementation of manifest model
#   input is a matrix of scores, rows are individual, columns are tasks

# facManCond. Conditional implementation of manifest model
#   input is a matrix of scores, rows are individual, columns are tasks

# facH. Conditional implementation of hierarchical model with a contrast.
#   input is dataframe with sub, task, cond, y

# facH0cfa. Conditional implementation of hierarchical model without contrast
#   input is dataframe with sub, task, y
#   must provide a matrix free that indicates which loadings are to be estimated
#   use with care to prevent rotation / mirror problems.

# cov. Conditional implementation of hierarchical model with a contrast.
#   input is dataframe with sub, task, cond, y
#   output is Sigma, prior is scaled inverse Wishert




require(R2jags)

facManMarg.mod="
model{
  for (j in 1:J){
    mu[j]~dnorm(mu.m, pow(mu.s, -2))
    pDelta2[j]~dgamma(.5,.5*pow(tuneDelta,2))
    delta2[j] = 1/pDelta2[j]
    for (d in 1:D){
      lambda[j,d] ~ dnorm(0, pow(tuneLambda, -2))}} 
  for (jx in 1:J) {
    for (jy in 1:J){
      DM[jx, jy] = equals(jx,jy) * delta2[jx]
      cross_product[jx, jy] <- inprod(lambda[jx,], lambda[jy,])}}
  Sigma = cross_product + DM
  for (i in 1:I){    
    y[i,1:J] ~ dmnorm.vcov(mu, Sigma)}
  }"

#Prior is list with mu.m, mu.s, tuneDelta, tuneLambda
facManMarg.run=function(y,numFactors,M=500,prior){
  I=dim(y)[1]  
  J=dim(y)[2]
  setup = list(
    "y" = y,
    "I" = I,
    "J" = J,
    "D" = numFactors)
  pars = c("Sigma","lambda", "mu","delta2")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(facManMarg.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}


facManCond.mod="
model{
  for (j in 1:J){
    mu[j]~dnorm(mu.m, pow(mu.s, -2))
    pDelta2[j]~dgamma(.5,.5*pow(tuneDelta,2))}
  for (d in 1:D){
    for (j in 1:J){
      lambda[j,d] ~ dnorm(0, pow(tuneLambda, -2))}
    for (i in 1:I){
      eta[i,d] ~ dnorm(0,1)}}
    
  for (i in 1:I){
    for (j in 1:J){
      for (d in 1:D){
        temp[i,j,d] <- lambda[j,d]*eta[i,d]}
      center[i,j] <- mu[j] + sum(temp[i,j,1:D])
      y[i,j] ~dnorm(center[i,j],pDelta2[j])}}
}"


facManCond.run=function(y,numFactors,M=500,prior){
  I=dim(y)[1]  
  J=dim(y)[2]
  setup = list(
    "y" = y,
    "I" = I,
    "J" = J,
    "D" = numFactors)
  pars = c("lambda", "eta","mu","pDelta2")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(facManCond.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}


facH.mod="
model{

  for (j in 1:J){
    pTau2[j] ~ dgamma(.5, .5*pow(tuneTau,2))
    mu[j]~dnorm(mu.m, pow(mu.s, -2))
    pDelta2[j]~dgamma(.5,.5*pow(tuneDelta,2))}
  for (d in 1:D){
    for (j in 1:J){
      lambda[j,d] ~ dnorm(0, pow(tuneLambda, -2))}
    for (i in 1:I){
      eta[i,d] ~ dnorm(0,1)}}
    
  for (i in 1:I){
    for (j in 1:J){
      for (d in 1:D){
        temp[i,j,d] <- lambda[j,d]*eta[i,d]}
      centerTheta[i,j] <- mu[j] + sum(temp[i,j,1:D])
      theta[i,j] ~dnorm(centerTheta[i,j],pDelta2[j])
      alpha[i,j] ~dnorm(alpha.m,pow(alpha.s, -2))}}
  for (n in 1:N){
    center[n] = alpha[sub[n],task[n]]+(cond[n]-1.5)*theta[sub[n],task[n]]
    y[n] ~ dnorm(center[n], pTau2[task[n]])}
}"


facH.run=function(dat,numFactors,M=500,prior){
  setup = list(
    "y" = dat$y,
    "sub" = dat$sub,
    "task" = dat$task,
    "cond" = dat$cond,
    "I" = max(dat$sub),
    "J" = max(dat$task),
    "N" = length(dat$y),
    "D" = numFactors)
  pars = c("lambda", "eta","mu","pDelta2","theta","pTau2")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(facH.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}


# free = 0, 1, 2: 0=fixed to zero, 1=estimate as real, 2=estimate as positive

facH0cfa.mod="
model{
  for (j in 1:J){
    pTau2[j] ~ dgamma(.5, .5*pow(tuneTau,2))
    mu[j]~dnorm(mu.m, pow(mu.s, -2))
    pDelta2[j]~dgamma(.5,.5*pow(tuneDelta,2))}
  for (d in 1:D){
    for (j in 1:J){
      prec[j,d]=pow(tuneLambda,-2)*ifelse(free[j,d],1,1e8)
      lambda[j,d] ~ dnorm(0, prec[j,d]) T(lower[pos[j,d]],)}
    for (i in 1:I){
      eta[i,d] ~ dnorm(0,1)}}
    
  for (i in 1:I){
    for (j in 1:J){
      for (d in 1:D){
        temp[i,j,d] <- lambda[j,d]*eta[i,d]}
      centerTheta[i,j] <- mu[j] + sum(temp[i,j,1:D])
      theta[i,j] ~dnorm(centerTheta[i,j],pDelta2[j])}}
  for (n in 1:N){
    y[n] ~ dnorm(theta[sub[n],task[n]], pTau2[task[n]])}
}"

# free = 0, 1, 2: 0=fixed to zero, 1=estimate as real, 2=estimate as positive

facH0cfa.run=function(dat,M=500,prior,free){
  pos=ifelse(free>1,2,1)
  setup = list(
    "y" = dat$y,
    "sub" = dat$sub,
    "task" = dat$task,
    "I" = max(dat$sub),
    "J" = max(dat$task),
    "N" = length(dat$y),
    "D" = dim(free)[2],
    "free" = free>0,
    "pos" = pos,
    "lower"=c(-1,0))
  pars = c("lambda", "eta","mu","pDelta2","theta","pTau2")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(facH0cfa.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}


cov.mod ="
model{
  for (j in 1:J){
    pTau2[j] ~ dgamma(.5, .5*(tau.s^2))
    mu[j]~dnorm(mu.m, pow(mu.s, -2))}
  for (i in 1:I){
    theta[i,1:J] ~ dmnorm(mu, Sigma)}
  for (i in 1:I){
    for (j in 1:J){
      alpha[i,j] ~ dnorm(alpha.m,pow(alpha.s,-2))}}
  for (n in 1:N){
    center[n] = alpha[sub[n],task[n]]+(cond[n]-1.5)*theta[sub[n],task[n]]
    y[n] ~ dnorm(center[n], pTau2[task[n]])}
  Sigma~dscaled.wishart(S,2)
}
"
cov.run <- function(dat,M=200,prior){
  I=length(unique(dat$sub))
  J=length(unique(dat$task))
  setup=list(
    "y" = dat$y,
    "task" = dat$task,
    "sub" = dat$sub,
    "cond" = dat$cond,
    "I" = I,
    "J" = J,
    "N" = nrow(dat))
  pars = c("Sigma","mu","theta")
  out=jags(data=c(setup,prior), 
           parameters=pars, 
           model.file = textConnection(cov.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)  
  return(out)}