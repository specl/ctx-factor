source('facMod.R')
source('facExtra.R')

fitEnkavi=function(dat,prior,numFactors){
  dat$oTask=dat$task
  dat$task=as.integer(as.factor(dat$task))
  dat$sub=as.integer(as.factor(dat$sub))
  dat$y=1000*dat$rt
  out=facH.run(dat,numFactors=numFactors,M=5000,prior=prior)
  samples=out$BUGSoutput$sims.list
  aligned=samples
  if (numFactors>1) aligned=makePositive(align(samples))
  return(aligned)
}

fitEnkaviCov=function(dat,prior){
  dat$oTask=dat$task
  dat$task=as.integer(as.factor(dat$task))
  dat$sub=as.integer(as.factor(dat$sub))
  dat$y=1000*dat$rt
  out=cov.run(dat,M=5000,prior=prior)
  return(out$BUGSoutput$sims.list)
}

myCorPlot=function(cormat,limits=c(-1,1),main=""){
  J=dim(cormat)[1]
  cn <- paste("T",1:J,sep="")
  rownames(cormat)=rep("",J)
  colnames(cormat)=rep("",J)
  corrplot(cormat,method="color",
           cl.pos='n',col.lim=limits,main=main,
           mar=c(0,0,1,0))}


enkaviDecomp=function(std){
  
  M=dim(std$lambda)[1]
  J=dim(std$lambda)[2]
  prop=array(dim=c(5,M,J,J))
  prop[1,,,]=std$rho
  for (m in 1:M){
    prop[2,m,,]=crossprod(t(std$lambda[m,,1]))
    prop[3,m,,]=crossprod(t(std$lambda[m,,2]))
    prop[4,m,,]=crossprod(t(std$lambda[m,,3]))
    prop[5,m,,]=diag(std$delta2[m,])
  }
  return(prop)
}