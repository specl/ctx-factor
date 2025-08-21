require(corrplot)
require(mvtnorm)
require(lavaan)
source('facMod.R')
source('facExtra.R')

makeFlatDesign=function(I,J,L){
  sub=rep(1:I,each=J*2*L)
  task=rep(rep(1:J,each=2*L),I)
  cond=rep(rep(1:2,each=L),I*J)
  return(data.frame(sub,task,cond))}

makeAlpha=function(I,J,sigVal=200){
  mu = rep(800,J)
  cor=matrix(nrow=J,ncol=J,.5)+.5*diag(J)
  sig= rep(sigVal,J)
  cov=diag(sig)%*%cor%*%diag(sig)
  vals=rmvnorm(I,mu,cov)
  return(vals)}

addObs=function(dat,truth){
  subtask=cbind(dat$sub,dat$task)
  m=truth$alpha[subtask]+(dat$cond-3/2)*truth$theta[subtask]
  dat$y=round(rnorm(length(m),m,sqrt(truth$tau2[dat$task])),1)
  return(dat)
}

makeTheta=function(I,tMu,tLambda,tEta,tDel2){
  I=dim(tEta)[1]
  J=dim(tLambda)[1]
  center=t((tLambda)%*%t(tEta)+tMu)
  tTheta=matrix(nrow=I,ncol=J,rnorm(I*J,center,outer(rep(1,I),sqrt(tDel2))))
  return(tTheta)}

makeTruth=function(I=200,J=8,trialNoise=200){
  D=2
  up=30
  down=7
  tLambda=matrix(nrow=J,ncol=D,0)
  tLambda[,1]=seq(up,down,length=J)
  tLambda[,2]=seq(down,up,length=J)
  tEta=matrix(nrow=I,ncol=D,rnorm(I*D))
  tDel2=rep(20^2,J)
  tMu=rep(70,J)
  tTau2=rep(trialNoise^2,J)
  tAlpha=makeAlpha(I,J)
  tTheta=makeTheta(I,tMu,tLambda,tEta,tDel2)
  tSigma=crossprod(t(tLambda))+diag(tDel2)
  return(list(
    "lambda"=tLambda,
    "eta"=tEta,
    "delta2"=tDel2,
    "mu"=tMu,
    "tau2"=tTau2,
    "alpha"=tAlpha,
    "theta"=tTheta,
    "Sigma"=tSigma))
}

makeData=function(truth,L){
  I=dim(truth$eta)[1]
  J=dim(truth$lambda)[1]
  dat=makeFlatDesign(I,J,L)
  dat=addObs(dat,truth)
  return(dat)
}


plotFacDecomp=function(truth,dat,samples){

  myCovPlot=function(cov,limits,main=""){
      cn <- paste("T",1:J,sep="")
      rownames(cov)=rep("",J)
      colnames(cov)=rep("",J)
      corrplot(cov,is.corr = F,method="color",
               cl.pos='n',col.lim=limits,main=main,
               mar=c(0,0,1,0))}
      

  myRowPlot=function(covs,labels=F,rowTitle="Test"){
    plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
    text(.6,.5,rowTitle,adj=.5,cex=1.2,srt=90)
    l=rep("",4)
    if (labels==T) l=c("Total Variance","Factor 1","Factor 2","Residual")
    myLim=c(0,max(covs))
    myCovPlot(covs[1,,],limits=myLim,main=l[1])
    plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
    text(.5,.5,"=",adj=.5,cex=2)
    myCovPlot(covs[2,,],limits=myLim,main=l[2])
    plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
    text(.5,.5,"+",adj=.5,cex=2)
    myCovPlot(covs[3,,],limits=myLim,main=l[3])
    plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
    text(.5,.5,"+",adj=.5,cex=2)
    myCovPlot(covs[4,,],limits=myLim,main=l[4])
    return(apply(covs, 1, function(mat) sum(diag(mat))))}
  
  J=dim(truth$lambda)[1]
  ssd=matrix(ncol=4,nrow=4)
  layout(matrix(nrow=4,ncol=8,1:32,byrow=T),widths=c(.3,1,.3,1,.3,1,.3,1))
  par(mar=c(1,1,1,1),mgp=c(1,1,0))
  covs=array(dim=c(4,J,J))
  covs[1,,]=truth$Sigma
  covs[2,,]=crossprod(t(truth$lambda[,1]))
  covs[3,,]=crossprod(t(truth$lambda[,2]))
  covs[4,,]=diag(truth$delta2)
  ssd[1,]=myRowPlot(covs,labels=T,rowTitle="Population")
  
  tTheta=data.frame(truth$theta)
  colnames(tTheta)=paste("T",1:8,sep="")
  lavaanMod <- 
    'efa("efa")*f1 + efa("efa")*f2 =~ T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8'
  fit <- cfa(lavaanMod,data=tTheta,rotation="varimax")
  loadings <- parameterEstimates(fit)
  lambda=matrix(ncol=2,loadings$est[1:16])
  resVar=loadings$est[17:24]
  #par(mfrow=c(1,4))
  covs[1,,]=cov(tTheta)
  covs[2,,]=crossprod(t(lambda[,1]))
  covs[3,,]=crossprod(t(lambda[,2]))
  covs[4,,]=diag(resVar)
  ssd[2,]=myRowPlot(covs,labels=F,rowTitle = "Individual")
  
  m=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
  effect=data.frame(m[,,2]-m[,,1])
  colnames(effect)=paste("T",1:J,sep="")
  fit <- cfa(lavaanMod,data=effect,rotation="varimax")
  loadings <- parameterEstimates(fit)
  lambda=matrix(ncol=2,loadings$est[1:16])
  resVar=loadings$est[17:24]
  
  covs[1,,]=cov(effect)
  covs[2,,]=crossprod(t(lambda[,1]))
  covs[3,,]=crossprod(t(lambda[,2]))
  covs[4,,]=diag(resVar)
  ssd[3,]=myRowPlot(covs,labels=F,rowTitle = "Trial Noise")


  aligned=align(samples)
  alignedPos=makePositive(aligned)
  order=alignToTruth(alignedPos$lambda,truth$lambda)
  pLambda=apply(alignedPos$lambda,2:3,mean)[,order]
  s2=apply(1/samples$pDelta2,2,mean)
  covs[1,,]=crossprod(t(pLambda))+diag(s2)
  covs[2,,]=crossprod(t(pLambda[,1]))
  covs[3,,]=crossprod(t(pLambda[,2]))
  covs[4,,]=diag(s2)
  ssd[4,]=myRowPlot(covs,labels=F,rowTitle = "Model")
  return(ssd)}

