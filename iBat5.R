loadiBat5=function(fn){
  dat=read.csv(fn)
  dat$taskN=as.integer(as.factor(dat$task))
  dat$tv=(dat$taskN-1)*2+dat$version
  labelsTV=paste(rep(levels(as.factor(dat$task)),each=2),rep(1:2,5),sep="")
  dat$sub=as.integer(as.factor(dat$sub))
  I=max(dat$sub)
  v1=tapply(dat$y,list(dat$sub,dat$tv),var)
  divisor=sqrt(apply(v1,2,mean))
  dat$yc=dat$y/divisor[dat$tv]
  out=data.frame(dat$sub,dat$tv,dat$yc)
  colnames(out)=c('sub','task','y')
  return(out)}


free=matrix(1,nrow=10,ncol=7)
free[,1:5]=0
free[9:10,1]=c(2,1)
free[7:8,2]=c(2,1)
free[5:6,3]=c(2,1)
free[3:4,4]=c(2,1)
free[1:2,5]=c(2,1)



fitiBat5=function(dat,prior){

  source('facMod.R')
  source('facExtra.R')
  out=facH0cfa.run(dat,prior=prior,free=free,M=5000)
  samples=out$BUGSoutput$sims.list
  samplesFree=list(
    "lambda"=samples$lambda[,,6:7],
    "eta"=samples$eta[,,6:7])
  alignedFree=eigenMax(align(samplesFree))
  myLambda=abind(samples$lambda[,,1:5],alignedFree$lambda,along=3)
  myEta=abind(samples$eta[,,1:5],alignedFree$eta,along=3)
  myPars=list(
    "lambda"=myLambda,
    "eta"=myEta)
  samples2=makePositive(myPars)
  samples$lambda=samples2$lambda
  samples$eta=samples2$eta
  return(samples)
}


myCorPlot=function(cormat,limits=c(-1,1),main=""){
  J=dim(cormat)[1]
  cn <- paste("T",1:J,sep="")
  rownames(cormat)=rep("",J)
  colnames(cormat)=rep("",J)
  corrplot(cormat,method="color",
           cl.pos='n',col.lim=limits,main=main,
           mar=c(0,0,1,0))}

iBat5decomp=function(std){

  M=dim(std$lambda)[1]
  J=dim(std$lambda)[2]
  cm=array(dim=c(5,J,J))
  prop=array(dim=c(5,M,J,J))
  prop[1,,,]=std$rho
  efac=cfac1=cfac2=resid=array(dim=c(M,J,J))
  for (m in 1:M){
    prop[2,m,,]=crossprod(t(std$lambda[m,,1:5]))
    prop[3,m,,]=crossprod(t(std$lambda[m,,6]))
    prop[4,m,,]=crossprod(t(std$lambda[m,,7]))
    prop[5,m,,]=diag(std$delta2[m,])
  }
  cm=apply(prop,c(1,3,4),mean)
  layout(matrix(nrow=1,ncol=9,1:9,byrow=T),widths=c(1,.3,1,.3,1,.3,1,.3,1))
  par(mar=c(1,1,1,1))
  myCorPlot(cm[1,,],main="Correlation")
  plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
  text(.5,.5,"=",adj=.5,cex=2)
  myCorPlot(cm[2,,],main="Nuisance")
  plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
  text(.5,.5,"+",adj=.5,cex=2)
  myCorPlot(cm[3,,],main="Factor 6")
  plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
  text(.5,.5,"+",adj=.5,cex=2)
  myCorPlot(cm[4,,],main="Factor 7")
  plot(0:1,0:1,ty='n',axes=F,xlab="",ylab="")
  text(.5,.5,"+",adj=.5,cex=2)
  myCorPlot(cm[5,,],main="Residual")
  return(prop)
}

br=function(xc=.5,yc=.5,len=.25,chevron=.1,dir=1,add=F){
  if (add==F) plot(NA,NA,xlim = 0:1,ylim =0:1,typ='n',axes=F,xlab="",ylab="")
  segments(xc-len-dir*chevron,yc-chevron*c(1,-1),xc-len,yc,col='black',lwd=2)
  segments(xc+len-dir*chevron,yc-chevron*c(1,-1),xc+len,yc,col='black',lwd=2)
  segments(xc+dir*chevron,yc-chevron*c(1,-1),xc,yc,col='black',lwd=2)
  segments(xc-len,yc,xc,yc,col='red',lwd=2)
  segments(xc+len,yc,xc,yc,col='blue',lwd=2)}

zol=function(xc=.5,yc=.5){
  y=.15*seq(-2,2,1)+yc
  x=seq(.1,.9,.1)
  par(lwd=2)
  plot(NA,NA,xlim = 0:1,ylim =0:1,typ='n',axes=F,xlab="",ylab="")
  segments(.1,y,.9,y)
  for (i in 1:5) {
    mult=c(-1,1)[i%%2+1]
    segments(x-mult*.05,
             y[i]-1*.05,
             x+mult*.05,
             y[i]+1*.05)}
}