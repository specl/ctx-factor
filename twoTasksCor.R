require(abind)
source('facMod.R')
source('facExtra.R')

draw2TaskFigure <- function(I,fileName){
  fileNameMod=paste(fileName,"Mod",".Rds",sep= "")
  fileNameCovg=paste(fileName,"Covg",".Rds",sep= "")
  
  makeData=function(I,L=100){
    sub=rep(1:I,each=4*L)
    task=rep(rep(1:2,each=2*L),I)
    cond=rep(rep(1:2,each=L),I*2)
    tAlpha=runif(I,500,1000)
    tTheta=rmvnorm(I,c(70,70),matrix(ncol=2,c(1,.5,.5,1)*25^2))
    y=rnorm(I*L*4,tAlpha[sub]+(cond-1.5)*tTheta[cbind(sub,task)],175)
    dat=data.frame(sub,task,cond,y)
    return(list("dat"=dat,"tTheta"=tTheta))}

  conventional=function(I,L=100){
    dat=makeData(I)$dat
    ybar=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
    score=ybar[,,2]-ybar[,,1]
    vals=cor.test(score[,1],score[,2])
    out=c(vals$estimate,vals$conf.int)
    names(out) = c("est","lo","hi")
    return(out)}

  oneIval=function(Ivalue,R=100){
    I=rep(Ivalue,R)
    vals=data.frame(I,t(sapply(I,conventional)))
    o=order(vals$est)
    return(vals=vals[o,])}

  
  x=makeData(I)
  ybar=tapply(x$dat$y,list(x$dat$sub,x$dat$task,x$dat$cond),mean)
  score=ybar[,,2]-ybar[,,1]
  tc=cor.test(score[,1],score[,2])

  runChain=!file.exists(fileNameCovg)  
  if (runChain) {
    g=lapply(c(20,100,500),oneIval)
    saveRDS(g,fileNameCovg)}
  g=readRDS(fileNameCovg)
    
  
  prior=list(
    "mu.m"=70,
    "mu.s"=100,
    "S"=rep(30,2),
    "tau.s"=200,
    "alpha.m"=700,
    "alpha.s"=1000)
  if (runChain) {
    out=cov.run(x$dat,M=2000,prior=prior)
    saveRDS(out,fileNameMod)}
  out=readRDS(fileNameMod)
  Sigma=out$BUGSoutput$sims.list$Sigma
  M=dim(Sigma)[1]
  rho=1:M
  for (m in 1:M) rho[m]=cov2cor(solve(Sigma[m,,]))[1,2]

  range=range(score)
  par(mfrow=c(2,2),mar=c(4,4,2,2),mgp=c(2,1,0))
  plot(x$tTheta,ylab="True Score Task 2 (ms)",
       xlab="True Score Task 1 (ms)",
       pch=19,col='darkblue')
  mtext(side=3,adj=0,cex=1.1,"A")
  plot(score,ylab="Observed Score Task 2 (ms)",
       xlab="Observed Score Task 1 (ms)",
       pch=19,col='darkred',xlim=range,ylim=range)
  string1=paste("r = ",round(tc$estimate,2),sep="")
  string2=paste("[",paste(round(tc$conf.int,2),collapse=","),"]",collapse="")
  #text(150,0,string1,adj=1)
  #text(150,-25,string2,adj=1)
  mtext(side=3,adj=0,cex=1.1,"B")
  
  vals=abind(g,along=1)
  M=dim(vals)[1]
  plot(1:M,vals[,2],typ='n',ylim=c(-1,1),axes=F,
       ylab="Observed Correlation",xlab="")
  axis(2)
  lightCols=c('red','blue','green')[as.integer(as.factor(vals[,1]))]
  darkCols=paste("dark",lightCols,sep="")
  abline(lwd=2,h=.5)
  segments(1:M,vals[,3],1:M,vals[,4],col=lightCols)
  points(1:M,vals[,2],pch=19,cex=.5,col=darkCols)
  text(M*c(1/6,1/2,5/6),.9,
       paste("I=",unique(vals[,1]),sep=""),
       col=unique(darkCols))
  mtext(side=3,adj=0,cex=1.1,"C")
  
  hist(rho,col='lightblue',prob=T,main="",xlab="Correlation Coefficient",
       ylab="Density",xlim=c(-1.1,1.1),ylim=c(0,3.25))
  xp=c(-1.1,-1,-1,1,1,1.1)
  yp=.5*c(0,0,1,1,0,0)
  lines(xp,yp)
  arrows(tc$conf.int[1],3,tc$conf.int[2],3,code=3,
         angle=90,length = .1,col='darkred')
  points(tc$estimate,3,pch=19,col='darkred')
  abline(v=.5,lwd=2,col='darkblue')
  mtext(side=3,adj=0,cex=1.1,"D")
}