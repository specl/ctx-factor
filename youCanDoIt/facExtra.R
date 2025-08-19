require(infinitefactor)
require(abind)
require(GPArotation)


makePositive=function(samples,I=0){
  out=samples
  dimMean=apply(samples$lambda,3,mean)
  D=length(dimMean)
  if (is.null(samples$eta)) samples$eta=matrix(0,nrow=I,ncol=D)
  for (d in 1:D) {
    out$lambda[,,d]=sign(dimMean[d])*samples$lambda[,,d]
    out$eta[,,d]=sign(dimMean[d])*samples$eta[,,d]}
  return(out)}

align=function(samples,I=0){
  out=samples
  M=dim(samples$lambda)[1]
  numFactors=dim(samples$lambda)[3]
  if (is.null(samples$eta)) samples$eta=matrix(0,nrow=I,ncol=numFactors)
  lambdaList=lapply(1:M,function(x) samples$lambda[x,,])
  etaList=lapply(1:M,function(x) samples$eta[x,,])
  aligned=jointRot(lambda = lambdaList, eta = etaList)
  out$lambda=aperm(abind(aligned$lambda, along=3),c(3,1,2))
  out$eta=aperm(abind(aligned$eta, along=3),c(3,1,2))
  return(out)
}

alignToTruth=function(lambda,tLambda){
  pLambda=apply(lambda,2:3,mean)
  a=cor(tLambda,pLambda)
  D=dim(a)[1]
  M=dim(a)[2]
  o=1:D
  for (i in 1:D)
    o[i]=which(a[i,]==max(a[i,]))
  return(c(o,(1:M)[-o]))}

alignToLargestVar=function(samples){
  v=apply(samples$lambda,1,function(x) diag(crossprod((x))))
  tot=apply(v,1,mean)
  return(order(tot,decreasing=T))
}

eigenMax=function(samples){
  M=dim(samples$lambda)[1]
  pLam=apply(samples$lambda,2:3,mean)
  rotSamp=orderSamp=samples
  rot=tandemI(pLam)$Th
  for (m in 1:M) {
    rotSamp$lambda[m,,]=samples$lambda[m,,]%*%rot
    rotSamp$eta[m,,]=samples$eta[m,,]%*%rot}
  o=alignToLargestVar(rotSamp)
  orderSamp$lambda=rotSamp$lambda[,,o] 
  orderSamp$eta=rotSamp$eta[,,o] 
  return(orderSamp)}


standardize=function(samples){
  M=dim(samples$lambda)[1]
  J=dim(samples$lambda)[2]
  sigma=array(dim=c(M,J))
  lambda=samples$lambda
  delta2=samples$pDelta2
  Sigma=rho=array(dim=c(M,J,J))
  for (m in 1:M){
    Sigma[m,,]=crossprod(t(samples$lambda[m,,]))+diag(1/samples$pDelta2[m,])
    rho[m,,]=cov2cor(Sigma[m,,])
    sigma[m,]=sqrt(diag(Sigma[m,,]))
    lambda[m,,]=samples$lambda[m,,]/sigma[m,]
    delta2[m,]=(1/samples$pDelta2[m,])/sigma[m,]^2}
  samples$lambda=lambda
  samples$delta2=delta2
  samples$Sigma=Sigma
  samples$rho=rho
  return(samples)}
