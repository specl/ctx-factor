br=function(xc=.5,yc=.5,len=.25,chevron=.1,dir=1,add=F){
  if (add==F) plot(NA,NA,xlim = 0:1,ylim =0:1,typ='n',axes=F)
  segments(xc-len-dir*chevron,yc-chevron*c(1,-1),xc-len,yc,col='black',lwd=2)
  segments(xc+len-dir*chevron,yc-chevron*c(1,-1),xc+len,yc,col='black',lwd=2)
  segments(xc+dir*chevron,yc-chevron*c(1,-1),xc,yc,col='black',lwd=2)
  segments(xc-len,yc,xc,yc,col='red',lwd=2)
  segments(xc+len,yc,xc,yc,col='blue',lwd=2)}

zol=function(){
  y=.15*seq(-2,2,1)+yc
  x=seq(.1,.9,.1)
  par(mar=c(0,0,0,0))
  par(lwd=2)
  plot(NA,NA,xlim = 0:1,ylim =0:1,typ='n',axes=F)
  segments(.1,y,.9,y)
  for (i in 1:5) {
    mult=c(-1,1)[i%%2+1]
    segments(x-mult*.05,
             y[i]-1*.05,
             x+mult*.05,
             y[i]+1*.05)}
}

