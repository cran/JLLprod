"JLL.plot" <-
function(m){
b <- rbind(cbind(1,1,2,2,3,3),cbind(1,1,2,2,3,3),cbind(4,4,4,5,5,5),cbind(4,4,4,5,5,5))
#win.graph()
layout(b)
GG <- m$Ghat-mean(m$Ghat); FF <- m$Fhat-mean(m$Fhat);
MM <- diag(m$Mhat); d<-cbind(MM,m$hhat); d<-d[order(MM),]
plot(m$x,GG,type="l",lty=1,xlab="ln(K/L)",ylab="",lwd=3,main="ln(G(K/L))")
plot(m$z,FF,type="l",lty=1,xlab="ln(L)",ylab="",lwd=3,main="ln(F(L))")
plot(d[,1],d[,2],type="l",lty=1,xlab="ln(M)",ylab="",lwd=3,main="H(M)")
persp(m$z,m$x,t(m$Rhat),axes=TRUE,lty=1,lwd=0.7,xlab="ln(L)",
      expand=0.5,ylab="ln(K/L)", zlab="",main="r(K/L,L)", font.lab=1, font.main=1,theta= 320, phi=17,cex.main=0.95)
persp(m$z,m$x,t(m$Mhat),axes=TRUE,lty=1,lwd=0.7,xlab="ln(L)",
      expand=0.5,ylab="ln(K/L)", zlab="",main="M(K/L,L)", font.lab=1, font.main=1,theta= 320, phi=17,cex.main=0.95)
      }

