\name{LLef}
\alias{LLef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Homothetic Production Function: Most Efficient Estimator }
\description{
This function implements the most efficient version of Lewbel & Linton's (2005) estimator. In general it estimates
the model Y=r(x,z)+e, imposing the following structure r(x,z)=E[Y|X=x,Z=z]=h[g(x,z)], and g(bx,bz)=b*g(x,z). The unknown function g is assumed
to be smooth and h is assumed to be a strictly monotonic smooth function. 
}
\usage{
LLef(xx, zz, yy, LLob, h0 = NULL, kernel0=NULL , kernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Numerical: Nx1 vector.}
  \item{zz}{Numerical: Nx1 vector.}
  \item{yy}{Numerical: Nx1 vector.}
  \item{LLob}{LL object}
  \item{h0}{Numerical: Bandwidth for smoothing. Default is the Silverman's rule of thumb. }
  \item{kernel0}{Kernel function for smoothing. Default is `gauss'.}
  \item{kernel}{Kernel function for all steps. Default is `gauss'.}
}
\details{
User may choose a variety of kernel functions. For example `uniform', `triangular', `quartic', `epanech', `triweight' or `gauss', see Yatchew (2003), pp 33. Another choice may be `order34', `order56' or `order78',
 which are third, fifth and seventh (gauss based) order kernel functions, see Pagan and Ullah (1999), pp 55.
}
\value{
  \item{gef}{N x 1 vector: Efficient Nonparametric component g (see above) evaluated at data points, i.e. g(xxi,zzi).}
  \item{hef}{N x 1 vector: Efficient Nonparametric component h (see above) evaluated at data points, i.e. h[g(xxi,zzi)].}
  \item{hdef}{N x 1 vector: Efficient Nonparametric first derivative of h (see above) evaluated at data points, i.e. h'[g(xxi,zzi)].}
}
\references{Lewbel, A., and Linton, O.B. (2005) Nonparametric Matching and Efficient Estimation of Homothetically Separable Functions. Unpublished manuscript.

Yatchew, A. (2003) Semiparametric Regression for the Applied Econometrician. Cambridge University Press.

Pagan, A. and Ullah, A. (1999) Nonparametric Econometrics. Cambridge Universtiy Press.
}
\author{ David Tomás Jacho-Chávez }


\seealso{\code{\link{JLL}}, \code{\link{LL}} , \code{\link{locpoly}}, \code{\link{Blocc}}}
\examples{

data(ecu)
##This part simply does some data sorting & trimming
xlnK <- ecu$lnk
xlnL <- ecu$lnl
xlnY <- ecu$lny
xqKL <- quantile(exp(xlnK)/exp(xlnL),  probs=c(2.5,97.5)/100)
yx <- cbind(xlnY,xlnK,xlnL)
tlnklnl <- yx[((exp(yx[,2])/exp(yx[,3]))>=xqKL[1]) 
              & ((exp(yx[,2])/exp(yx[,3]))<=xqKL[2]),]
Y<-tlnklnl[,1]
K<-exp(tlnklnl[,2])/median(exp(tlnklnl[,2]))
L<-exp(tlnklnl[,3])/median(exp(tlnklnl[,3]))

LLb<-LL(xx=K,zz=L,yy=Y,xxo=median(K),zzo=median(L),k=80,j=100)
LLbef <- LLef(xx=K,zz=L,yy=Y,h0=1,LLob=LLb)

#win.graph()
nf <- layout(matrix(c(1,2,1,2),2,2, byrow=TRUE),respect=TRUE)
plot(log(K)-log(L),log(LLbef$gef)-log(L),pch=3,xlab="ln(K/L)"
     ,ylab="ln(g(K/L,1))",main="Homogeneous Component g")
plot(log(LLbef$gef),LLbef$hef,xlab="ln(g)",pch=3,ylab="h(g)"
     ,main="Nonhomogeneous Component h",ylim=c(min(min(LLbef$hef)
     ,min(LLb$r)),max(max(LLbef$hef),max(LLb$r))))
points(log(LLbef$gef),LLb$r,type="p",pch=1,col="blue",lwd=2)
legend(-0.5,15.6,c("Nonparametric","Kernel Regression")
       ,merge=TRUE,lty=c(1,-1),pch=c(3,1),lwd=c(1,2)
       ,col=c("black","blue"),cex=0.95)
}
\keyword{smooth}
\keyword{regression}
