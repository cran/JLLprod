\name{LL}
\alias{LL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Homothetic Production Function: Main Estimator }
\description{
This function implements the Lewbel & Linton's (2003) estimator. In general it estimates
the model Y=r(x,z)+e, imposing the following structure r(x,z)=E[Y|X=x,Z=z]=h[g(x,z)], and g(bx,bz)=b*g(x,z). The unknown function g is assumed
to be smooth and h is assumed to be a strictly monotonic smooth function.
}
\usage{
LL(xx, zz, yy, xxo, zzo, Vmin=NULL, Vmax = NULL,
   k, j, h = NULL, kernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Numerical: Nx1 vector.}
  \item{zz}{Numerical: Nx1 vector.}
  \item{yy}{Numerical: Nx1 vector.}
  \item{xxo}{Scalar: Normalization in xx direction.}
  \item{zzo}{Scalar: Normalization in zz direction. }
  \item{Vmin}{Scalar: Minimum value for Vk's, where Vk's are elements in the interval [Vmin,Vmax]. Default value is -1.}
  \item{Vmax}{Scalar: Maximum value for Vk's, where Vk's are elements in the interval [Vmin,Vmax]. Default value is 3.}
  \item{k}{Scalar: See Lewbel & Linton (2003). There is NO default, you must provide a number, i.e. 80}
  \item{j}{Scalar: See Lewbel & Linton (2003). There is NO default, you must provide a number, i.e. 100}
  \item{h}{Numerical: 2x1 vector of bandwidths, [hxx,hzz], used in the estimation. Default is the Silverman's rule of thumb in each direction.}
  \item{kernel}{Kernel function in all steps. Default is `gauss'.}
}
\details{
User may choose a variety of kernel functions. For example `uniform', `triangular', `quartic', `epanech', `triweight' or `gauss', see Yatchew (2003), pp 33. Another choice may be `order34', `order56' or `order78',
 which are third, fifth and seventh (gauss based) order kernel functions, see Pagan and Ullah (1999), pp 55.

Vmax should be chosen wisely. The user should make sure that Vmax*xxo belongs to the interior of the observed support of xx and similarly Vmax*zzo belongs to the interior
of the observed support of zz.
}
\value{
  \item{r}{N x 1 vector: Unrestricted Nonparametric r (see above) evaluated at data points, i.e. r(xxi,zzi).}
  \item{g}{N x 1 vector: Nonparametric component g (see above) evaluated at data points, i.e. g(xxi,zzi).}
  \item{h}{N x 1 vector: Nonparametric component h (see above) evaluated at data points, i.e. h[g(xxi,zzi)].}
  \item{hd}{N x 1 vector: Nonparametric first derivative of h (see above) evaluated at data points, i.e. h'[g(xxi,zzi)].}
}
\references{Lewbel, A., and Linton, O.B. (2003) Nonparametric Estimation of Homothetic and 
Homothetically Separable Functions. Unpublished manuscript.

Yatchew, A. (2003) Semiparametric Regression for the Applied Econometrician. Cambridge University Press.

Pagan, A. and Ullah, A. (1999) Nonparametric Econometrics. Cambridge Universtiy Press.
}
\author{ David Tomás Jacho-Chávez }
\section{Warning }{Simple and fast kernel regression is used in each step for computational time gain. However, it could take several minutes to complete for sample
sizes bigger than 300 observations.

k & j should be chosen accordingly with Vmin and Vmax. Try keeping them below 100.
}

\seealso{\code{\link{JLL}}, \code{\link{LLef}} , \code{\link{locpoly}}, \code{\link{Blocc}}}
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

#win.graph()
nf <- layout(matrix(c(1,2,1,2),2,2, byrow=TRUE),respect=TRUE)
plot(log(K)-log(L),log(LLb$g)-log(L),pch=3,xlab="ln(K/L)"
     ,ylab="ln(g(K/L,1))",main="Homogeneous Component g")
plot(log(LLb$g),LLb$h,xlab="ln(g)",pch=3,ylab="h(g)"
     ,main="Nonhomogeneous Component h",ylim=c(min(min(LLb$h)
     ,min(LLb$r)),max(max(LLb$h),max(LLb$r))))
points(log(LLb$g),LLb$r,type="p",pch=1,col="blue",lwd=2)
legend(-0.5,15.5,c("Nonparametric","Kernel Regression")
       ,merge=TRUE,lty=c(1,-1),pch=c(3,1),lwd=c(1,2)
       ,col=c("black","blue"),cex=0.95)
}
\keyword{smooth}
\keyword{regression}
