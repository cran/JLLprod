\name{Blocc}
\alias{Blocc}

\title{Fast Local Constant Regression for Bivariate Data}
\description{
This procedure is a fast implementation of the Nadaraya-Watson estimator for conditional mean function such as
r(x,z)=E[y|X=x,Z=z].
}
\usage{
Blocc(xx, zz, yy, kernel = NULL, ev = NULL, h = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Numerical: Nx1 vector. XX direction.}
  \item{zz}{Numerical: Nx1 vector. ZZ direction.}
  \item{yy}{Numerical: Nx1 vector. Dependent variable.}
  \item{kernel}{Kernel function. Default is `gauss'.}
  \item{ev}{Numerical: Mx2 matrix of evaluation points to get smoothed values at. The use of ev is highly recommended, as it drastically reduces computation time.
 Default is a 40x2 matrix covering the entire observed support.}
  \item{h}{Numerical: 2x1 vector of bandwidths, [hxx,hzz], used in the estimation. Default is the Silverman's rule of thumb in each direction.}
}
\details{User may also choose a variety of kernel functions. For example `uniform', `triangular', `quartic', `epanech', `triweight' or `gauss', see Yatchew (2003), pp 33. Another choice may be `order34',
 `order56' or `order78', which are third, fifth and seventh 
(gauss based) order kernel functions, see Pagan and Ullah (1999), pp 55.
}
\value{
  \item{r }{Numerical: MxM (i,j) matrix of nonparametric estimates or r(xxe[i],zze[j]).}
  \item{xxe }{Numerical: Mx1 vector of evaluation points in the xx direction.}
  \item{zze }{Numerical: Mx1 vector of evaluation points in the zz direction.}
}
\references{Yatchew, A. (2003) Semiparametric Regression for the Applied Econometrician. Cambridge University Press.

Pagan, A. and Ullah, A. (1999) Nonparametric Econometrics. Cambridge Universtiy Press.}
\author{David Tomás Jacho-Chávez}
\note{This function may fail for very big values of N or M, because it uses matrices and NO loops.}

\seealso{\code{\link{locpoly}}}
\examples{
#A very simple case
n <- 1000
x <- runif(n); z <- runif(n); e <- rnorm(n,sd=0.2)
G <- function(x){(1/2)*sin(2*pi*x)}
F <- function(z){-1/3+2*z-2*(z^2)}
y <- G(x)+F(z) + e
xgrid <- seq(0,1,length=30); zgrid <- seq(0,1,length=30)

m <- Blocc(xx=x,zz=z,yy=y,ev=cbind(xgrid,zgrid))
GF <- matrix(G(m$xxe),nr=30,nc=30,byrow=FALSE)+matrix(F(m$zze)
             ,nr=30,nc=30,byrow=TRUE)

#win.graph()
layout(matrix(c(1,2),nr=1,nc=2,byrow=TRUE))
persp(x=m$zze,y=m$xxe,z=t(GF),theta= 320, phi=17,xlab="z"
      ,ylab="x",zlab="",main="True G(x)+F(z)")
persp(x=m$zze,y=m$xxe,z=t(m$r),theta= 320, phi=17,xlab="z"
      ,ylab="x",zlab="",main="Estimated G(x)+F(z)")
}
\keyword{ smooth }
\keyword{ regression }
