\name{JLL}
\alias{JLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Homothetic Production Function: Main Estimator}
\description{
This function implements the Jacho-Chávez, Lewbel & Linton's (2005) estimator. In general it estimates
the model Y=r(x,z)+e, imposing the following structure r(x,z)=E[Y|X=x,Z=z]=H[M(x,z)], and M(x,z)=G(x)+F(z). The unknown functions G and F are assumed
to be smooth and H is assumed to be strictly monotonic.
}
\usage{
JLL(xx, zz, yy, z0, r0, ngrid = NULL, h1 = NULL, h2 = NULL,
    hstar = NULL, k1 = NULL, k2 = NULL, kstar = NULL, p1 = NULL,
    p2 = NULL, pstar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{Numerical: Nx1 vector.}
  \item{zz}{Numerical: Nx1 vector.}
  \item{yy}{Numerical: Nx1 vector.}
  \item{z0}{Scalar: Scale and sign normalization, i.e. F'(z0)=f(z0)=1.}
  \item{r0}{Scalar: Location normalization, i.e. r0=H[0].}
  \item{ngrid}{scalar: number of equally spaced grid points in each direction at which the model will be estimated.}
  \item{h1}{Numerical: 2x1 vector in the form c(hxx,hzz). Scalar: same bandwidth in each direction. First step bandwidth. Default is the Silverman's rule of thumb in each direction.}
  \item{h2}{Numerical: 2x1 vector in the form c(hr,hzz). Scalar: same bandwidth in each direction. Second step bandwidth. Default is the Silverman's rule of thumb in each direction.}
  \item{hstar}{Scalar: bandwidth used in the estimation of H.}
  \item{k1}{Kernel function in step 1. Default is `gauss'.}
  \item{k2}{Kernel function in step 2. Default is `gauss'.}
  \item{kstar}{Kernel function in step 3. Default is `gauss'.}
  \item{p1}{Order of polynomial to be used in step 1. From 0 up to 3.}
  \item{p2}{Order of polynomial to be used in step 2. From 0 up to 3.}
  \item{pstar}{Order of polynomial to be used in step 3. From 0 up to 3.}
}
\details{
User may choose a variety of kernel functions. For example `uniform', `triangular', `quartic', `epanech', `triweight' or `gauss' when p[l]=1,2 or 3, and
l=1,2,star, see Yatchew (2003), pp 33. Another choice may be `order34', `order56' or `order78' when p[ll]=0 for ll=2,star, which are third, fifth and seventh 
(gauss based) order kernel functions, see Pagan and Ullah (1999), pp 55.
}
\value{
  \item{Rhat}{ngrid x ngrid matrix: Unrestricted nonparametric fit at grid points.}
  \item{Mhat}{ngrid x ngrid matrix: Nonparametric component M (see above) evaluated at grid points.}
  \item{Ghat}{ngrid x 1 vector: Nonparametric component G (see above) evaluated at grid points.}
  \item{Fhat}{ngrid x 1 vector: Nonparametric component F (see above) evaluated at grid points.}
  \item{x}{ngrid x 1 vector: grid points in the xx direction.}
  \item{z}{ngrid x 1 vector: grid points in the zz direction.}
  \item{hhat}{ngrid x 1 vector: Nonparametric component H (see above) evaluated at cbind(x,z).}
}
\references{Jacho-Chávez, D.T., Lewbel, A., and Linton, O.B. (2005) Identification and Nonparametric Estimation
of a Transformed Additively Separable Model. Unpublished manuscript.

Yatchew, A. (2003) Semiparametric Regression for the Applied Econometrician. Cambridge University Press.

Pagan, A. and Ullah, A. (1999) Nonparametric Econometrics. Cambridge Universtiy Press.
}
\author{ David Tomás Jacho-Chávez }
\section{Warning }{Results seem to be more sensitive to the choice of z0 to the choice of r0. Try choosing values of z0 in the interior of the support of zz.
}

\seealso{\code{\link{JLL.plot}}, \code{\link{LL}}, \code{\link{LLef}} , \code{\link{locpoly}}}
\examples{
#A simple simulation
d<- 1;  n=400
xx <- runif(n^d); zz <- runif(n^d);
x <- cbind(xx,zz)
e <- rnorm(n^d,sd=1)
G <- function(x){(1/2)*sin(2*pi*x)}
F <- function(z){-2*(z^2)+2*z-1/3}
M <- function(x,z){G(x)+F(z)}
H <- function(m){1+(16/7)*m}
z0<-1/4; r0<-H(0);
y <- H(M(x[,1],x[,2])) + e
m <- JLL(xx=x[,1],zz=x[,2],yy=y,p1=3,p2=1,pstar=1,h1=0.15,h2=0.7,ngrid=10,
         z0=z0,r0=r0)

#win.graph()
layout(matrix(c(1,2,1,2),2,2, byrow=TRUE),respect=TRUE)
plot(m$x,m$Ghat,type="l",lty=1)
lines(m$x,G(m$x),type="l",lty=2,col="blue")
plot(m$z,m$Fhat,type="l",lty=1)
lines(m$z,F(m$z),type="l",lty=2,col="blue")
}
\keyword{smooth}
\keyword{regression}
