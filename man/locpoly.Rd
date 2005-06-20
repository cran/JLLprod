\name{locpoly}
\alias{locpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Local Polynomial Regression}
\description{
This procedure performs Local Polynomial Regression up to order 3 for the model Y=m(x) + e, where m(x)=E[Y|X=x], where x is a kx1 vector
of covariates.
}
\usage{
locpoly(y, x, h = NULL, p = NULL, targmat = NULL, der = 0,
        nobmin = NULL, kernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y }{Numerical: Nx1 vector. Dependent variable.}
  \item{x }{Numerical: Nxk matrix. Covariates.}
  \item{h }{Numerical: h is a common scalar bandwidth for each variable, or a kx1 vector of bandwidths for each variable, 
or a Nxk vector of local bandwidths.}
  \item{p }{Scalar: scalar indicating the degree of the polynomial to fit locally.  This must be less than or equal to 3.
 p = 1 or 2 are recommended for estimating the function.}
  \item{targmat }{Numerical: mxk matrix of target points to get smoothed values at, or if targmat=0, then y will be smoothed at all data points.
 For large datasets, the use of targmat is highly recommended, as it drastically reduces computation time.}
  \item{der }{Scalar: It indicates the function of interest to estimate. if der=0, then actual function estimated. 
If der=1 then matrix of 1st derivatives is estimated. If der=2 then matrix of 2nd derivatives is estimated, corresponding to column 1, column2 etc. and then unique cross derivatives in the order x[1]x[2],..x[1]x[k],x[2]x[3],..x[2]x[k],..}
  \item{nobmin }{Scalar: It specifies minimum number of observations required to
 estimate local polynomial at each point.}
  \item{kernel }{Kernel function. Default is `gauss'.}
}
\details{
If in a region with sparse data, having nobmin set to a certain number will cause the bandwidth to be locally cranked up to capture more observations. Use with care.

User may also choose a variety of kernel functions. For example `uniform', `triangular', `quartic', `epanech', `triweight' or `gauss' when p[l]=1,2 or 3, and
l=1,2,star, see Yatchew (2003), pp 33. Another choice may be `order34', `order56' or `order78' when p[ll]=0 for ll=2,star, which are third, fifth and seventh 
(gauss based) order kernel functions, see Pagan and Ullah (1999), pp 55.
}
\value{
  \item{yhat }{N(m)x1 vector of smoothed values for function of interest given by der. if der=0, mhat estimated, if der=1 an estimate of the 1st derivatives is provided. if der=2 an estimate of the second derivatives is provided, etc.}
  \item{varhat }{Second component of estimated variance function.}
}
\references{Yatchew, A. (2003) Semiparametric Regression for the Applied Econometrician. Cambridge University Press.

Pagan, A. and Ullah, A. (1999) Nonparametric Econometrics. Cambridge Universtiy Press.}
\author{David Tomás Jacho-Chávez}
\note{This function is an adaptation of that provided by Cameron Rookley (1999), written for GAUSS.}

\section{Warning }{This procedures does a local polynomial fitting at each evaluation point, so for very bid datasets it is highly recommended that the user makes use of the argument targmat.
The smaller the number of evaluation points, the faster this procedure will be.}

\seealso{\code{\link{Blocc}}}
\examples{
#A very simple case
n <- 5000
x <- runif(n); e <- rnorm(n,sd=1)
G <- function(x){(1/2)*sin(2*pi*x)}
y <- G(x) + e
xgrid <- seq(0,1,length=30)
m <- locpoly(y=y,x=x,h=0.08,p=1,targmat=xgrid,der=0,nobmin=5,kernel="gauss")
plot(xgrid,G(xgrid),xlab="x",ylab="G",main="",col="black",type="l",lty=1)
lines(xgrid,m$yhat,col="red",type="l",lty=2)
}
\keyword{ smooth }
\keyword{ regression }
