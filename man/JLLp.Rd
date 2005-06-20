\name{JLLp}
\alias{JLLp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parametric Generalized Homothetic Production Function }
\description{
This function simply fits 2 parametric models, P2 and P3, as described in Jacho-Chávez, Lewbel and Linton (2005).
}
\usage{
JLLp(lnY,lnK,lnL,theta,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{lnY}{Ln of Y, output.}
    \item{lnK}{Ln of K, capital.}
    \item{lnL}{Ln of L, labour.}
    \item{theta}{A list of starting values of the form list(a,b0,b1,b2,g) if model=2 or a list of the form list(a,b0,b1,b2) if model=3.}
    \item{model}{Scalar, 2 or 3.}
}
\details{
This function is a simple call to `nls' to fit specific parametric models for production by NonLinear Least Squares.
}
\value{
  A `nls' object.
}
\references{Jacho-Chávez, D.T., Lewbel, A., and Linton, O.B. (2005) Identification and Nonparametric Estimation
of a Transformed Additively Separable Model. Unpublished manuscript.}
\author{ David Tomás Jacho-Chávez }
\section{Warning }{This function is very sensible to theta. It may fail most of the time. If this is the case, the user is adviced to
use their own call to `nls' or `optim' which is most likely to work.}

\seealso{\code{\link{JLL}}}
\examples{
library(JLLprod)
data(ecu)
##This part simply does some data sorting & trimming
xlnK <- ecu$lnk
xlnL <- ecu$lnl
xlnY <- ecu$lny
xqKL <- quantile(xlnK-xlnL,  probs=c(2.5,97.5)/100)
yx <- cbind(xlnY,xlnK,xlnL)
tlnklnl <- yx[((yx[,2]-yx[,3])>=xqKL[1]) & ((yx[,2]-yx[,3])<=xqKL[2]),]
tlnklnl[,2]<-tlnklnl[,2]-tlnklnl[,3]

bb<-list(b0=11,b1=1,b2=0,g=-0.15,a=0.4)
Y <- tlnklnl[,1]; K <- tlnklnl[,2]; L <- tlnklnl[,3]
pJLL<-JLLp(Y,K,L,theta=bb,model=2)
print(summary(pJLL))
}
\keyword{regression}

