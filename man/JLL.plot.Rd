\name{JLL.plot}
\alias{JLL.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots of Nonparametric Elements of a Generalized Homothetic Production Function }
\description{
This function will return a self-explanatory graph of each component in the following production function model: E[Y|k=K/L,L=L]=r(K/L,L)=H[G(K/L)F(L)],
where Y is log of output, K is capital and L is labour.
}
\usage{
JLL.plot(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m }{ JLL object created by JLL function.}
}
\details{
This is a generic function which will allow the user to visualize each component in a Generalized Homothetic Production Function.
}
\value{
  A graph with 5 self-explanatory plots with labels.
}
\references{Jacho-Chávez, D.T., Lewbel, A., and Linton, O.B. (2005) Identification and Nonparametric Estimation
of a Transformed Additively Separable Model. Unpublished manuscript.}
\author{ David Tomás Jacho-Chávez }
\section{Warning }{The user should take a look at the reference below for a description of each plane in outcome window.} 

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
g<-seq(min(tlnklnl[,3]),max(tlnklnl[,3]),length=29)

#This part does the actual calculations & plots
m.JLL <- JLL(xx=tlnklnl[,2],zz=tlnklnl[,3],yy=tlnklnl[,1],p1=2,p2=2,
             pstar=0,h1=3.125,h2=5.25,ngrid=30,z0=g[25],r0=7.5)
JLL.plot(m.JLL)
}
\keyword{smooth}

