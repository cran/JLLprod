\name{ecu}
\alias{ecu}
\docType{data}
\title{Production Data: Ecuador (2002)}
\description{
Production data for the Petroleum, Chemical & Plastics industry in Ecuador for the year 2002.
}
\usage{data(ecu)}
\format{
  A data frame with 406 observations on the following 4 variables.
  \describe{
    \item{lny}{Ln of Y, output in thousands of current USD.}
    \item{lnk}{Ln of K, capital in thousands of current USD.}
    \item{lnl}{Ln of L, the average number of employees.}
    \item{ind}{    2: Petroleum, Chemical & Plastics}
  }
}
\details{
We would like to thank The Central Bank of Ecuador for providing us with this data set, specially to Miguel Acosta.
}
\source{
INEC (National Institute of Statistic & Census of Ecuador)
}
\references{
Raw data prepared by: Sebastián Zambrano and Diego Benítez (dbenitez@uio.bce.fin.ec)
}
\examples{
data(ecu)
#some plots
#win.graph()
layout(matrix(c(1,2,1,2),2,2, byrow=TRUE),respect=TRUE)
plot(ecu$lnk-ecu$lnl,ecu$lny,xlab="ln(K/L)",ylab="ln(y)",main="")
plot(ecu$lnl,ecu$lny,xlab="ln(L)",ylab="ln(y)",main="")
par(oma=c(1,1,1.85,1),new=TRUE,font=2,cex=1)
mtext(outer=TRUE,"Petroleum, Chemical & Plastics - 2002",side=3,cex=1)
}
\keyword{datasets}
