"LLef" <-
function(xx,zz,yy,LLob,h0=NULL,kernel0=NULL,kernel=NULL){
if (length(xx)!=length(zz)){stop("JLLProd: lengths of x and z are NOT the same, I won't continue!")
} else {}

if (is.null(kernel0)) ker <- function(u){dnorm(u,mean=0,sd=1)
} else if (kernel0=="uniform") ker <- function(u){(1/2)*(-1 <= u)*(u <= 1)
} else if (kernel0=="triangular") ker <- function(u){(1-abs(u))*(-1 <= u)*(u <= 1)
} else if (kernel0=="quartic") ker <- function(u){(15/16)*((1-u^2)^2)*(-1 <= u)*(u <= 1)
} else if (kernel0=="epanech") ker <- function(u){(3/4)*(1-u^2)*(-1 <= u)*(u <= 1)
} else if (kernel0=="triweight") ker <- function(u){(35/32)*((1-u^2)^3)*(-1 <= u)*(u <= 1)
} else if (kernel0=="gauss") ker <- function(u){dnorm(u,mean=0,sd=1)
} else {}

if (is.null(kernel)){ kernel<-"gauss"
} else {kernel<-kernel}

r <- sqrt(xx^2+zz^2)
t <- atan(xx/zz)
g <- LLob$g; h <- LLob$h; hd <- LLob$hd; rhat <- LLob$r
G <- g/r

if (is.null(h0)){st <- 1.06*length(G)^-0.2
} else {st<-h0}

m = (as.matrix(t)%*%matrix(1,nr=1,nc=length(t)) - matrix(1,nr=length(t),nc=1)%*%t(as.matrix(t)))/st
Ktheta = ker(m)/(length(t)*st)
seta <- (rhat-h)*hd*r - mean((rhat-h)*hd*r*G)

Gef <- G - (Ktheta%*%seta)/(Ktheta%*%((hd^2)*(r^2)))
gef <- r*Gef
H <- locpoly(y=yy,x=gef,h=1.06*sd(gef)*length(gef)^-0.2,p=1,targmat=gef,der=0,kernel=kernel);
Hd <- locpoly(y=yy,x=gef,h=1.06*sd(gef)*length(gef)^-0.2,p=1,targmat=gef,der=1,kernel=kernel);
hef <- H$yhat; hdef <- Hd$yhat
return(list(gef=gef,hef=hef,hdef=hdef))
}

