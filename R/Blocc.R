"Blocc" <-
function(xx,zz,yy,kernel=NULL,ev=NULL,h=NULL){
if (is.null(ev)){xxe<- seq(min(xx),max(xx),length=40); zze<- seq(min(zz),max(zz),length=40)}
else {xxe <- ev[,1]; zze <- ev[,2]}

if (length(xx)!=length(zz)){stop("JLLprod:lengths of x and z are NOT the same, I won't continue!")}
else

if (is.null(h)){hx=1.06*sd(xx)*(length(xx)^-0.2); hz=1.06*sd(zz)*(length(zz)^-0.2)}
else if (length(h)==1){hx=h; hz=h}
else if (length(h)==2){hx=h[1];hz=h[2]}
else {stop("JLLprod: Vector of bandwidths is NOT either scalar or a 2x1 vector, I won't continue!")}

if (is.null(kernel)) k <- function(u){dnorm(u,mean=0,sd=1)}       
else if (kernel=="uniform") k <- function(u){(1/2)*(-1 <= u)*(u <= 1)}
else if (kernel=="triangular") k <- function(u){(1-abs(u))*(-1 <= u)*(u <= 1)}
else if (kernel=="quartic") k <- function(u){(15/16)*((1-u^2)^2)*(-1 <= u)*(u <= 1)}
else if (kernel=="epanech") k <- function(u){(3/4)*(1-u^2)*(-1 <= u)*(u <= 1)}
else if (kernel=="triweight") k <- function(u){(35/32)*((1-u^2)^3)*(-1 <= u)*(u <= 1)}
else if (kernel=="gauss") k <- function(u){dnorm(u,mean=0,sd=1)}
else if (kernel=="order34") k <- function(u){(1/2)*(3-u^2)*dnorm(u)}
else if (kernel=="order56") k <- function(u){(1/8)*(15-10*u^2+u^4)*dnorm(u)}
else if (kernel=="order78") k <- function(u){(1/48)*(105-105*u^2+21*u^4-u^6)*dnorm(u)}
else {}

   m1 = (as.matrix(xxe)%*%matrix(1,nr=1,nc=length(xx)) - matrix(1,nr=length(xxe),nc=1)%*%t(as.matrix(xx)))/hx
   K1 = k(m1)/hx
   m2 = (as.matrix(zze)%*%matrix(1,nr=1,nc=length(zz)) - matrix(1,nr=length(zze),nc=1)%*%t(as.matrix(zz)))/hz
   K2 = k(m2)/hz
   S1 = K1%*%t(K2)/length(xx); T1 = (K1*matrix(yy,nr=length(xxe),nc=length(xx),byrow=TRUE))%*%t(K2)/length(xx);
   return(list(r=T1/S1,xxe=xxe,zze=zze))
}

