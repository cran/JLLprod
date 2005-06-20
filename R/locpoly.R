"locpoly" <-
function(y,x,h=NULL,p=NULL,targmat=NULL,der=0,nobmin=NULL,kernel=NULL){
if (!is.matrix(y)){y<-as.matrix(y);n = nrow(y)}else{y<-y;n = nrow(y)}
if (ncol(y)!=1){stop("LOCPOLY: y must be a Nx1 vector only")}else{}
if (!is.matrix(x)){x<-as.matrix(x);k = ncol(x);}else{x<-x;k = ncol(x);}
if (n!=nrow(x)){stop("LOCPOLY: y must have the same number of rows")}else{}
if (is.null(targmat)){targmat=x; m=n}else{targmat=as.matrix(targmat);m=nrow(as.matrix(targmat))}
if (is.null(h)){stop("LOCPOLY: I'm not that clever, so please specify a bandwidth")}else{}
if (is.null(p)){stop("LOCPOLY: I'm not that clever, so please specify a polynomial degree")}else{}
if (is.vector(h)){h <- as.matrix(h)}else{h <- as.matrix(h)}

if ((nrow(h)==1)&(nrow(h)==1)){h=matrix(h,nr=m,nc=k);
}else if ((nrow(h)==1)&(ncol(h)==k)){h=matrix(1,nr=m,nc=k)%*%h;
}else if ((nrow(h)==k)&(ncol(h)==1)){h=matrix(1,nr=m,nc=1)%*%t(h);
}else if ((nrow(h)==m)&(ncol(h)==ncol(targmat))){
}else if ((nrow(h)==m)&(ncol(h)==1)){h=h%*%matrix(1,nr=1,nc=k)
}else{stop("LOCPOLY: nrow of bandwidth vector inconsistent with ncol of x")}

if (is.null(nobmin)){nobmin=10
}else {nobmin=nobmin}

if (is.null(kernel)){kern <- function(u){dnorm(u,mean=0,sd=1)}
}else if (kernel=="uniform"){kern <- function(u){(1/2)*(-1 <= u)*(u <= 1)}
}else if (kernel=="triangular"){kern <- function(u){(1-abs(u))*(-1 <= u)*(u <= 1)}
}else if (kernel=="quartic"){kern <- function(u){(15/16)*((1-u^2)^2)*(-1 <= u)*(u <= 1)}
}else if (kernel=="epanech"){kern <- function(u){(3/4)*(1-u^2)*(-1 <= u)*(u <= 1)}
}else if (kernel=="triweight"){kern <- function(u){(35/32)*((1-u^2)^3)*(-1 <= u)*(u <= 1)}
}else if (kernel=="gauss"){kern <- function(u){dnorm(u,mean=0,sd=1)}
}else if (kernel=="order34"){kern <- function(u){(1/2)*(3-u^2)*dnorm(u)}
}else if (kernel=="order56"){kern <- function(u){(1/8)*(15-10*u^2+u^4)*dnorm(u)}
}else if (kernel=="order78"){kern <- function(u){(1/48)*(105-105*u^2+21*u^4-u^6)*dnorm(u)}
}else {}


if (is.null(kernel)){kerncut=4
}else if (kernel=="gauss"){kerncut=4
}else{kerncut=1}

if (p>3){stop("LOCPOLY: I cannot deal with p>3")}else{}

if (der==0){
   yhat=matrix(NA,nr=m,nc=1);
   if (p>=2){b2=as.matrix(seq(from=(k+2),to=(k+1+0.5*(k^2+k))))}else{}
}else if (der==1){
   if (p<1) {print("LOCPOLY: Resetting polynomial to 1, to estimate 1st derivitive");p=1}else{}
   if (k==1){b1=as.matrix(2)}else{b1=as.matrix(seq(2,length=k))}
   yhat=matrix(NA,nr=m,nc=k)
}else if (der==2){
   if (p<2){print("LOCPOLY: Resetting polynomial to 2, to estimate 2nd derivitive");p=2}else{}
   if (k==1){b2=as.matrix(3)}else{b2=as.matrix(seq(from=(k+2),to=(k+1+0.5*(k^2+k))));}
   yhat=matrix(NA,nr=m,nc=nrow(b2));
}else {stop("LOCPOLY: I cannot deal with p>3 or its derivatives")}

varhat=yhat; heff=matrix(NA,nr=m,nc=k);

if ((p==0)&(der==0)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        yhat[i]=beta[1,];
        varhat[i,]=varmat[1,1];
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else if ((p==1)&(der==0)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        yhat[i]=beta[1,];
        varhat[i,]=varmat[1,1];
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else if ((p==2)&(der==0)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=cbind(ind,xdiff,(xdiff^2));
      pc=1;
      while (pc <= k){sc=pc+1; 
                     while (sc <= k){ind=cbind(ind,xdiff[,pc]*xdiff[,sc]); sc=sc+1;}
                     pc=pc+1
                     }
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        yhat[i]=beta[1,];
        varhat[i,]=varmat[1,1];
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else if ((p==1)&(der==1)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=cbind(ind,xdiff);
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        if (k==1){varhat[i,]=varmat[b1,b1];
                  yhat[i]=beta[b1,];
                 }else{yhat[i,]=t(beta[b1,]);
                       varhat[i,]=t(diag(varmat[b1,b1]));}
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else if ((p==2)&(der==1)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=cbind(ind,xdiff,(xdiff^2));
      pc=1;
      while (pc <= k){sc=pc+1; 
                     while (sc <= k){ind=cbind(ind,xdiff[,pc]*xdiff[,sc]); sc=sc+1;}
                     pc=pc+1
                     }
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        if (k==1){varhat[i,]=varmat[b1,b1];
                  yhat[i]=beta[b1,];
                 }else{yhat[i,]=t(beta[b1,]);
                       varhat[i,]=t(diag(varmat[b1,b1]));}
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else if ((p==2)&(der==2)){
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      ind=cbind(ind,xdiff,(xdiff^2));
      pc=1;
      while (pc <= k){sc=pc+1; 
                     while (sc <= k){ind=cbind(ind,xdiff[,pc]*xdiff[,sc]); sc=sc+1;}
                     pc=pc+1
                      }
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        if (k==1){yhat[i,]=2*t(beta[b2,]);
                 varhat[i,]=4*varmat[b2,b2];
                 }else{yhat[i,]=t(beta[b2]);
                      varhat[i,]=t((diag(varmat[b2,b2])))}
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }
}
return(list(yhat=yhat,varhat=varhat));
}else{
for (i in 1:m){
  curh=h[i,];
  cont=1;
  while (cont==1){
   kernarg=(x-matrix(targmat[i,],nr=nrow(x),nc=k,byrow=TRUE))*matrix((curh)^-1,nr=nrow(x),nc=k,byrow=TRUE);
   mkernarg=apply(abs(kernarg),1,max);
   keepind=which((mkernarg > -kerncut)&(mkernarg<kerncut),arr.ind=TRUE);
   rk=length(keepind);
   if (rk >= nobmin){
      xdiff=x[keepind,]-matrix(targmat[i,],nr=rk,nc=k,byrow=TRUE);
      ki=apply(kern(xdiff*matrix((curh)^-1,nr=rk,nc=ncol(h),byrow=TRUE)),1,prod)
      w=sqrt(ki);
      ind=matrix(1,nr=rk,nc=1);
      if (p>=1){ind=cbind(ind,xdiff)}else{}
      if (p>=2){ind=cbind(ind,(xdiff)^2);
                pc=1;
                while (pc <= k){sc=pc+1;
                                while (sc <= k){ind=cbind(ind,xdiff[,pc]*xdiff[,sc]); sc=sc+1}
                                pc=pc+1
                                }
               }else{}
      if (p>=3){ind=cbind(ind,(xdiff)^3);
                pc=1;
                while (pc <= k){sc=1;
                                while (sc <= k){if (sc != pc){
                                                ind=cbind(ind,(xdiff[,pc]^2)*xdiff[,sc]); sc=sc+1
                                                }else{sc=sc+1}
                                                }
                                pc=pc+1
                                }
               }else{}
      if (k==3){ind=cbind(ind,apply(xdiff,1,prod))}else{}
      ind=ind*w;
      ipi=t(ind)%*%ind;
      if (abs(det(ipi)) > (10^(-20))){
        invipi=solve(ipi);
        smmat=invipi%*%t(ind*w);
        heff[i,]=t(curh);
        dep=y[keepind];
        beta=smmat%*%dep;
        varmat=invipi%*%t(ind*w)%*%(ind*w)%*%invipi;
        if (der==0){yhat[i]=beta[1]; varhat[i]=varmat[1,1]
        }else if (der==1){if (k==1){varhat[i,]=varmat[b1,b1];
                                    yhat[i]=beta[b1,];
                                    }else{yhat[i,]=t(beta[b1,]);
                                          varhat[i,]=t(diag(varmat[b1,b1]));
                                          }
        }else if (der==2){if (k==1){yhat[i,]=2*t(beta[b2,]);
                                    varhat[i,]=4*varmat[b2,b2];
                                    }else{varhat[i,]=4*t((diag(varmat[b2,b2])));
                                          yhat[i,]=t(beta[b2]);
                                          varhat[i,]=t((diag(varmat[b2,b2])))
                                          }
        }else{stop("LOCPOLY: I cannot deal with p>3 or its derivatives")}
        cont=0;
      }else{curh=curh*1.1;cont=1}
   }else{
     curh=curh*1.1;
     cont=1;
   }
  }

}
return(list(yhat=yhat,varhat=varhat));
}
}

