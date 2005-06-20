"JLL" <-
function(xx,zz,yy,z0,r0,ngrid=NULL,h1=NULL,h2=NULL,hstar=NULL,
                k1=NULL,k2=NULL,kstar=NULL,p1=NULL,p2=NULL,pstar=NULL){

if (is.null(ngrid)){xxe<- seq(min(xx),max(xx),length=15); zze<- unique(sort(c(seq(min(zz),max(zz),length=14),z0)))
}else {xxe <- seq(min(xx),max(xx),length=ngrid); zze <- unique(sort(c(seq(min(zz),max(zz),length=ngrid-1),z0)))}

if (length(xx)!=length(zz)){stop("JLLProd: lengths of x and z are NOT the same, I won't continue!")
}else{}

if (is.null(h1)){hx=1.06*sd(xx)*(length(xx)^-0.2); hz=1.06*sd(zz)*(length(zz)^-0.2)
}else if (length(h1)==1){hx=h1; hz=h1
}else if (length(h1)==2){hx=h1[1];hz=h1[2]
}else {stop("JLLprod: Vector of bandwidths h1 is NOT either scalar or a 2x1 vector, I won't continue!")}

klist<-list("uniform","triangular","quartic","epanech","triweight","gauss","order34","order56","order78")
if (is.null(k1)){k1="gauss"
}else if (any(klist==k1)){k1=k1
}else {stop("JLLProd: You must select a VALID k1 kernel, see the manual")}

R <- locpoly(y=yy,x=cbind(xx,zz),h=cbind(hx,hz),p=p1,targmat=expand.grid(xxe,zze),der=0,kernel=k1);
S <- locpoly(y=yy,x=cbind(xx,zz),h=cbind(hx,hz),p=p1,targmat=expand.grid(xxe,zze),der=1,kernel=k1);

if (is.null(h2)){hr=1.06*sd(R$yhat)*(nrow(R$yhat)^-0.2)
}else if (length(h2)==1){hr=h2; hz=h2
}else if (length(h2)==2){hr=h2[1];hz=h2[2]
}else {stop("JLLprod: Vector of bandwidths h1 is NOT either scalar or a 2x1 vector, I won't continue!")}

if (is.null(k2)){k2="gauss"
}else if (any(klist==k2)){k2=k2
}else {stop("JLLProd: You must select a VALID k2 kernel, see the manual")}

Q <- locpoly(y=S$yhat[,2],x=cbind(R$yhat,expand.grid(xxe,zze)[,2]),p=p2,h=cbind(hr,hz),der=0,kernel=k2)
grid <- matrix(R$yhat,nr=length(xxe),nc=length(zze)); q <- matrix(Q$yhat,nr=length(xxe),nc=length(zze));
grid <- grid[,(zze==z0)]; qinvs <- (q[,(zze==z0)])^(-1)

prhat<-matrix(R$yhat,nr=length(xxe),nc=length(zze))

dat<-cbind(grid,qinvs); o <- order(dat[,1])
dat<-dat[o,]; grid <- dat[,1]; qinvs <- dat[,2]
ng<-length(qinvs)
mm1 = cbind(qinvs[1:(ng-1)],qinvs[2:ng]);

int1 = apply(mm1,1,max)*(grid[2:ng]-grid[1:(ng-1)]); int1 =c(int1,0);
int2 = apply(mm1,1,min)*(grid[2:ng]-grid[1:(ng-1)]); int2 = c(int2,0);

m<-(qinvs[2]-qinvs[1])/(grid[2]-grid[1]); mm<-(qinvs[ng]-qinvs[ng-1])/(grid[ng]-grid[ng-1])

Mhat<-matrix(NA,nr=length(xxe),nc=length(zze))

for (i in 1:nrow(prhat)){
      for (k in 1:ncol(prhat)){
          if (r0 < grid[1]){ # Case A#
                    if (prhat[i,k] <= r0){ # Case A-1#
                       Mhat[i,k] = -((r0-prhat[i,k])*(qinvs[1]-m*(grid[1]-r0))-(m/2)*(r0 - prhat[i,k])^2)
                       }
                    else if ((r0 < prhat[i,k])&(prhat[i,k] <= grid[1])){ # Case A-2#
                           Mhat[i,k] = (prhat[i,k]-r0)*(qinvs[1]-m*(grid[1]-prhat[i,k]))-(m/2)*(prhat[i,k]-r0)^2 ;
                           }
                    else if ((grid[1] < prhat[i,k])&(prhat[i,k] <= grid[ng])){ # Case A-3#
                           j = sum((grid <= prhat[i,k])) ;
                           Mhat[i,k] = (grid[1]-r0)*qinvs[1] - (m/2)*(grid[1]-r0)^2 + sum((int1[1:j]+int2[1:j])/2) ;
                           }                      
                    else {# Case A-4#
                           Mhat[i,k] = (grid[1]-r0)*qinvs[1] - (m/2)*(grid[1]-r0)^2 + sum((int1[1:ng]+int2[1:ng])/2) + (prhat[i,k]-grid[ng])*qinvs[ng] + (mm/2)*(prhat[i,k]-grid[ng])^2 ;
                           }
                    }
          else if ((grid[1] <= r0)&(r0 <= grid[ng])){ # Case B#
                 jj = sum((grid <= r0)) ;
                    if (prhat[i,k] < grid[1]){ # Case B-1#
                       Mhat[i,k] = -((grid[1]-prhat[i,k])*qinvs[1] - (m/2)*(grid[1]-prhat[i,k])^2 + sum((int1[1:jj]+int2[1:jj])/2)) ;
                       }
                    else if ((grid[1] <= prhat[i,k])&(prhat[i,k] < r0)){ j = sum((grid <= prhat[i,k])) ;   # Case B-2#
                           Mhat[i,k] = - sum((int1[j:jj]+int2[j:jj])/2) ;
                           }
                    else if ((r0 <= prhat[i,k])&(prhat[i,k] <= grid[ng])){ j = sum((grid <= prhat[i,k])) ;   # Case B-3#
                           Mhat[i,k] = sum((int1[jj:j]+int2[jj:j])/2) ;
                           }
                    else { # Case B-4#
                           Mhat[i,k] = (prhat[i,k]-grid[ng])*qinvs[ng] + (mm/2)*(prhat[i,k]-grid[ng])^2 + sum((int1[jj:ng]+int2[jj:ng])/2) ;
                           }
                    }

          else {   # Case C#
                    if (prhat[i,k] < grid[1]){ # Case C-1#
                       Mhat[i,k] = -((grid[1]-prhat[i,k])*qinvs[1] - (m/2)*(grid[1]-prhat[i,k])^2 + (r0-grid[ng])*qinvs[ng] + (mm/2)*(r0-grid[ng])^2 + sum((int1[1:ng]+int2[1:ng])/2)) ;
                       }
                    else if ((grid[1] <= prhat[i,k])&(prhat[i,k] < grid[ng])){ j=sum((grid <= prhat[i,k])) ;   # Case C-2#
                           Mhat[i,k] = -(sum((int1[j:ng]+int2[j:ng])/2)+(r0-grid[ng])*qinvs[ng]+(mm/2)*(r0-grid[ng])^2) ;
                           }
                    else if ((grid[ng] <= prhat[i,k])&(prhat[i,k] < r0)){ # Case C-3#
                           Mhat[i,k] = -((r0-prhat[i,k])*(qinvs[ng] + mm*(prhat[i,k]-grid[ng])) + (mm/2)*(r0-prhat[i,k])^2) ;
                           }
                    else {  # Case C-4#
                           Mhat[i,k] = (prhat[i,k]-r0)*(qinvs[ng]+(mm/2)*(r0-grid[ng])) + (mm/2)*(prhat[i,k]-r0)^2 ;
                         }

               }    
                       
      }
      }
Ghat=apply(Mhat,1,mean); Fhat=apply(Mhat,2,mean)

if (is.null(hstar)){hh<-1.06*sd(diag(Mhat))*(length(diag(Mhat))^-0.2)
}else if (length(hstar)==1){hh=hstar
}else {stop("JLLProd: hstar must be a scalar")}

if (is.null(kstar)){kstar="gauss"
}else if (any(klist==kstar)){kstar=kstar
}else {stop("JLLProd: You must select a VALID kstar kernel, see the manual")}

H <- locpoly(y=diag(prhat),x=diag(Mhat),h=hh,targmat=diag(Mhat),p=pstar,der=0,kernel=kstar)

return(list(Rhat=prhat,Mhat=Mhat,Ghat=Ghat,Fhat=Fhat,hhat=H$yhat,x=xxe,z=zze))
}

