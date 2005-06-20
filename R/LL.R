"LL" <-
function(xx,zz,yy,xxo,zzo,Vmin=NULL,Vmax=NULL,k,j,h=NULL,kernel=NULL){

if (length(xx)!=length(zz)){stop("JLLProd: lengths of x and z are NOT the same, I won't continue!")
} else {}

if (is.null(kernel)){kernel<-"gauss"
} else {kernel<-kernel}

if (is.null(Vmin)){Vmin <- -1
} else {Vmin<-Vmin}

if (is.null(Vmax)){Vm <- 3
} else {Vm<-Vmax}


vx=sqrt(xx^2+zz^2)/sqrt(xxo^2+zzo^2)
qx=cbind(xx,zz)/vx
Vk=seq(Vmin,Vm,length=k); Vj=seq(Vmin,Vm,length=j)

if (is.null(h)){hx=1.06*sd(Vj*xxo)*(length(Vj*xxo)^-0.2); hz=1.06*sd(Vj*zzo)*(length(Vj*zzo)^-0.2)
} else if (length(h)==1){hx=h; hz=h
} else if (length(h)==2){hx=h[1];hz=h[2]
} else {stop("JLLprod: Vector of bandwidths h is NOT either scalar or a 2x1 vector, I won't continue!")}

mj <- Blocc(xx,zz,yy,ev=cbind(Vj*xxo,Vj*zzo),h=cbind(hx,hz))
rVjxo <- matrix(diag(mj$r),nrow=k,ncol=j,byrow=TRUE)
sij<- matrix(NA,nr=length(yy),nc=j)

for (i in 1:nrow(qx)){
mk <- Blocc(xx,zz,yy,ev=cbind(Vk*qx[i,1],Vk*qx[i,2]),kernel=kernel)
d <- abs(matrix(diag(mk$r),nr=k,nc=j)-rVjxo)
sij[i,]<-matrix(Vk,nr=k,nc=j)[(d==matrix(apply(d,2,min),nr=k,nc=j,byrow=TRUE))]/vx[i]
}

ghat <- apply(matrix(Vj,nr=nrow(qx),nc=j,byrow=TRUE)/sij,1,mean)
mr<-Blocc(xx,zz,yy,ev=cbind(xx,zz),kernel=kernel)
H <- locpoly(y=yy,x=ghat,h=1.06*sd(ghat)*length(ghat)^-0.2,p=1,targmat=ghat,der=0,kernel=kernel);
Hd <- locpoly(y=yy,x=ghat,h=1.06*sd(ghat)*length(ghat)^-0.2,p=1,targmat=ghat,der=1,kernel=kernel);
hhat <- H$yhat; hdhat <- Hd$yhat

return(list(g=ghat,h=hhat,hd=hdhat,r=diag(mr$r)))
}

