"JLLp" <- 
function(lnY,lnK,lnL,theta,model){
lnk<-lnK-lnL
if (model==2){result <- nls(lnY~b0+b1*(a*lnk+log(exp(lnL)+g))+b2*(a*lnk+log(exp(lnL)+g))^2,start=theta)
} else if (model==3){result <- nls(lnY~b0+b1*(a*lnk+lnL)+b2*(a*lnk+lnL)^2,start=theta)
} else{stop("JLLprod: Your model must be either 2 or 3")}
return(result)
}
