require(NISTnls)
L1 <- Lanczos1
L1[,"y"] <- round(L1[,"y"],2)
L1

lanc.res <- function(b,data){
  res <- rep(NA, length(data$x))
  res <- b[4]*exp(-b[1]*data$x) + b[5]*exp(-b[2]*data$x) + b[6]*exp(-b[3]*data$x) - data$y
}

lanc.jac <- function(b,data){
  
  expr3 <- exp(-b[1]*data$x)
  expr7 <- exp(-b[2]*data$x)
  expr12 <- exp(-b[3]*data$x)
  J <- matrix(0,nrow=length(data$x),ncol = length(b))
  J[,4] <- expr3
  J[,1] <- -(b[4]*(expr3*data$x))
  J[,5] <- expr7
  J[,2] <- -(b[5]*(expr7*data$x))
  J[,6] <- expr12
  J[,3] <- -(b[6]*(expr12*data$x))
  J
}

bb1 <- c(1,2,3,4,5,6)
bb2 <- c(1,1.1,1.2,0.1,0.1,0.1)
names(bb1) <- c("a1","a2","a3","c1","c2","c3")
names(bb2) <- c("a1","a2","a3","c1","c2","c3")

cat("Sumsquares at bb1 =", as.numeric(crossprod(lanc.res(bb1,data=L1))))

require(numDeriv)
JJ <- lanc.jac(bb1, data=L1)
#JJn <- jacobian(lanc.res,bb1,data=L1)

lancexp <- "y ~ (c1*exp(-a1*x) + c2*exp(-a2*x) + c3*exp(-a3*x))"

tnls1 <- try(nls(lancexp,start=bb1,data=Lanczos1,trace=FALSE))
strwrap(tnls1)

require(nlmrt)
tn1xb1 <- try(nlxb(lancexp,start=bb1,data=Lanczos1,trace=FALSE))
tn1xb1