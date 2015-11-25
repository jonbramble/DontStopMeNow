## Taken from example in the minpack guide
library(minpack.lm) #for non-linear least squares

#f <- function(TT, A, B, k1, k2){
#  expr <- expression(B + (A-B)*(1-exp(-1*k1*TT))+(A-B)*(1-exp(-1*k2*TT)))
#  eval(expr)
#}

#j <- function(TT, A, B, k1, k2){
#  expr <- expression(B + (A-B)*(1-exp(-1*k1*TT))+(A-B)*(1-exp(-1*k2*TT)))
#  c(eval(D(expr,"A")),eval(D(expr,"B")),eval(D(expr,"k1")),eval(D(expr,"k2")))
#}

#f <- function(TT, A, k1){
#  expr <- expression(A*(1-exp(-1*k1*TT)))
#  eval(expr)
#}

#j <- function(TT, A, k1){
#  expr <- expression(A*(1-exp(-1*k1*TT)))
#  c(eval(D(expr,"A")),eval(D(expr,"k1")))
#}

f <- function(TT, A, B, k1, k2){
  expr <- expression(A*(1-exp(-1*k1*TT))+B*(1-exp(-1*k2*TT)))
  eval(expr)
}

j <- function(TT, A, B, k1, k2){
  expr <- expression(A*(1-exp(-1*k1*TT))+B*(1-exp(-1*k2*TT)))
  c(eval(D(expr,"A")),eval(D(expr,"B")),eval(D(expr,"k1")),eval(D(expr,"k2")))
}

fcn <- function(p,TT,N,fcall,jcall){
  (N-do.call("fcall", c(list(TT = TT), as.list(p))))
}

fcn.jac <- function(p,TT,N,fcall,jcall){
  -do.call("jcall", c(list(TT = TT), as.list(p)))
}

TT=seq(0,60,length.out=1000)

#initial input data
p <- c(A=8,B=1,k1=0.05,k2=-0.01)
Ndet <- do.call("f", c(list(TT = TT), as.list(p)))

N <- Ndet + rnorm(length(Ndet), sd=.01*max(Ndet))
par(mfrow=c(2,1), mar = c(3,5,2,1))
plot(TT, N, bg = "black", cex = 0.5, main="data")

guess <- c(A=8,B=1,k1=0.1,k2=-0.01)

out <- nls.lm(par = guess, fn = fcn, jac = fcn.jac,
              fcall = f, jcall = j,
              TT = TT, N = N, control = nls.lm.control(nprint=1,maxiter=300))

N1 <- do.call("f", c(list(TT = TT), out$par))
lines(TT, N1, col="blue", lwd=2)

plot(1:(out$niter+1), log(out$rsstrace), type="b",
     main="log residual sum of squares vs. iteration number",
     xlab="iteration", ylab="log residual sum of squares", pch=21,bg=2)
## get information regarding standard errors
summary(out)


