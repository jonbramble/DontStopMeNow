pp = list(A=9,B=8,k1=0.05)
simDNoisy <- exp1(t,pp) #+ rnorm(1000,sd=.01)
plot(t,simDNoisy,type="l")
#rs <- residFun(list(A=10,k1=-0.01),simDNoisy,t)