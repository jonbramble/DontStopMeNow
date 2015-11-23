library(ggplot2)
library(plyr)
library(reshape2)
library(minpack.lm)

#base_dir = "/media/jon/kelima/jon/Data/StoppedFlow/Durham_140515"
base_dir = "/media/mbzjpb/data/Experimental/StoppedFlow/Durham_140515"
data.sources = list.files(base_dir,pattern="*.csv",full.names=TRUE)

#function to read in the data files, ignoring the first 34 rows of csv file
#adds the buffer, lipid, temp and conc from the file names only
annotate <- function(f){
  dat <- read.table(f,skip=34,nrows=1000,sep=",")
  #read in more information for dataframe
  name <- basename(f)
  name_list <- unlist(strsplit(name,"['.',' ']")[[1]])
  #print(name_list)
  dat$buffer <- as.factor(name_list[1])
  dat$lipid <- as.factor(name_list[2])
  dat$temp <- as.factor(name_list[3])
  dat$conc <- as.factor(name_list[4])
  return(dat)
}

#the trace function grouped by concentration
plot_trace <- function(df,abuffer,alipid,atemp){
  plot_data <- subset(df, buffer == abuffer & lipid == alipid & temp == atemp)
  buffer_str <- paste("Buffer:",abuffer)
  lipid_str <- paste("Lipid:", alipid)
  temp_str <- paste("Temp:", atemp, "C")
  title_str <- paste(buffer_str,lipid_str,temp_str, sep=", ")
  print(title_str)
  p <- ggplot( plot_data, aes(x=Time,y=Avg, group=conc, color=conc)) + geom_line() + ggtitle(title_str)+ theme_minimal(base_size=22) + xlab("Time (s)") + ylab("Fluorescence Intensity") + scale_color_discrete(name="Concentration")
}

yb <- ldply(data.sources,annotate)
#average the repeats
yb$avg <- apply(yb[2:7],1,mean)
#calc the sd of the repeats
yb$sd <- apply(yb[2:7],1,sd)
#here we correct the colnames
colnames(yb) <- c("Time","R0","R1","R2","R3","R4","R5","buffer","lipid","temp","conc","Avg","sd")

#To inspect the individual repeats need a melt version
myb <- melt(yb,id.vars=c("Time","buffer","lipid","temp","conc","Avg","sd"), variable.name="R",value.name="DF")

#use the plot_function
p <- plot_trace(yb,"PBS","POPCPOPS",25)
p

#how does the sd look over time
syb <- subset(yb, buffer == "PBS" & lipid == "POPCPOPS" & temp == 25)
r <- ggplot( syb, aes(x=Time,y=Avg, group=conc, color=conc)) + geom_line() + theme_minimal(base_size=22) + xlab("Time (s)") + ylab("Var") + scale_color_discrete(name="Concentration")
r + geom_ribbon(aes(ymin=Avg-sd,ymax=Avg+sd),alpha=0.1)

#look at repeats in detail
tset <- subset(myb, buffer == "PBS" & lipid == "POPCPOPS" & temp == 25 & conc==100 )
q <- ggplot(tset,aes(x=Time,y=DF,group=R,color=R))
q + geom_point() + theme_minimal(base_size=18)



#simple fitting

#define a suitable test function
exp1 <- function(tt,params){
  params$B + (params$A-params$B)*(1-exp(-1*params$k1*tt))
}

exp2 <- function(tt,params){
  params$B + (params$A-params$B)*(1-exp(-1*params$k1*tt))+(params$A-params$B)*(1-exp(-1*params$k2*tt))
}

#calculate the residuals
residFun <- function(p,observed,tt){
  observed - exp2(tt,p)
}

#run nls inside function

fitTrace <- function(data,params){
  nls.out <- nls.lm(par=params, fn = residFun, observed = data$DF, tt = t)
  summary(nls.out)
  ret <- c(nls.out$par$A,nls.out$par$B,nls.out$par$k1,nls.out$par$k2,nls.out$deviance) #add all other params
}

fitSet <- function(data,params,column){
  obs <- subset(data, R==column)
  nls.out <- nls.lm(par=params, fn = residFun, observed = obs$DF, tt = t)
  summary(nls.out)
  ret <- c(nls.out$par$A,nls.out$par$B,nls.out$par$k1,nls.out$par$k2,nls.out$deviance) #add all other params
}

t=seq(0,60,length.out=1000)

#initial guess
parStart <- list(A=9,B=8,k1=0.01,k2=0.001)

#list of columns
datacols <- c("R0","R1","R2","R3","R4","R5")
fitData <- lapply(datacols,fitSet,data=a,params=parStart)
fitData2 <- lapply(datacols,fitSet,data=tset,params=parStart)

nls.out <- fitSet(tset,parStart,"R5")  # ONLY WORKS IF OUTPUT IS NLS.OUT

#plotting
fitLine <- sapply(t,exp2,nls.out$par)
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,1.4), respect = FALSE)
par(mar = c(0, 4.1, 3, 2.1))
obs <- subset(tset, R=="R5")
plot(obs$Time,obs$DF,xaxt = 'n',ylab="F")
lines(t,fitLine,col="#FF0000")
par(mar = c(4.1, 4.1, 0.3, 2.1))
plot(t,nls.out$fvec,ylab = "Residuals", xlab = "Time (s)")

## Apply to all samples
#make a grid over the possible options
buffer <- c("Tris","PBS")
lipid <- c("POPC","POPCPOPS")
temp <- c(12,25,37,45)
conc <- c(6,12,25,50,100)
repeats <-  c("R0","R1","R2","R3","R4","R5")
samples <- data.frame(expand.grid(buffer,lipid,temp,conc,repeats))
colnames(samples) <- c("buffer","lipid","temp","conc","repeats")

#subset the data from the melt version and run the fit over each sample
#dataSet <- function(params,data){
 #wds <- subset(myb, buffer==params$buffer & lipid==params$lipid & temp==params$temp & conc==params$conc & R==as.character(params$repeats))
 #wds <- subset(data,buffer==params["buffer"])
 #print(dim(wds))
 #fitData <- fitTrace(wds,parStart)
#}

#wps <- dataSet(samples[456,],myb)

#b <- apply(samples,1,dataSet,myb) #eventually apply over data frame of possible values
# problem is that apply and direct call dont seem to exist inside function is the same form. 

# not correct either
processSet <- function(data,buffer,lipid,temp,conc,repeats){
  wds <- subset(data,buffer==buffer & lipid==lipid & conc==conc & temp==temp)
}
abc <- processSet(myb,c(samples[456,]))
