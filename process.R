library(ggplot2)
library(plyr)

base_dir = "/media/jon/kelima/jon/Data/StoppedFlow/Durham_140515"
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
myb <- melt(yb,id.vars=c("Time","buffer","lipid","temp","conc","Avg","Var"), variable.name="R",value.name="DF")

#use the plot_function
p <- plot_trace(yb,"PBS","POPCPOPS",37)
p

#how does the sd look over time
syb <- subset(yb, buffer == "PBS" & lipid == "POPCPOPS" & temp == 45)
r <- ggplot( syb, aes(x=Time,y=Avg, group=conc, color=conc)) + geom_line() + theme_minimal(base_size=22) + xlab("Time (s)") + ylab("Var") + scale_color_discrete(name="Concentration")
r + geom_ribbon(aes(ymin=Avg-sd,ymax=Avg+sd),alpha=0.1)

#look at repeats in detail
tset <- subset(myb, buffer == "PBS" & lipid == "POPCPOPS" & temp == 25 & conc==12 )
q <- ggplot(tset,aes(x=Time,y=DF,group=R,color=R))
q + geom_line() + theme_minimal(base_size=18)
