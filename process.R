library(ggplot2)
library(plyr)

base_dir = "/media/jon/kelima/jon/Data/StoppedFlow/Durham_140515"
data.sources = list.files(base_dir,pattern="*.csv",full.names=TRUE)

annotate <- function(f){
  dat <- read.table(f,skip=34,nrows=1000,sep=",")
  #read in more information for dataframe
  name <- basename(f)
  name_list <- unlist(strsplit(name,"['.',' ']")[[1]])
  print(name_list)
  dat$buffer <- as.factor(name_list[1])
  dat$lipid <- as.factor(name_list[2])
  dat$temp <- as.factor(name_list[3])
  dat$conc <- as.factor(name_list[4])
  return(dat)
}

plot_trace <- function(df,abuffer,alipid,atemp){
  print(abuffer)
  print(alipid)
  print(atemp)
  plot_data <- subset(df, buffer == abuffer & lipid == alipid & temp == atemp)
  print(nrow(plot_data))
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
colnames(yb) <- c("Time","R0","R1","R2","R3","R4","R5","buffer","lipid","temp","conc","Avg")


p <- plot_trace(yb,"PBS","POPCPOPS",37)
p

tset <- subset(yb, buffer == "PBS" & lipid == "POPCPOPS" & temp == 25 & conc==25 )
plot(tset$Time,tset$R0)
