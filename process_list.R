library(plyr)

base_dir = "/media/mbzjpb/data/Experimental/StoppedFlow/Durham_140515"
data.sources = list.files(base_dir,pattern="*.csv",full.names=TRUE)

#function to read in the data files, ignoring the first 34 rows of csv file
#adds the buffer, lipid, temp and conc from the file names only
annotate <- function(f){
  dat <- read.table(f,skip=34,nrows=1000,sep=",")
  #read in more information for dataframe
  name <- basename(f)
  name_list <- unlist(strsplit(name,"['.',' ']")[[1]])
  dat$buffer <- as.factor(name_list[1])
  dat$lipid <- as.factor(name_list[2])
  dat$temp <- as.factor(name_list[3])
  dat$conc <- as.factor(name_list[4])
  return(dat)
}

#does the combination have data?
isData <- function(data,params){
  pl = as.list(params)
  wds <- subset(data,lipid==pl$lipid & conc==as.numeric(pl$conc) & temp==as.numeric(pl$temp) & buffer==pl$buffer)
  setdim <- dim(wds)[1] # find the data length which should be 1000 in this test case
  ret <- FALSE
  if(setdim>0){
    ret <- TRUE
  }
  else {
    ret <- FALSE
  }
  ret
}

setOutput <- function(data,set){
  #extract params
  params=sets[set,]
  pl = as.list(params)
  select <- c(pl$R0,pl$R1,pl$R2,pl$R3,pl$R4,pl$R5)
  print(select)
  
  #subset the data based on params
  wds <- subset(data,lipid==pl$lipid & conc==as.numeric(pl$conc) & temp==as.numeric(pl$temp) & buffer==pl$buffer)
  
  #select the columns needed
  cols <- c("R0","R1","R2","R3","R4","R5")
  dat <- wds[cols]
  
  #weighted mean
  save_data <- wds[c("Time")]
  save_data$avg <- apply(dat,1,weighted.mean, w=select)
  
  #setup labels for file names
  buffer_str <- paste("Buffer:",pl$buffer)
  lipid_str <- paste("Lipid:", pl$lipid)
  temp_str <- paste("Temp:", pl$temp, "C")
  conc_str <- paste("Concentration:", pl$conc)
  title_str <- paste(buffer_str,lipid_str,temp_str, conc_str, sep=", ")
  save_str_name <- paste(pl$row,pl$buffer,pl$lipid,pl$conc,pl$temp,"avg",sep="_")
  save_str <- paste(save_str_name,".csv")
  
  #pre write string for ProK-IV import otherwise it crashes - rem these lines otherwise
  fileConn<-file(save_str)
  writeLines(",1",fileConn)
  close(fileConn)
  
  #write data to file
  write.table(save_data,save_str,append= TRUE, row.names = FALSE, sep=",", col.names = FALSE)
}

yb <- ldply(data.sources,annotate)
#here we correct the colnames
colnames(yb) <- c("Time","R0","R1","R2","R3","R4","R5","buffer","lipid","temp","conc")
myb <- melt(yb,id.vars=c("Time","buffer","lipid","temp","conc"), variable.name="R",value.name="DF")


## Apply to all samples
#make a grid over the possible options
buffer <- c("Tris","PBS")
lipid <- c("POPC","POPCPOPS")
temp <- c(12,25,37,45)
conc <- c(6,12,25,50,100)
repeats <-  c("R0","R1","R2","R3","R4","R5")
samples <- data.frame(expand.grid(repeats,conc,temp,lipid,buffer))
sets <- data.frame(expand.grid(conc,temp,lipid,buffer))
colnames(samples) <- c("repeats","conc","temp","lipid","buffer")
colnames(sets) <- c("conc","temp","lipid","buffer")
sets$row <- seq(1,dim(sets)[1])

#puts the results in the sets df
sets$data <- apply(sets,1, isData, data=myb)

### selection process - manual! 
# default selected
sets$R0 <- sets$data
sets$R1 <- sets$data
sets$R2 <- sets$data
sets$R3 <- sets$data
sets$R4 <- sets$data
sets$R5 <- sets$data

# list of rejected traces
sets[6,]$R0 <- FALSE

sets[7,]$R2 <- FALSE
sets[7,]$R5 <- FALSE

sets[10,]$R2 <- FALSE

sets[13,]$R5 <- FALSE

sets[23,]$R0 <- FALSE

sets[25,]$R0 <- FALSE

sets[31,]$R1 <- FALSE

sets[32,]$R0 <- FALSE

sets[33,]$R0 <- FALSE

sets[34,]$R0 <- FALSE

sets[39,]$R0 <- FALSE

sets[42,]$R0 <- FALSE

sets[51,]$R1 <- FALSE

sets[66,]$R0 <- FALSE

sets[68,]$R0 <- FALSE

sets[71,]$R2 <- FALSE

sets[72,]$R2 <- FALSE

sets[73,]$R1 <- FALSE

sets[76,]$R0 <- FALSE
sets[76,]$R1 <- FALSE
sets[76,]$R2 <- FALSE

sets[78,]$R1 <- FALSE

sets[80,]$R1 <- FALSE


# run the averaging on selected data and write to file
toprocess <- c(c(6,7,8,9,10,13,14,15),seq(21,40,by=1),seq(42,80,by=1))
lapply(toprocess,setOutput,data=yb)





