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
  
  #pre write string for ProK-IV import otherwise it crashes
  fileConn<-file(save_str)
  writeLines(",1",fileConn)
  close(fileConn)
  
  #write data to file
  write.table(save_data,save_str,append= TRUE, row.names = FALSE, sep=",", col.names = FALSE)
}


toprocess <- c(c(6,7,8,9,10,13,14,15),seq(21,40,by=1),seq(42,80,by=1))
#toprocess <- c(6)

lapply(toprocess,setOutput,data=yb)

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


