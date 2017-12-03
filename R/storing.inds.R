storing.inds <-
  function (folder, channels=NULL, fourier=TRUE, saturated=TRUE, lets.pullup=TRUE, plotting=FALSE, rawPlot=FALSE) {
    #if (!require("seqinr")) {
    #  install.packages("seqinr")
    #  require("seqinr")
    #}
    #storing the names of the files
    coli <-  c("cornflowerblue", "chartreuse4", "gold2", "red", 
               "orange", "purple")
    setwd(paste(folder))
    listp2 <- dir(folder, "*.fsa")
    FSA <- TRUE
    
    if(length(listp2)==0){
      listp2 <- dir(folder, "*.txt")
      if(length(listp2)==0){
        stop(paste("We have not found files with extension .fsa or .txt. Please\nmake sure this is the right folder:", folder,"\n"), call. = FALSE)
      }
      FSA=FALSE
    }
    # list to store all matrices of 4 columns
    all.inds.mats <- list(NA)
    
    if(FSA){
      print("Reading FSA files")
    }else{
      print("Reading files")
    }
    ####################
    ## initialize the progress bar
    count <- 0
    tot <- length(listp2)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
    #####################
    
    if(FSA){ # if fsa files were read
      for(i in 1:length(listp2)){
        ###################
        count <- count + 1
        ###################
        fsaFile<-read.abif(listp2[i]) # file read
        lens <- lapply(fsaFile$Data, length)
        aaa <- table(unlist(lens)) # length of vectors in data
        #### CONDITION IF CHANNELS WAS NOT SPECIFIED
        if(is.null(channels)){ # find number of channels
          channels <- as.vector(aaa[which(as.numeric(names(aaa)) > 1000 & as.numeric(names(aaa)) < 20000)])
        }
        real.len <- as.numeric(names(aaa)[which(aaa == channels & as.numeric(names(aaa)) > 1000 & as.numeric(names(aaa)) < 20000 )]) # length of the data elements that meet the requirement to be a data file
        v <- as.vector(which(unlist(lens) == real.len)) # length to find, elements storing the data found
        reads<-list(NA) # to store info
        for(j in 1:channels){ # for each channel
          v2 <- v[j]
          reads[[j]] <- fsaFile$Data[[v2]]
        }
        all.inds.mats[[i]] <- matrix(unlist(reads),ncol=channels)
        names(all.inds.mats)[i] <- as.character(listp2[i])
        #attributes(all.inds.mats[[i]]) <- list(FSA=FALSE)
        #attributes(list.ladders[[t]]) <- list(mycomm=names(list.ladders)[t])
        ################################
        setTxtProgressBar(pb, count/tot)### keep filling the progress bar
        ################################
      }
      #attributes(all.inds.mats) <- list(FSA=TRUE)
    }else{ # cqs files
      for(i in 1:length(listp2)){
        ###################
        count <- count + 1
        ###################
        
        # read the 7th line where the line with column names is located
        dodo <-scan(listp2[i],what="a",nlines=7,skip=6,quiet=TRUE) 
        chaco <- grep("FILTER",dodo)# which have name filter
        nchaco <- length(chaco) # number of channels run
        topick <- chaco[1]:(chaco[1] + (nchaco-1)) # where are the lines to read
        
        ds <- read.table(listp2[i], header = FALSE, skip=84)
        all.inds.mats[[i]] <- as.matrix(ds[,topick])
        names(all.inds.mats)[i] <- as.character(listp2[i])
        #attributes(all.inds.mats[[i]]) <- list(FSA=FALSE)
        ################################
        setTxtProgressBar(pb, count/tot)### keep filling the progress bar
        ################################
      }
      #attributes(all.inds.mats[[i]]) <- list(FSA=FALSE)
    }
    
    close(pb) # close the progress bar
    if(fourier == TRUE){ #FOURIER
      print("Applying Fourier tranformation for smoothing...")
      all.inds.mats <- lapply_pb(all.inds.mats, function(x){apply(x, 2, transfft)})
    }
    if(saturated == TRUE){ #SATURATION
      print("Checking and correcting for saturated peaks...")
      all.inds.mats <- lapply_pb(all.inds.mats, function(x){apply(x, 2, saturate)})
    }
    if(FSA){
      if(lets.pullup==TRUE){#PULL UP
        print("Applying pull up correction to the samples to decrease noise from channel to channel")
        if(plotting == TRUE){
          all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels, plotting=TRUE)
        }
        all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels)
      }
    }
    ### ------------------
    ### ------------------
    
    if(rawPlot == TRUE){
      layout(matrix(1:2,2,1))
      coli <- cfp <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
      naname <- c("FAM","HEX","NED","ROX","LIZ")
      
      print("Plotting raw data")
      
      ####################
      ## initialize the progress bar
      count <- 0
      tot <- length(listp2)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
      #####################
      for(i in 1:dim(all.inds.mats[[1]])[2]){
        plot(all.inds.mats[[1]][,i], col=transp(coli[i],0.6), type="l", ylab="RFU", main=paste(naname[i],"channel. Plot",i,"of",dim(all.inds.mats[[1]])[2]), las=2)
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey75")
        
        if(length(all.inds.mats) > 1){
          for(j in 1:length(all.inds.mats)){
            ###################
            count <- count + (1/dim(all.inds.mats[[1]])[2])
            ###################
            lines(all.inds.mats[[j]][,i], col=transp(coli[i],0.4), lwd=.2)
            ################################
            setTxtProgressBar(pb, count/tot)### keep filling the progress bar
            ################################
          }
        }
      }
      close(pb) # close the progress bar
    }
    
    
    layout(matrix(1,1,1))
    return(all.inds.mats)
  }

