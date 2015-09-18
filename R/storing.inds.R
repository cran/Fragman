storing.inds <-
  function (folder, channels=NULL, fourier=TRUE, saturated=TRUE, lets.pullup=TRUE, plotting=FALSE) {
    #if (!require("seqinr")) {
    #  install.packages("seqinr")
    #  require("seqinr")
    #}
    #storing the names of the files
    setwd(paste(folder))
    listp2 <- dir(folder, "*.fsa")
    # list to store all matrices of 4 columns
    all.inds.mats <- list(NA)
    
    print("Reading FSA files")
    
    ####################
    ## initialize the progress bar
    count <- 0
    tot <- length(listp2)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
    #####################
    
    for(i in 1:length(listp2)){
      ###################
      count <- count + 1
      ###################
      fsaFile<-read.abif(listp2[i]) # file read
      lens <- lapply(fsaFile$Data, length)
      aaa <- table(unlist(lens)) # length of vectors in data
      #### CONDITION IF CHANNELS WAS NOT SPECIFIED
      if(is.null(channels)){
        channels <- as.vector(aaa[which(as.numeric(names(aaa)) > 1000)])
      }
      real.len <- as.numeric(names(aaa)[which(aaa == channels & as.numeric(names(aaa)) > 1000 )]) # length of the data elements that meet the requirement to be a data file
      v <- as.vector(which(unlist(lens) == real.len)) # length to find, elements storing the data found
      reads<-list(NA) # to store info
      for(j in 1:channels){
        v2 <- v[j]
        reads[[j]] <- fsaFile$Data[[v2]]
      }
      all.inds.mats[[i]] <- matrix(unlist(reads),ncol=channels)
      names(all.inds.mats)[i] <- as.character(listp2[i])
      ################################
      setTxtProgressBar(pb, count/tot)### keep filling the progress bar
      ################################
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
    if(lets.pullup==TRUE){#PULL UP
      print("Applying pull up correction to the samples to decrease noise from channel to channel")
      if(plotting == TRUE){
        all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels, plotting=TRUE)
      }
      all.inds.mats <- lapply_pb(all.inds.mats, pullup, channel=channels)
    }
    return(all.inds.mats)
  }

