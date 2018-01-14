ladder.corrector <- function(stored, to.correct, ladder, thresh=200, 
                             env = parent.frame(), ...){
  list.data <- env$list.data.covarrubias
  
  vvv <- which(names(list.data) %in% to.correct)
  chan <- dim(stored[[1]])[2]
  ###
  #--
  ###
  for(j in 1:length(vvv)){
    cat(paste("\nSTARTING SAMPLE",j,""))
    www <- vvv[j]
    x <- stored[[www]][,chan]
    roxy <- big.peaks.col(x, thresh)
    
    outlier_values <- boxplot.stats(roxy$hei)$out
    if(length(outlier_values)>0){
      ylimi <- max(sort(roxy$hei[which(roxy$hei != outlier_values)], decreasing=TRUE)) 
    }else{
      ylimi <- max(sort(roxy$hei, decreasing=TRUE)) 
    }
        
    layout(matrix(1,1,1))
    plot(x, type="l", ylim=c(0,ylimi), main=paste("Provitional plot for\n",to.correct[j]), las=2, cex.main=.8, ...);grid()
    abline(h=c(20,50,100,200), col="red",lty=3)
    axis(2, at=c(20,50,100,200),labels = c(20,50,100,200), las=1, col.ticks ="red")
    cat("\nPlease specify if you prefer manual or automatic scoring. \nType the character 'm' or 'a' (autom is preferred):\n")
    type <- readline()
    type <- gsub("\'","",type) # in case the user gave character
    type <- gsub("\"","",type)
    
    basicc <- c("a","m")
    choco <- length(intersect(basicc,type))
    while((choco)==0){
      cat("\nType the character 'm' or 'a':\n")
      type <- readline()
      type <- gsub("\'","",type) # in case the user gave character
      type <- gsub("\"","",type)
      choco <- length(intersect(basicc,type))
    }
    
    if(type=="a"){
      cat("Please provide the threshold value for the ladder in RFUs (number):")
      hisladder <- as.numeric(readline())
      cat("Should we reduce the search in the x-axis to certain peaks?.If none type the 'Esc' key\notherwise click over the plot where where the start is (one click) and \nthe end (another click), and finally type the 'Esc' key")
      #tolook <- as.numeric(readline())
      mine <- locator(type="p", pch=20, col="red")$x
      if(is.null(mine)){
        reducing <- NULL
      }else{
        tolook <- round(min(mine)):round(max(mine))
        reducing <- tolook
      }
      
      provi <- find.ladder(x, ladder=ladder,init.thresh = hisladder, 
                           draw = FALSE, method="iter2", reducing = reducing)
      abline(v=provi$pos, col="red",lty=3)
      list.data[[www]]$pos <- provi$pos
      list.data[[www]]$hei <- provi$hei
      list.data[[www]]$wei <- provi$wei
      list.data[[www]]$corr <- provi$corr
      cat("\nSample adjusted.")
    }else if(type=="m"){
      cat("Please provide the threshold value for the ladder in RFUs (number):")
      hisladder <- as.numeric(readline())
      x <- stored[[www]][,chan]
      roxy <- big.peaks.col(x, hisladder)
      cat(paste("\nYou need you to click over peaks (",length(ladder),") corresponding to the following ladder: \n\n", paste(ladder, collapse = " "), "\n\nPress 'esc' when you are done. \n"))
      good <- locator(type="p", pch=20, col="red")$x
      #############################################
      ## now use such values to correct the model
      corrected <- numeric()
      for(k in 1:length(good)){
        yoyo <- abs(roxy$pos - good[k])
        corrected[k] <- roxy$pos[which(yoyo == min(yoyo))]
      }
      corrected.hei <- roxy$hei[which(roxy$pos %in% corrected)]
      if(length(corrected) != length(ladder)){
        cat("You did not select the same number of peaks than the number of elements \nof your ladder. Getting more likely peaks. Please check the resulting peaks \nand repeat if needed. ")
      }else{
        list.data[[www]]$pos <- corrected
        list.data[[www]]$hei <- corrected.hei
        list.data[[www]]$wei <- ladder
        list.data[[www]]$corr <- cor(ladder, corrected)
      }
      ##############################################
    }
    ###
  }
  ###
  #--
  ###
  layout(matrix(1:3,nrow = 3,ncol=1))
  for(l in 1:length(vvv)){
    www <- vvv[l]
    x <- stored[[www]][,chan]
    limi <- sort(list.data[[www]]$hei, decreasing = TRUE)
    plot(x, type="l", xaxt="n", ylim=c(0,(limi[3]+100)), cex.axis=0.6, las=2,  col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, cex.main=.8,main=paste("Adjusted ladder for\n",to.correct[l]))#xlim=c((min(list.data[[www]]$pos)-100),(max(list.data[[www]]$pos)+100)),
    axis(1, at=list.data[[www]]$pos, labels=list.data[[www]]$wei, cex.axis=0.6)
    points(x=list.data[[www]]$pos, y=list.data[[www]]$hei,cex=1.1, col=transp("black",0.85))
    points(x=list.data[[www]]$pos, y=list.data[[www]]$hei, pch=20, col=transp("red",0.7))
    legend("topleft", legend=paste("Correlation:",round(list.data[[www]]$corr, digits=4), sep=""), bty="n")
    legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
    
  }
  env$list.data.covarrubias <- list.data
  layout(matrix(1,1,1))
  cat("\nJOB DONE!!! ladder has been adjusted for the samples provided")
}