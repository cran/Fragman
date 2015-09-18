find.ladder <-
  function(x, ladder, ci.upp=1.96, ci.low=1.96, draw=TRUE, dev=50, warn=TRUE, init.thresh=250, sep.index=8, method="red", avoid=1500, who="sample"){
    ##### stablish a minimum threshold
    thresh=init.thresh
    roxy <- big.peaks.col(x, thresh) # find all peaks with that threshold
    nnn <- length(roxy$pos)
    while(nnn < length(ladder)){
      warning(paste("Warning: Your ladder for individual", who, "is bad, very low intensity was found, we are reducing the threshold 2x trying to find it, don't worry the analysis will keep going anyways"))
      thresh=thresh/2
      roxy <- big.peaks.col(x, thresh) # find all peaks with that threshold
      nnn <- length(roxy$pos)
    }
    ####
    #use the indicator peak at the beggining which is usually the tallest peak in the ladder
    #qw <- which(roxy$hei == max(roxy$hei))
    #roxy <- list(pos=roxy$pos[-qw], hei=roxy$hei[-qw])
    # get rid of very bad peaks noising the ladder
    roxy <- separate(roxy, shift=sep.index, type="pos") 
    #####################
    if(method == "cicor"){
      roxy <- separate(roxy, shift=40, type="pos")
      nnn <- length(ladder)
      mm <- median(roxy$hei)
      vvv <- which(roxy$hei == max(roxy$hei))
      se <- sd(roxy$hei[-vvv])/sqrt(length(roxy$hei[-vvv]))
      reduced <- roxy$pos[which(roxy$hei < mm+(1*se) & roxy$pos > mm-(1*se) & roxy$hei > init.thresh)]
      reduced2 <- roxy$pos[which(roxy$pos >= min(reduced))]
      
      ## if there's still just toomany roxy after the reduced search
      if(length(reduced2) >= nnn & length(reduced2) < nnn+8){
        all.combs <- combn(reduced2, m=length(ladder))
        cors <- apply(all.combs, 2, function(x,y){cor(x,y)}, y=ladder)
        # positions of roxy found
        found <- all.combs[,which(cors == max(cors))]
        roxy <- list(pos=found, hei=x[found], wei=ladder)
        if(draw == TRUE){
          plot(x, type="l", xaxt="n", ylab="Intensity")
          axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
          abline(v=found, col="red", lty=3)
        }
        # reduced search to maximum correlations
      }else{# there's actually no ladder so no good correlation was found
        roxy <- list(pos=seq(1,nnn) + rnorm(nnn,0,1), hei=seq(1,nnn)+ rnorm(nnn,0,1), wei= ladder)
        print("Friend I was not able to find a ladder in this sample or you used the wrong ladder, look the plot")
        plot(x, type="l")
      }
    }
    #####################
    if(method == "red"){
      nnn <- length(ladder)
      
      #### get the length of the peak region
      le <- length(roxy$pos[1]:(roxy$pos[length(roxy$pos)]))
      #### create a vector of expected indexes according to our ladder
      tra <- (le * ladder) / max(ladder)
      #### number of possible tries
      tries <- length(seq(tra[length(tra)], roxy$pos[length(roxy$pos)], 1))+100
      
      # define function returning absolute values
      abso <- function(test1, pos){
        test1 <- matrix(test1, nrow=1)
        res <- apply(test1,2,function(x, y){
          xx1 <- abs(as.vector(x) - y)
          z <- y[which(xx1 == min(xx1))[1]]
          return(z)}, y=pos)
        return(res)
      }
      
      # get all possible tests
      all.tests <- apply(matrix(1:tries,nrow=1), 2, function(q1, q2){q1+q2},q2=tra)
      # get all possible absolute differences
      all.abso <- apply(all.tests, 2, abso, pos=roxy$pos)
      # get all possible correlations
      all.cors <- apply(all.abso, 2, function(x,y){cor(x,y)}, y=ladder)
      # get all sum of squares
      step1 <- (all.tests - all.abso)^2
      all.ss2 <- apply(step1, 2, sum)
      # get all standarized variances and final response
      sss2 <- abs ( (all.ss2 - mean(all.ss2) )/ sd(all.ss2)) # the smaller the better
      response <- all.cors/sss2
      
      vv4 <- which(response >= min(sort(response, decreasing=TRUE)[1:5]))
      #vv4 <- which(response >= .99)
      reduced <- sort(unique(as.vector(all.tests[,vv4])))
      reduced2 <- roxy$pos[which(roxy$pos >= min(reduced))]
      
      ## if there's still just toomany roxy after the reduced search
      if(length(reduced2) >= nnn & length(reduced2) < nnn+8){
        all.combs <- combn(reduced2, m=length(ladder))
        cors <- apply(all.combs, 2, function(x,y){cor(x,y)}, y=ladder)
        # positions of roxy found
        found <- all.combs[,which(cors == max(cors))]
        roxy <- list(pos=found, hei=x[found], wei=ladder)
        if(draw == TRUE){
          
          plot(x, type="l", xaxt="n", ylab="Intensity", ylim=c(0,(max(roxy$hei)+1000)), cex.axis=0.6, las=2)
          axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
          abline(v=found, col="red", lty=3)
          legend("topleft", legend=paste("Correlation:",round(max(cors), digits=4), sep=""), bty="n")
          legend("topright", legend=c("Peak found"), col=c("red"), bty = "n", lty=c(3), lwd=c(1), cex=0.85)
        }
        #which(sss2 == min(sss2))
        # reduced search to maximum correlations
      }else{# there's actually no ladder so no good correlation was found
        roxy <- list(pos=seq(1,nnn) + rnorm(nnn,0,1), hei=seq(1,nnn)+ rnorm(nnn,0,1), wei= ladder)
        print("Friend I don't think you have ladder in this sample or you used the wrong ladder, look the plot")
        plot(x, type="l")
      }
    }
    ### end of method "red"
    #####################
    #####################
    #####################
    if(method == "cor"){ #EXHAUSTIVE CORRELATION
      if(length(roxy$pos) > (length(ladder)+10)){
        print(paste("WOOW too many peaks in this",who,"!! low thresholds throw too many noisy peaks, consider increasing the initial thresold for your ladder, the number of possible combinations is too high to be computed, we will have to do 15,000 samples and get the most likely sizing, you better double check this sample"))
        
        #thresh = init.thresh
        #roxy <- big.peaks.col(x, thresh)
        nono <- which(roxy$pos < avoid)
        roxy <- list(pos = roxy$pos[ -nono], hei = roxy$hei[-nono])
        roxy <- separate(roxy, shift=sep.index, type="pos")
        pos.mod <- matrix(0, ncol=15000, nrow=length(ladder))
        for(k in 1:15000){pos.mod[,k] <- sort(sample(roxy$pos, size=length(ladder), replace=FALSE), decreasing=FALSE)}
        dd <- apply(pos.mod, 2, function(x5, ladder) {cor(x5, ladder)}, ladder)
        v <- which(dd == max(dd))[1]
        v2 <- which(roxy$pos %in% pos.mod[, v])
        roxy <- list(pos = pos.mod[, v], hei = roxy$hei[v2], wei = ladder)
      }else{
        mi=length(ladder)
        if(mi > length(roxy$pos)){
          print(paste("ERROR!! using the initial threshold you specified we did not find enough peaks for", who,", please try reducing the 'init.thresh' or check the plot for this plant using 'detect.ladder' function, maybe this sample did not have ladder :)"))
          roxy <- list(pos=seq(1,length(ladder)) + rnorm(length(ladder),0,1), hei=seq(1,length(ladder))+ rnorm(length(ladder),0,1), wei= ladder)
        }else{
          pos.mod <- combn(roxy$pos, m=mi)
          dd <- apply(pos.mod, 2, function(x5,ladder){cor(x5,ladder)}, ladder)
          v <- which(dd == max(dd))
          v2 <- which(roxy$pos %in% pos.mod[,v])
          roxy <- list(pos=pos.mod[,v], hei=roxy$hei[v2], wei=ladder) 
        }
      }
    }
    # END OF METHOD "cor"
    #######################
    #####################
    #####################
    #####################
    if(method == "ci"){
      z <- which(roxy$hei ==  max(roxy$hei))
      mm <- median(roxy$hei[-z]) # get the median height of the peaks
      se2.low <- (sd(roxy$hei[-z])/sqrt(length(roxy$hei[-z]))) * ci.low # produce the confidnce interval
      se2.upp <- (sd(roxy$hei[-z])/sqrt(length(roxy$hei[-z]))) * ci.upp 
      v <- which( (roxy$hei > (mm-se2.low)) & (roxy$hei < (mm+se2.upp))) # reduce the ladder to the peaks inside the 95% confidence interval
      roxy <- list(pos=roxy$pos[v], hei=roxy$hei[v])
      ## keep looking
      vv <- which(diff(roxy$pos) < dev) 
      vv2 <- vv + 1
      # start condition
      while(length(vv) > 0){
        keep <- numeric()
        for(h in 1:length(vv)){
          a1 <- vv[h]
          a2 <- vv2[h]
          a3 <- c(roxy$hei[a1],roxy$hei[a2])
          a4 <- c(a1,a2)
          keep[h] <- (a4[which(a3 == max(a3))])[1]
        }
        keep <- unique(keep)
        '%!in%' <- function(x,y)!('%in%'(x,y))
        keep2 <- unique(c(vv,vv2)[which(c(vv,vv2) %!in% keep)])
        roxy <- list(pos=roxy$pos[-keep2], hei=roxy$hei[-keep2])
        # check again
        vv <- which(diff(roxy$pos) < dev) 
        vv2 <- vv + 1
      }
      
      s1 <- length(roxy$pos)- (length(ladder) - 1)
      s2 <- length(roxy$pos)
      if(s1 <= 0){
        if(warn==TRUE){
          print("Are you sure this is a ladder channel, we did not find a clear pattern, stop a minute to check the plot")
        }
        if(draw != F){
          plot(x, ylim=c(0,mm+se2.upp), type="l")
        }
        #roxy <- list(pos=roxy$pos, hei=roxy$hei, wei= ladder)
        roxy <- list(pos=seq(1,length(ladder)) + rnorm(length(ladder),0,1), hei=seq(1,length(ladder))+ rnorm(length(ladder),0,1), wei= ladder)
      }else{roxy <- list(pos=roxy$pos[s1:s2], hei=roxy$hei[s1:s2], wei= ladder)}
    }
    ################
    # once ladder is found
    ##############################################################
    if(method == "cor" | method == "ci"){
      if(draw == TRUE){
        xx <- roxy$wei
        yy <- roxy$pos
        mod <- lm(yy~xx)
        xlabels <- as.vector(predict(mod, newdata=data.frame(xx=seq(0,max(ladder),by=25))))
        plot(x, type="l", main=paste("Ladder in channel provided",sep=""), xlab="", ylab="Intensity", xaxt="n")
        if(method == "cor"){
          legend("topleft", legend=paste("Correlation=",round(max(dd)), sep=""), bty="n") 
        }
        axis(1, at=xlabels, labels=seq(0,max(ladder),by=25), cex.axis=0.9)
        abline(v=roxy[[1]], col="red", lty=3)
        if(method == "ci"){
          abline(h=mm, col="blue", lty=3)
          abline(h=(mm+se2.upp), col="blue", lty=3)
          abline(h=(mm-se2.low), col="blue", lty=3)
          legend("topright", legend=c("90% CI", "Peaks found"), col=c("blue", "red"), bty = "n", lty=c(3,3), cex=1)
        }else{legend("topright", legend=c("Peaks found"), col=c("red"), bty = "n", lty=c(3), cex=1)}
        text(x=roxy[[1]], rep(-200, length(ladder)), labels=ladder, cex=0.6)
        #text(x=0, mm, labels="90% CI", cex=0.6)
      }
    }
    return(roxy)
    #########
  }
