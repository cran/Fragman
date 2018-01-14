find.ladder <-function(x, ladder, draw=TRUE, dev=50, warn=TRUE, init.thresh=NULL, sep.index=8, method=NULL, reducing=NULL, who="sample", attempt=10, cex.title=0.8){
  #roundUP <- function(x) 10^ceiling(log10(x))
  
  hohoho <- big.peaks.col(x[1:length(x)], 100)
  
  # plot(x, type="l", ylim = c(0,2000))
  # abline(v=hohoho$pos, lty=3, col="red")
  
  if(is.null(method)){
    ## check if the ladder passes the length of x of is a small ladder
    tototo <- big.peaks.col(x[1:length(x)], median(hohoho$hei)*2)
    close <- (length(x) - tototo$pos[length(tototo$hei)])
    if(close < 100){
      method="iter"
    }else{
      method="iter2"
    }
  }
  
  if(is.null(init.thresh)){
    if(mean(x) > 1000){
      #hohoho2 <- big.peaks.col(x[1:length(x)], 2*median(hohoho$hei))
      init.thresh <- quantile(hohoho$hei,.95)
    }else{
      init.thresh <- median(hohoho$hei)/2.2#quantile(hohoho$hei,.32)
    }
  }
  
  MSE <- function(x,y){
    X <- cbind(1, x)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% y
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X%*%beta
    res <- list(mse=sum((y-fit)^2), beta=beta[2])
    return(res)
  }
  MSE2 <- function(x,y){
    X <- cbind(1, x)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% y
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X%*%beta
    mse <- sum((y-fit)^2)# is actually SSE
    sst <- sum((y-mean(y))^2)
    r2 <- 1 - (mse/sst) #cor(x,y)^2#mse/sum((y-mean(y,na.rm=TRUE))^2)
    res <- list(mse=mse, beta=beta, r2=r2)
    return(res)
  }
  MSE3 <- function(mix,miy){
    X <- cbind(1, mix)
    qr.X <- qr(X)
    b <- t(qr.Q(qr.X)) %*% miy
    R <- qr.R(qr.X)
    beta <- as.vector(backsolve(R, b))
    fit <- X%*%beta
    mse <- sum((miy-fit)^2)# is actually SSE
    sst <- sum((miy-mean(miy))^2)
    r2 <- 1 - (mse/sst) #cor(x,y)^2#mse/sum((y-mean(y,na.rm=TRUE))^2)
    res <- list(mse=mse, beta=beta, r2=r2)
    return(res)
  }
  ##### stablish a minimum threshold
  thresh=init.thresh
  roxy <- big.peaks.col(x[1:length(x)], thresh) # find all peaks with that threshold
  
  if(!is.null(reducing)){
    nono <- which(roxy$pos %in% reducing)
    roxy <- lapply(roxy, function(x){x[nono]})
    #roxy <- list(pos = roxy$pos[ -nono], hei = roxy$hei[-nono])
  }
  
  nnn <- length(roxy$pos)
  
  # plot(x, type="l", ylim = c(0,2000))
  # abline(v=roxy$pos, lty=3, col="red")
  
  fdsa=1
  while(nnn < length(ladder)){
    #print(paste("Warning: Your ladder for the sample",attr(x,"name") ,"is bad, very low intensity in RFU was found, we are reducing the threshold 2x trying to find it, don't worry too much the analysis will keep going anyways"))
    
    if(fdsa == 1){
      cat("\nReducing threshold 2x to find ladder \n")
      #cat("There's too many ladder peaks expected and not many peaks found in this sample\nPlease make sure all your ladder peaks specified show up in the sample\nOtherwise provide a smaller ladder so all peaks can be found\n")
    }
    thresh=thresh/2
    roxy <- big.peaks.col(x[1:length(x)], thresh) # find all peaks with that threshold
    nnn <- length(roxy$pos)
    fdsa<-fdsa+1
  }
  ####
  #use the indicator peak at the beggining which is usually the tallest peak in the ladder
  # get rid of very bad peaks noising the ladder
  whot <- length(roxy$pos)*.2
  what <- which(roxy$hei == max(roxy$hei))
  #what <- which(roxy$hei == sort(roxy$hei, decreasing = TRUE)[2])
  # if you tallest peak is far beyond don't adjust
  #if(what > whot){
  #  roxy <- roxy
  #}else{ max(roxy$hei)
  roxy <- separate(roxy, shift=sep.index, type="pos")
  
  
  #} plot(x, type="l",xlim=c(2000,7000)); abline(v=roxy$pos, col="red", lty=3)
  #####################
  if (method == "iter") {
    ii <- which(roxy$hei == max(roxy$hei)) + 1
    iii <- length(roxy$hei)
    roxy <- list(pos = roxy$pos[ii:iii], hei = roxy$hei[ii:iii])
    
    step1 <- combn(roxy$pos[1:attempt], 3)
    step2 <- apply(step1/10, 2, MSE, y = ladder[1:3])
    mse <- unlist(lapply(step2, function(x) {
      x$mse
    }))
    covs <- apply(step1, 2, function(x, y) {
      cov(x, y)
    }, y = ladder[1:3])
    step2 <- mse * covs
    step3 <- step1[, which(step2 < sort(step2, decreasing = FALSE)[20])]
    step4 <- apply(step3, 2, function(x, y) {
      which(y %in% x)
    }, y = roxy$pos)
    caller <- function(roxy, www, ladder.call, x) {
      threshold <- length(x)
      posi <- numeric()
      fact2 <- length(ladder.call)
      expect <- roxy$pos[www]
      xxx <- ladder.call[c(1:3)]
      modx <- lm(expect ~ poly(xxx, degree = 1))
      expecto <- predict(modx, data.frame(xxx = ladder.call))
      ladder.call <- ladder.call[which(expecto < threshold * 
                                         0.85)]
      available <- length(roxy$pos) - length(ladder.call)
      ava2 <- length(ladder.call) - abs(available)
      if (available > 0) {
        if ((length(ladder.call) - 1) < 3) {
          tope <- length(ladder.call)
        }
        else {
          tope <- length(ladder.call) - 1
        }
      }
      else {
        tope <- ava2 - 2
      }
      expect <- rep(NA, tope + 1)
      for (i in 3:tope) {
        if (i == 3 & i != tope) {
          expect[1:3] <- roxy$pos[www]
          xxx <- ladder.call[c(1:3)]
          mod <- MSE2(xxx, expect[1:3])
          beta <- (mod)$beta
          expecto <- as.vector(beta[1] + matrix(ladder.call) %*% 
                                 beta[-1])
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo,na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i > 3 & i <= 5) {
          xx <- ladder.call[c(1:i)]
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          expecto <- as.vector(beta[1] + matrix(ladder.call) %*% 
                                 beta[-1])
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          #print(min(yoyo))
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i > 5) {
          xx <- cbind(ladder.call[c(1:i)], ladder.call[c(1:i)]^2, 
                      ladder.call[c(1:i)]^3, ladder.call[c(1:i)]^4)
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          if (length(which(is.na(beta))) > 0) {
            beta[which(is.na(beta))] <- 0
          }
          toto <- cbind(matrix(ladder.call), matrix(ladder.call)^2, 
                        matrix(ladder.call)^3, matrix(ladder.call)^4)
          expecto <- cbind(rep(1, dim(toto)[1]), toto) %*% 
            beta
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo,na.rm = TRUE))
          expect[i + 1] <- act[good]
          if (is.na(mod$r2)) {
            mod$r2 <- 0.1
          }
          if (mod$r2 < 0.9) {
            i = tope
          }
        }
        if (i == tope & i != 3) {
          if (i < 5) {
            expect[1:3] <- roxy$pos[www]
            xx <- ladder.call[c(1:i)]
          }
          else {
            xx <- cbind(ladder.call[c(1:i)], ladder.call[c(1:i)]^2, 
                        ladder.call[c(1:i)]^3, ladder.call[c(1:i)]^4)
          }
          mod <- MSE2(xx, expect[1:i])
          beta <- (mod)$beta
          if (length(which(is.na(beta))) > 0) {
            beta[which(is.na(beta))] <- 0
          }
          if (i < 5) {
            toto <- cbind(matrix(ladder.call))
          }
          else {
            toto <- cbind(matrix(ladder.call), matrix(ladder.call)^2, 
                          matrix(ladder.call)^3, matrix(ladder.call)^4)
          }
          expecto <- cbind(rep(1, dim(toto)[1]), toto) %*% 
            beta
          act <- roxy$pos[-which(roxy$pos %in% expect)]
          yoyo <- abs(expecto[i + 1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          expect[i + 1] <- act[good]
        }
        if (i == tope & i == 3) {
          expect[1:3] <- roxy$pos[www]
        }
      }
      posi <- expect
      tutu <- abs(length(x) - posi)
      posi <- posi[1:which(tutu == min(tutu, na.rm = TRUE))]
      heii <- roxy$hei[which(roxy$pos %in% posi)]
      fact3 <- length(posi)/fact2
      if (length((posi)) < 6) {
        fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                             poly(posi, degree = length((posi)) - 1)))$r.squared * 
          fact3
      }
      else {
        fact <- summary(lm(ladder.call[1:length(posi)] ~ 
                             poly(posi, degree = 5)))$r.squared * fact3
      }
      roxy2 <- list(pos = posi, hei = heii, wei = ladder.call[1:length(posi)], 
                    corr = abs(cor(ladder.call[1:length(posi)], posi)), 
                    error = fact)
      return(roxy2)
    }
    rt <- apply(data.frame(step4), 2, FUN = caller, roxy = roxy, 
                ladder.call = ladder, x = x)
    corrs3 <- unlist(lapply(rt, function(x) {
      x$error
    }))
    roxy3 <- rt[[which(corrs3 == max(corrs3))]]
    if (draw == TRUE) {
      limi <- sort(roxy3$hei, decreasing = TRUE)
      plot(x, type = "l", xaxt = "n", ylim = c(0, (limi[3] + 
                                                     1000)), cex.axis = 0.6, las = 2, xlim = c((min(roxy3$pos) - 
                                                                                                  100), (max(roxy3$pos) + 100)), col = transp("grey35", 
                                                                                                                                              0.7), ylab = "RFU", xlab = "", lwd = 2, main = attributes(x)$mycomm, 
           cex.main = cex.title)
      axis(1, at = roxy3$pos, labels = roxy3$wei, cex.axis = 0.6)
      points(x = roxy3$pos, y = roxy3$hei, cex = 1.1, col = transp("black", 
                                                                   0.85))
      points(x = roxy3$pos, y = roxy3$hei, pch = 20, col = transp("red", 
                                                                  0.7))
      legend("topleft", legend = paste("Correlation:", 
                                       round(roxy3$corr, digits = 4), sep = ""), bty = "n")
      legend("topright", legend = c("Peaks selected"), 
             col = c("red"), bty = "n", pch = c(20), cex = 0.85)
    }
    roxy <- roxy3
  }
  if(method == "iter2"){
    #caller
    last <- length(roxy$pos)
    lld <- length(ladder)
    if((last-attempt) < 0){
      #roxy2 <- roxy#list(pos=posi, hei=heii, wei=ladder.call[1:length(posi)], corr=abs(cor(ladder.call[1:length(posi)],posi)), error=fact)#sum(error, na.rm=TRUE))
      roxy$wei <- ladder[1:length(roxy$pos)]
      roxy$corr <- 0
      roxy$error <- 0
      if(draw == TRUE){
        limi <- sort(roxy$hei, decreasing = TRUE)
        plot(x, type="l", xaxt="n", ylim=c(0,(limi[3]+1000)), cex.axis=0.6, las=2,  col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, main=attributes(x)$mycomm, cex.main=cex.title)
        axis(1, at=roxy$pos, labels=roxy$wei, cex.axis=0.6)
        points(x=roxy$pos, y=roxy$hei,cex=1.1, col=transp("black",0.85))
        points(x=roxy$pos, y=roxy$hei, pch=20, col=transp("red",0.7))
        legend("topleft", legend=paste("Correlation:",round(roxy$corr, digits=4), sep=""), bty="n")
        legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
        legend("center", legend=c("Intensity too low to be detected"), col=c("red"), bty = "n", cex=0.85)
      }
    }else{
    # by default attempt=10 making all combinations of 1st 10 peaks
    step1 <- combn(roxy$pos[last:(last-attempt)],3)
    #### check MSE
    step2 <- apply(step1/10, 2, MSE, y=ladder[lld:(lld-2)]) # function(x,y){sum((summary(lm(I(x)~y))$residuals)^2)}create a custom function with matrices to extract residuals and calculate MSE
    #### mean squared errors
    mse <- unlist(lapply(step2, function(x){x$mse}))
    ##### covariances
    covs <- apply(step1, 2, function(x,y){cov(x,y)}, y=ladder[lld:(lld-2)])
    ##### scores for selection
    step2 <- mse #* covs 
    ##### selected based on scores
    step3 <- step1[,which(step2 < sort(step2, decreasing=FALSE)[20])] # 15 models with least MSE
    # winner = 1
    # for(i in 1:dim(step3)[2]){ # see selected
    #   plot(x,type="l", main=i)
    #   abline(v=step3[,i], lty=3, col="red")
    # }
    # position of peaks-combinations selected 
    step4 <- apply(step3,2,function(x,y){sort(which(y %in% x), decreasing = TRUE)}, y=roxy$pos) # which peaks are
    #caller  6
    ############
    #-----------
    ############
    caller <- function(roxy,www, ladder.call,x){
      threshold <- length(x)
      
      last3 <- length(ladder.call):(length(ladder.call)-2)
      posi <- numeric()
      fact2 <- length(ladder.call)
      # www <- step4[,5]; ladder.call=ladder
      ############################################################
      ## short initial model to avoid people using a ladder.call too long, we cut the ladder.call if necessary
      expect <- roxy$pos[www] # initial positions for this combo (last)
      xxx <- ladder.call[last3] # ladder positions assumed (last)
      modx <- lm(expect~poly(xxx, degree=1))
      ## assuming model is correct predict where should be the rest of the ladder
      expecto <- predict(modx, data.frame(xxx=ladder.call))# + facto
      #ladder.call <- ladder.call[which(expecto < threshold*.85)]
      
      available <- length(roxy$pos) - length(ladder.call) # there's x extra peaks
      ava2 <-  length(ladder.call) - abs(available) # is it more than the user specified??
      if(available < 0){ # if there is less peaks available than existing in ladder BAD USER
        #if((length(ladder.call)-1) < 3){
        #  tope <- length(ladder.call)
        #}else{tope <- length(ladder.call)-1}
        tope <- 3 # is bad sample
        
      }else{
        tope <- length(ladder.call)
      } # if there is more peaks available than expecting in ladder GOOD USER
      
      #if(available > 0){tope <- length(ladder.call)-1}
      #}else{tope <- abs(available)-1}
      #####
      #tope <- length(ladder.call) #- 3
      expect <- rep(NA,tope)# (length(ladder.call)-maxo)
      
      lenlad <- length(ladder.call)
      
      for(i in 3:(tope-1)){ # for all possible peaks
        if(i == 3 & i != tope){# initial step
          
          expect[tope:(tope-2)] <- roxy$pos[www]
          xxx <- ladder.call[last3]
          #mod1 <- lm(expect[1:3]~poly(xxx, degree=1))
          #expecto <- predict(mod1, data.frame(xxx=ladder.call))# + facto
          mod <- MSE3(mix=xxx, miy=expect[tope:(tope-2)])
          beta <- (mod)$beta # extract intercept and slope
          ## predict the rest of the ladder
          expecto <- as.vector(beta[1] + matrix(sort(ladder.call, decreasing = TRUE)) %*% beta[-1])
          # remove missing data
          condo <- sort(expect[which(!is.na(expect))], decreasing = TRUE)
          # keep roxy values that were picked already
          act <- sort(roxy$pos[-which(roxy$pos %in% condo)], decreasing = TRUE)
          # see which is the next most likely value
          yoyo <- abs(expecto[i+1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          
          qwer <- i
          qwer2 <- length(expect)-qwer
          expect[qwer2] <- act[good]
          # --------------------------------
          #if(summary(mod1)$r.squared < .9){
          if(mod$r2 < .9){
            i= tope #length(ladder.call) - maxo
          }
          # --------------------------------
        }else if(i > 3 & i <= 5 ){
          
          xx <- ladder.call[c(lenlad:(lenlad-(i-1)))]
          #mod1 <- lm(expect[1:i]~poly(xx, degree=1))
          #expecto <- predict(mod1, data.frame(xx=ladder.call))
          
          mod <- MSE3(mix=xx, miy=expect[tope:qwer2])
          beta <- (mod)$beta # extract intercept and slope
          ## predict the rest of the ladder
          expecto <- as.vector(beta[1] + matrix(sort(ladder.call, decreasing = TRUE)) %*% beta[-1])
          # remove missing data
          condo <- sort(expect[which(!is.na(expect))], decreasing = TRUE)
          # keep roxy values that were picked already
          act <- sort(roxy$pos[-which(roxy$pos %in% condo)], decreasing = TRUE)
          # see which is the next most likely value
          yoyo <- abs(expecto[i+1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          qwer <- i
          qwer2 <- length(expect)-qwer
          expect[qwer2] <- act[good]
          # --------------------------------
          #if(summary(mod1)$r.squared < .9){
          if(mod$r2 < .9){
            i= tope #length(ladder.call) - maxo
          }
          # --------------------------------
        }else if(i > 5 ){
          #xx <- ladder.call[c(1:i)]
          #mod1 <- lm(expect[1:i]~poly(xx, degree=4, raw=TRUE))
          #expecto <- predict(mod1, data.frame(xx=ladder.call))
          ladder.call2 <- sort(ladder.call, decreasing = TRUE)
          expect2 <- sort(expect, decreasing = TRUE)
          xx <- cbind(ladder.call2[c(1:i)])#,ladder.call2[c(1:i)]^2,ladder.call2[c(1:i)]^3,ladder.call2[c(1:i)]^4)
          mod <- MSE3(mix=xx, miy=expect2)
          beta <- (mod)$beta
          if(length(which(is.na(beta))) > 0){
            beta[which(is.na(beta))] <- 0
          }
          toto <- cbind(matrix(ladder.call2))#,matrix(ladder.call2)^2, matrix(ladder.call2)^3,matrix(ladder.call2)^4)
          expecto <- cbind(rep(1,dim(toto)[1]),toto) %*% beta
          
          # remove missing data
          condo <- sort(expect[which(!is.na(expect))], decreasing = TRUE)
          # keep roxy values that were picked already
          act <- sort(roxy$pos[-which(roxy$pos %in% condo)], decreasing = TRUE)
          # see which is the next most likely value
          yoyo <- abs(expecto[i+1] - act)
          good <- which(yoyo == min(yoyo, na.rm = TRUE))
          qwer <- i
          qwer2 <- length(expect)-qwer
          expect[qwer2] <- act[good]
          # --------------------------------
          #if(summary(mod1)$r.squared < .9){
          if(is.na(mod$r2)){ # if model is too bad just assign a zero value
            mod$r2 <- 0.1
          }
          if(mod$r2 < .9){
            i= tope #length(ladder.call) - maxo
          }
          # --------------------------------
        }else if(i == tope & i == 3){
          expect[1:3] <- roxy$pos[www]
        }
        
        ##
      } # end of for loop
      #plot(x, type="l",ylim=c(0,800)); abline(v=expect, col="red", lty=3)
      ###########################################
      posi <- expect
      ## get rid of selected peaks after reaching the maximum values
      #tutu <- abs(length(x) - posi)
      #if(length(posi) > 3){
      #posi <- posi[1:which(tutu == min(tutu,na.rm = TRUE))]
      #}
      heii <- roxy$hei[which(roxy$pos %in% posi)]
      
      fact3 <- length(posi)/fact2
      
      if(length((posi)) < 6){ # good user
        fact <- summary(lm(ladder.call[1:length(posi)]~poly(posi, degree=length((posi))-1)))$r.squared * fact3
      }else{
        fact <- summary(lm(ladder.call[1:length(posi)]~poly(posi, degree=5)))$r.squared * fact3
      }
    
      #plot(ladder.call~posi)
      roxy2 <- list(pos=posi, hei=heii, wei=ladder.call[1:length(posi)], corr=abs(cor(ladder.call[1:length(posi)],posi)), error=fact)#sum(error, na.rm=TRUE))
    
      return(roxy2)
    }
    ############
    #----------- END OF CALLER
    ############
    ## end of caller
    # www <- data.frame(step4)[,1]
    #for(g in 1:dim(step4)[2]){
    #  plot(x,type="l",ylim=c(0,1000))
    #  abline(v=caller(www=step4[,g],roxy=roxy, ladder.call=ladder,x=x)$pos,col="red")
    #}
    
    rt <- apply(data.frame(step4), 2, FUN=caller, roxy=roxy, ladder.call=ladder,x=x)#
    corrs3 <- unlist(lapply(rt, function(x){x$error})) #; dis[which(dis == Inf)] <- 1
    #corrs <- unlist(lapply(rt, function(x){x[[4]]})) # extract correlations
    roxy3 <- rt[[which(corrs3 == max(corrs3))[1]]]
    #roxy3 <- rt[[5]]
    # plot(x, type="l",xlim=c(2000,7000), ylim=c(0,20000)); abline(v=roxy3$pos, col="red")
    if(draw == TRUE){
      limi <- sort(roxy3$hei, decreasing = TRUE)
      plot(x, type="l", xaxt="n", ylim=c(0,(limi[3]+1000)), cex.axis=0.6, las=2,  col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2, main=attributes(x)$mycomm, cex.main=cex.title)
      axis(1, at=roxy3$pos, labels=roxy3$wei, cex.axis=0.6)
      points(x=roxy3$pos, y=roxy3$hei,cex=1.1, col=transp("black",0.85))
      points(x=roxy3$pos, y=roxy3$hei, pch=20, col=transp("red",0.7))
      legend("topleft", legend=paste("Correlation:",round(roxy3$corr, digits=4), sep=""), bty="n")
      legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
      
    }
    roxy <- roxy3
    }
  }
  #####################
  return(roxy)
  #########
}
