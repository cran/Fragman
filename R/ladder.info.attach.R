ladder.info.attach <- function(stored, ladder, channel.ladder=NULL, ci.upp=1.96, ci.low=1.96, dev=50, warn= FALSE, method="red", ladd.init.thresh=200, env = parent.frame(), prog=TRUE, draw=TRUE){
  
  if(is.null(channel.ladder)){
    channel.ladder <- dim(stored[[1]])[2]
  }else{channel.ladder <- channel.ladder}

  #####################
  ### this part extracts all the models for each single plant in my plants
  list.ladders <- lapply(stored, function(x){y <- x[,channel.ladder]; return(y)})
  # extract ladder channels for all plants
  if(prog==TRUE){
    res <- lapply_pb(list.ladders, find.ladder, ladder=ladder, ci.upp=ci.upp, ci.low=ci.low, dev=dev, warn=warn, method=method,init.thresh=ladd.init.thresh, draw=draw)
  }
  if(prog==FALSE){
    res <- lapply(list.ladders, find.ladder, ladder=ladder, ci.upp=ci.upp, ci.low=ci.low, dev=dev, warn=warn,  method=method,init.thresh=ladd.init.thresh, draw=draw)
  }
  env$list.data.covarrubias <- res
  ######################
}
