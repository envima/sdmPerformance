
# modified from PresenceAbsence::cmx
cmx<- function (DATA, threshold = 0.5, which.model = 1, na.rm = FALSE) 
{
  if (is.logical(na.rm) == FALSE) {
    stop("'na.rm' must be of logical type")
  }
  if (sum(is.na(DATA)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(DATA), 1, sum)
      warning(length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      DATA <- DATA[NA.rows == 0, ]
    }
    else {
      return(NA)
    }
  }
  if (length(which.model) != 1) {
    stop("this function will only work for a single model, 'which.model' must be of length one")
  }
  if (which.model < 1 || round(which.model) != which.model) {
    stop("'which.model' must be a positive integer")
  }
  if (which.model + 2 > ncol(DATA)) {
    stop("'which.model' must not be greater than number of models in 'DATA'")
  }
  if (length(threshold) != 1) {
    stop("'threshold' must be a single number between zero and one")
  }
  if (max(threshold) > 1) {
    stop("'threshold' must be a single number between zero and one")
  }
  if (min(threshold) < 0) {
    stop("'threshold' must be a single number between zero and one")
  }
  # OBS.ind <- DATA[, 2] > 0
  OBS.ind<- DATA[, 2]
  if (threshold == 0) {
    PRED.ind = DATA[, which.model + 2] >= threshold
  }
  else {
    #PRED.ind = DATA[, which.model + 2] > threshold
    PRED.ind = ifelse(DATA[, which.model + 2] > threshold, 1, 0)
  }
  C = data.frame(predicted=PRED.ind,
                 observed=OBS.ind)
  #C = as.table(matrix(C, nrow = 2))
  #dimnames(C) = list(predicted = c(1, 0), observed = c(1, 0))
  #storage.mode(C) = "double"
  return(C)
}

# modiefied from mecofun::evalSDM: https://gitup.uni-potsdam.de/macroecology/mecofun/-/blob/master/R/evalSDM.R
confusionMatrix <- function(observation, predictions, thresh=NULL, thresh.method='MaxSens+Spec', req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1, weigths=rep(1, length(observation))){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  if (is.null(thresh)) {
    thresh.mat <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, req.sens=req.sens, req.spec = req.spec, FPC=FPC, FNC=FNC)
    thresh <- thresh.mat[thresh.mat$Method==thresh.method,2]
  }
  
  cmx.opt <- cmx(DATA= thresh.dat, threshold=thresh)
  return(cmx.opt)
}



# modiefied from mecofun::evalSDM: https://gitup.uni-potsdam.de/macroecology/mecofun/-/blob/master/R/evalSDM.R
confusionMatrix_PO <- function(observation, predictions, thresh=NULL, prediction){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  thresh.dat=thresh.dat%>%dplyr::filter(obs==1)
  thresh.dat.PO <- thresh.dat
  
  bg=terra::spatSample(prediction, size=10000, as.df=T,na.rm=T)[[1]]
  
  thresh.dat.bg <- data.frame(ID=seq_len(length(bg)), 
                              obs = 0,
                              pred = bg)
  thresh.dat<- rbind(thresh.dat, thresh.dat.bg);rm(thresh.dat.bg,bg)
  thresh.dat$ID <- 1:nrow(thresh.dat)
  
  thresh <- bestThresholdFpb(thresh.dat = thresh.dat)
  
  n_pres <- sum(thresh.dat$obs == 1)
  # downsample absences to match presences
  #set.seed(123)  # for reproducibility
  thresh.dat <- thresh.dat %>%
    group_by(obs) %>%
    slice_sample(n = n_pres, replace = FALSE) %>%
    ungroup()
  threshBalanced <- bestThresholdFpb(thresh.dat = thresh.dat)
  thresh <- 0.3*threshBalanced+thresh*0.7
  
  
  cmx.opt <- cmx(DATA= thresh.dat.PO, threshold=thresh)
  return(cmx.opt)
}

bestThresholdFpb <- function(thresh.dat,thresholds=NULL){
  # optimize thrshold with fbp
  thr <- if (is.null(thresholds)) sort(unique(thresh.dat$pred)) else sort(unique(thresholds))
  best_Fbp <- -Inf
  bestThreshold <- -Inf
  for (t in thr) {
    
    # fpb_value <- Fbp(actual=thresh.dat$obs, predicted=thresh.dat$pred)
    
    pred_bin <- as.integer(thresh.dat$pred >= t)
    TP <- sum(pred_bin == 1 & thresh.dat$obs == 1)
    #TN <- sum(pred_bin == 0 & thresh.dat$obs == 0)
    FP <- sum(pred_bin == 1 & thresh.dat$obs == 0)
    FN <- sum(pred_bin == 0 & thresh.dat$obs == 1)
    
    # precision <- TP / (TP + fp)
    #recall <- TP / (TP + FN)
    #precision2 = precision/(1-precision)
    if (TP == 0) {
      fbp<-0
    } else {
      fbp=(2*TP)/(TP+FN+FP)
      
    }
   # print(paste0("Fbp: ",fbp, ". threshold: ",t,"."))
    if (fbp > best_Fbp) {
      best_Fbp <- fbp
      bestThreshold <- t 
    }
  }
  return(bestThreshold)
  
}


# modiefied from mecofun::evalSDM: https://gitup.uni-potsdam.de/macroecology/mecofun/-/blob/master/R/evalSDM.R
confusionMatrix_PO <- function(observation, predictions, thresh=NULL, prediction){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  thresh.dat=thresh.dat%>%dplyr::filter(obs==1)
  thresh.dat.PO <- thresh.dat
  
  bg=terra::spatSample(prediction, size=10000, as.df=T,na.rm=T)[[1]]
  
  thresh.dat.bg <- data.frame(ID=seq_len(length(bg)), 
                              obs = 0,
                              pred = bg)
  thresh.dat<- rbind(thresh.dat, thresh.dat.bg);rm(thresh.dat.bg,bg)
  thresh.dat$ID <- 1:nrow(thresh.dat)
  
  thresh <- bestThresholdFpb(thresh.dat = thresh.dat)
  
  n_pres <- sum(thresh.dat$obs == 1)
  # downsample absences to match presences
  #set.seed(123)  # for reproducibility
  thresh.dat <- thresh.dat %>%
    group_by(obs) %>%
    slice_sample(n = n_pres, replace = FALSE) %>%
    ungroup()
  threshBalanced <- bestThresholdFpb(thresh.dat = thresh.dat)
  thresh <- 0.3*threshBalanced+thresh*0.7
  
  
  cmx.opt <- cmx(DATA= thresh.dat.PO, threshold=thresh)
  return(cmx.opt)
}

bestThresholdFpb <- function(thresh.dat,thresholds=NULL){
  # optimize thrshold with fbp
  thr <- if (is.null(thresholds)) sort(unique(thresh.dat$pred)) else sort(unique(thresholds))
  best_Fbp <- -Inf
  bestThreshold <- -Inf
  for (t in thr) {
    
    # fpb_value <- Fbp(actual=thresh.dat$obs, predicted=thresh.dat$pred)
    
    pred_bin <- as.integer(thresh.dat$pred >= t)
    TP <- sum(pred_bin == 1 & thresh.dat$obs == 1)
    #TN <- sum(pred_bin == 0 & thresh.dat$obs == 0)
    FP <- sum(pred_bin == 1 & thresh.dat$obs == 0)
    FN <- sum(pred_bin == 0 & thresh.dat$obs == 1)
    
    # precision <- TP / (TP + fp)
    #recall <- TP / (TP + FN)
    #precision2 = precision/(1-precision)
    if (TP == 0) {
      fbp<-0
    } else {
      fbp=(2*TP)/(TP+FN+FP)
      
    }
   # print(paste0("Fbp: ",fbp, ". threshold: ",t,"."))
    if (fbp > best_Fbp) {
      best_Fbp <- fbp
      bestThreshold <- t 
    }
  }
  return(bestThreshold)
  
}



# modiefied from mecofun::evalSDM: https://gitup.uni-potsdam.de/macroecology/mecofun/-/blob/master/R/evalSDM.R
confusionMatrix_PBG <- function(observation, predictions, thresh=NULL, prediction,  thresh.method='MaxSens+Spec', req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1, weigths=rep(1, length(observation))){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  thresh.dat=thresh.dat%>%dplyr::filter(obs==1)
  thresh.dat.PO <- thresh.dat
  
  bg=terra::spatSample(prediction, size=10000, as.df=T,na.rm=T)[[1]]
  
  thresh.dat.bg <- data.frame(ID=seq_len(length(bg)), 
                              obs = 0,
                              pred = bg)
  thresh.dat<- rbind(thresh.dat, thresh.dat.bg);rm(thresh.dat.bg,bg)
  thresh.dat$ID <- 1:nrow(thresh.dat)
  
  
  if (is.null(thresh)) {
    thresh.mat <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, req.sens=req.sens, req.spec = req.spec, FPC=FPC, FNC=FNC)
    thresh <- thresh.mat[thresh.mat$Method==thresh.method,2]
  }
  
  
  cmx.opt <- cmx(DATA= thresh.dat.PO, threshold=thresh)
  return(cmx.opt)
}

