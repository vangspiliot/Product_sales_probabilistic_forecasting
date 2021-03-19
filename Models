pl <- function(x, y, yhat, q){
  # Computes scaled pinball loss for a given quantile q
  #x, y, and yhat correspond to the insample, outsample, and forecasts, respectively
  x<-as.numeric(x) ; y<-as.numeric(y) ; yhat<-as.numeric(yhat)
  
  temp_score <- 0
  for (k in 1:length(y)){
    if (y[k]>=yhat[k]){
      temp_score <- temp_score + q*(y[k]-yhat[k])
    }else{
      temp_score <- temp_score + (1-q)*(yhat[k]-y[k])
    }
  }
  temp_score <- (temp_score/length(y))/mean(abs(diff(x)))
  return(temp_score)
}

ses_PL <- function(x, h, q, groups){
  # PL optimization for ses_quantile()
  error <- c()
  for (a in seq(0.01,0.99,0.01)){
    model <- as.numeric(ses(x, h=1, alpha = a)$fitted)
    fitted = data_fitted <- c()
    for (i in 1:nrow(groups)){
      data_fitted <- c(data_fitted, groups[i,])
      fitted <- c(fitted, rep(model[i], h))
    }
    
    temp_score <- 0
    for (k in 1:length(data_fitted)){
      if (data_fitted[k]>=fitted[k]){
        temp_score <- temp_score + q*(data_fitted[k]-fitted[k])
      }else{
        temp_score <- temp_score + (1-q)*(fitted[k]-data_fitted[k])
      }
    }
    temp_score <- temp_score/length(data_fitted)
    error <- c(error, temp_score)
  }
  a <- seq(0.01,0.99,0.01)[which.min(error)]
  forecast <- rep(ses(x, h=1, alpha = a)$mean,h)
  return(forecast)
}

ses_quantile <- function(x, q, fh, opt, positive){
  
  # Applies SES on the empirical quantiles of the series instead of the original data
  
  # Create windows of size equal to fh
  x <- as.numeric(x)
  cut <- fh
  N <- floor(length(x)/cut)
  insampleCut <- tail(x, N*cut)
  insample_group <- NULL
  q_out <- matrix(NA, nrow = N, ncol = length(q))
  for (w_id in 1:N){
    temp <- insampleCut[(cut*(w_id-1)+1):(cut*w_id)]
    quants <- as.numeric(quantile(temp,q))
    for (q_id in 1:length(q)){
      q_out[w_id,q_id] <- as.numeric(quants[q_id]) #Series containing quantiles per window
    }
    insample_group <- rbind(insample_group, temp) #windows of data created
  }
  
  forc <- matrix(NA, nrow = fh, ncol = length(q))
  if (opt=="mse"){
    # MSE optimization
    for (q_id in 1:length(q)){
      forc[,q_id] <- as.numeric(ses(q_out[,q_id], h=fh)$mean)
    }
    method <- "SES-QEE-mse"
  }else{
    # PL optimization
    for (q_id in 1:length(q)){
      forc[,q_id] <- ses_PL(q_out[,q_id], h = fh, q=q[q_id], insample_group)
    }
    method <- "SES-QEE-pl"
  }
  forc <- data.frame(forc)
  colnames(forc) <- paste0("q",q)
  
  if (positive==T){
    for (qid in 1:ncol(forc)){
      for (fhid in 1:nrow(forc)){
        forc[fhid,qid] <- max(0,forc[fhid,qid])
      }
    }
  }
  
  return(forc)
  
}

res_sim_errors <- function(x, q, fh, Nboots, positive, type_m){
  
  # Computes SES with empirically computed residual errors (SES-emp) as well as
  #SES with simulated overlapping and non-overlapping forecast errors (SES-sim-o, SES-sim-no, SES-sim-no-fh, SES-sim-no-fh)
  
  x <- as.numeric(x)
  # Compute SES point forecasts
  ses_for <- as.numeric(ses(x, h = fh)$mean)
  forc <- matrix(NA, nrow = fh, ncol = length(q))
  
  if (type_m=="SES-emp"){
    
    # Empirical residual errors
    in_errors_ses <- x-as.numeric(ses(x, h = fh)$fitted)
    sq <- quantile(unlist(in_errors_ses), q)
    for (q_id in 1:length(q)){
      forc[,q_id] <- ses_for+sq[q_id]
    }
    
  }else{
    
    if ((type_m=="SES-sim-o")|(type_m=="SES-sim-o-fh")){
      # Create overlapping windows
      boots <- NULL
      for (i in 1:Nboots){
        if ((length(x)-fh-(i-1))>=7){
          bin <- x[1:c(length(x)-fh-(i-1))]
          bout <- x[c(length(x)-fh-(i-1)+1):c(length(x)-(i-1))]
          boots[length(boots)+1] <- list(list(x=bin, xx=bout))
        }
      }
    }else {
      # Create non-overlapping windows
      boots <- NULL
      for (i in 1:Nboots){
        if ((length(x)-fh-(i-1)*fh)>=7){
          bin <- x[1:c(length(x)-fh-(i-1)*fh)]
          bout <- x[c(length(x)-fh-(i-1)*fh+1):c(length(x)-(i-1)*fh)]
          boots[length(boots)+1] <- list(list(x=bin, xx=bout))
        }
      }
    }
    
    # Estimate forecast errors for all overlapping windows
    res_ses <- NULL
    for (i in 1:length(boots)){
      sesf <- as.numeric(ses(boots[[i]]$x, h = fh)$mean)
      res_ses <- rbind(res_ses, boots[[i]]$xx-sesf)
    }
    
    if ((type_m=="SES-sim-o") | (type_m=="SES-sim-no")){
      # Simulated forecast errors 
      sq <- quantile(unlist(res_ses), q)
      for (q_id in 1:length(q)){
        forc[,q_id] <- ses_for+sq[q_id]
      }
    }else{
      sq <- lapply(c(1:fh), function(x) quantile(res_ses[,x], q))
      for (q_id in 1:length(q)){
        forc[,q_id] <- unlist(lapply(c(1:fh), function(x) ses_for[x]+sq[[x]][q_id]))
      }
    }
    
  }
  
  forc <- data.frame(forc)
  colnames(forc) <- paste0("q",q)
  
  if (positive==T){
    for (qid in 1:ncol(forc)){
      for (fhid in 1:nrow(forc)){
        forc[fhid,qid] <- max(0,forc[fhid,qid])
      }
    }
  }
  
  return(forc)
}

m5_stat_benchmarks <- function(x, q, fh, positive, model){
  
  # These are the five official statistical benchmarks used in the M5 competition
  
  PIs <- unique((abs(1-q_targets*2)*100))
  PIs <- PIs[PIs>0]
  PIs <- PIs[order(-PIs)]
  
  if (model=="ses"){
    mf <- ses(x, h = fh, level = PIs)
  }else if (model=="naive"){
    mf <- naive(x, h = fh, level = PIs)
  }else if (model=="snaive"){
    mf <- snaive(x, h = fh, level = PIs)
  }else if (model=="ets"){
    mf <- forecast(ets(x), h = fh, level = PIs)
  }else{
    mf <- forecast(auto.arima(x), h = fh, level = PIs)
  }
  
  
  forc <- matrix(NA, nrow = fh, ncol = length(PIs)*2+1)
  
  if ((model!="naive")&(model!="snaive")){
    for (qid in 1:length(PIs)){
      forc[,length(PIs)-qid+1] <- as.numeric(mf$lower[,qid])
    }
    forc[,qid+1] <- as.numeric(mf$mean)
    for (qid in 1:length(PIs)){
      forc[,length(PIs)+1+qid] <- as.numeric(mf$upper[,qid])
    }
  }else{
    for (qid in 1:length(PIs)){
      forc[,qid] <- as.numeric(mf$lower[,qid])
    }
    forc[,qid+1] <- as.numeric(mf$mean)
    for (qid in 1:length(PIs)){
      forc[,(2*length(PIs)+1)-qid+1] <- as.numeric(mf$upper[,qid])
    }
  }
  
  
  forc <- data.frame(forc)
  colnames(forc) <- paste0("q",c(((100-PIs)/100)/2,0.5,rev(1-((100-PIs)/100)/2)))
  forc <- forc[,paste0("q",q)]
  
  if (positive==T){
    for (qid in 1:ncol(forc)){
      for (fhid in 1:nrow(forc)){
        forc[fhid,qid] <- max(0,forc[fhid,qid])
      }
    }
  }
  
  return(forc)
}
dist_calc <- function(x, q, fh, positive, distT){
  
  #Probabilistic forecasts based on empirical and theoretical data distributions
  
  # Quantile empirical estimation (QEE) corresponds to the Kernel benchmark of the competition
  #which was originally computed using type=7 (R's default) in the quantile function (estimates the mode)
  #However, here we also compute QEE for type=8 as quantile estimates are approximately median-unbiased 
  #regardless of the distribution of x
  x <- as.numeric(x)
  if (distT=="norm"){
    quants <- qnorm(q, mean = mean(x), sd = sd(x))
  }else if (distT=="pois"){
    quants <- qpois(q, lambda = mean(x))
  }else if (distT=="nbinom"){
    
    fit <- try(fitdist(x, "nbinom"))
    if("try-error" %in% class(fit)){
      quants <- qpois(q, lambda = mean(x))
    }else{
      quants <- qnbinom(q, size = fit$estimate[1], mu=fit$estimate[2])
    }
  }else if (distT=="emp-8"){
    quants <- as.numeric(quantile(x,q,type=8))
  }else{
    quants <- as.numeric(quantile(x,q,type=7))
  }
  
  forc <- data.frame(matrix(NA, nrow = fh, ncol = length(q)))
  colnames(forc) <- paste0("q",q)
  for (qid in 1:ncol(forc)){
    forc[,qid] <- quants[qid]
  }
  
  if (positive==T){
    for (qid in 1:ncol(forc)){
      for (fhid in 1:nrow(forc)){
        forc[fhid,qid] <- max(0,forc[fhid,qid])
      }
    }
  }
  
  return(forc)
  
}
