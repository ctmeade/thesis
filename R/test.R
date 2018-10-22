library(forecast)
library(Mcomp)
library(tidyverse)

ensemble <- function(ts, horizon){
  aa <- forecast(auto.arima(ts, approximation = F, stepwise = F), h = horizon)$mean
  nn <- forecast(nnetar(ts, scale.inputs = T, decay = .5), h = horizon)$mean
  ets <- forecast(ets(ts), h = horizon)$mean
  theta <- forecast(thetaf(ts), h = horizon)$mean
  snaive <- snaive(ts, h = horizon)$mean
  return(cbind(aa, nn, ets, theta, snaive))
}

coerce_stationary <- function(ts, method){
  if(method == "difference"){
    return(diff(ts))
  }
}

aaErrors <- c()
for(i in 1:200){
  cat(i)
  model <- auto.arima(M3[[i]]$x, approximation = F, stepwise = F)
  pred <- forecast(model, 6)$mean
  aaErrors[i] <- accuracy(pred, M3[[i]]$xx)[[5]]
}

nnErrors <- c()
for(i in 1:200){
  cat(i, "\n")
  pred <- forecast(nnetar(M3[[i]]$x, scale.inputs = T, decay = .5), h = 6)$mean
  nnErrors[i] <- accuracy(pred, M3[[i]]$xx)[[5]]
}

tErrors <- c()
for(i in 1:200){
  cat(i, "\n")
  pred <- forecast(thetaf(M3[[i]]$x), h = 6)$mean
  tErrors[i] <- accuracy(pred, M3[[i]]$xx)[[5]]
}

etsErrors <- c()
for(i in 1:200){
  cat(i, "\n")
  pred <- forecast(ets(M3[[i]]$x), h = 6)$mean
  etsErrors[i] <- accuracy(pred, M3[[i]]$xx)[[5]]
}

sErrors <- c()
for(i in 1:3000){
  cat(i, "\n")
  pred <- snaive(M3[[i]]$x, h = 6)$mean
  sErrors[i] <- accuracy(pred, M3[[i]]$xx)[[5]]
}

eErrors <- c()
for(i in 1:3000){
  cat(i)
  h = length(M3[[i]]$xx)
  model <- ensemble(M3[[i]]$x, horizon = h)
  eErrors[i] <- accuracy(rowMeans(model), M3[[i]]$xx)[[5]]
}

aaERR <- c()
nnERR <- c()
etsERR <- c()
thetaERR <- c()
snaiveERR <- c()
ensembleErrors <- c()
for(i in 1:length(M3)){
  cat(i, "\n")
  horizon = length(M3[[i]]$xx)
  
  aa <- forecast(auto.arima(M3[[i]]$x, approximation = F, stepwise = F), h = horizon)$mean
  aaERR[i] <- accuracy(aa, M3[[i]]$xx)[[5]]
  
  nn <- forecast(nnetar(M3[[i]]$x, scale.inputs = T, decay = .2), h = horizon)$mean
  nnERR[i] <- accuracy(nn, M3[[i]]$xx)[[5]]
  
  ets <- forecast(ets(M3[[i]]$x), h = horizon)$mean
  etsERR[i] <- accuracy(ets, M3[[i]]$xx)[[5]]
  
  theta <- forecast(thetaf(M3[[i]]$x), h = horizon)$mean
  thetaERR[i] <- accuracy(theta, M3[[i]]$xx)[[5]]
  
  snaive <- snaive(M3[[i]]$x, h = horizon)$mean
  snaiveERR[i] <- accuracy(snaive, M3[[i]]$xx)[[5]]
  
  ensemble <- rowMeans(cbind(aa, nn, ets, theta, snaive))
  ensembleErrors[i] <- accuracy(ensemble, M3[[i]]$xx)[[5]]
}
