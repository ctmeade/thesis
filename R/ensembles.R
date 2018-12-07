library(forecast)
library(xgboost)

nEnsemble <- function(ts, h){
  # ts is a time series object
  # h is the forecast horizon, the number of data points into the future to forecast
  
  # Fit the automatic models to the ts object
  # Arima Forecast
  aa <- forecast(auto.arima(ts, approximation = F, stepwise = F), h = h)$mean
  
  # Neural Network
  nn <- forecast(nnetar(ts, scale.inputs = T, decay = .2), h = h)$mean
  
  # Exponential Smoothing
  ets <- forecast(ets(ts), h = h)$mean
  
  # Theta forecast
  theta <- forecast(thetaf(ts, h = h), h = h)$mean
  
  # Combine each h by 1 forecast into a h by p matrix, where p is the number
  # of models fit to the data
  
  # Forecast Matrix
  fcMat <- as.matrix(cbind(as.numeric(aa), as.numeric(nn), as.numeric(ets), 
                           as.numeric(theta)))
  
  # Return the Forecast Matrix and the average forecast at each point on
  # the horizon as a list object
  out <- list(predictionMatrix = fcMat, forecast = rowMeans(fcMat))

  return(out)
}


xgStackEnsemble <- function(ts, h, validationSize = round(length(ts)/2)){
  # Demonstration of Ensemble Stacking Technique for Time series forecasts
  
  # Divide data into train and validation sets
  train <- head(ts, length(ts)-validationSize)
  valid <- tail(ts, validationSize)
  
  # Using training data, make prediction on validation set using nEnsemble
  valMat <- nEnsemble(train, validationSize)[[1]]
  cat("Validation Predications Made")
  
  # Use XGBoost algorithm to find best nonlinear transformation of forecasted validation
  # data points into actual validation observations
  xgStack <- xgboost(data = as.matrix(valMat), 
                     label = as.numeric(valid),
                     nrounds = 200,
                     eta = 0.1)
  cat("RF Stack completed")
  
  # Combine train and validation sets to make forecasts h steps into the horizon
  fcMat <- nEnsemble(ts, h)[[1]]
  cat("Forecast Matrix Created")
  # Apply the learned XGBoost nonlinear transformation
  xgPred <- predict(xgStack, fcMat)
  cat("Predictions Made")
  # Return Predictions
  return(xgPred)
}
