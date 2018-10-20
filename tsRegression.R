################################################################################
#
#            Regression with an aggregated response
#
#                Author : Pierre Masselot
#
#                     Date : 2017
#
#	         Dependencies : package 'forecast'
#
################################################################################

##############################################################
#
#  NWsmooth(x, Kernel, h, side = 2, standard = T)
#
#  Arguments
#      x: vector containing the time series to aggregate
#      Kernel: numeric vector, character or function. The weights to be used. See kernel functions below for Kernels
#      h: numeric. window size for the Kernel
#      side: 1 to use only past values, 2 to use a centered Kernel, 3 to use only future values.
#      standard: boolean. If TRUE (the default), the weights used to aggregate are standardized to unit sum.
#
#  Value
#      A list with the following elements
#          y: the aggregated series
#          wi: the weight vector used to aggregate x

NWsmooth <- function(x, Kernel, h, side = 2, standard = T)
{
    x <- as.vector(x)
    n <- length(x)
    if (inherits(Kernel,"numeric")){
       wi <- Kernel
    } else {
       if (inherits(Kernel,"character")){
          Kernel <- get(Kernel)
       }
       bd <- ifelse(side == 2, h / 2, h)
       wi <- Kernel(-n:n/bd)/bd       
    }
    wi <- wi[wi>0]
    nwi <- floor(length(wi)/2) 
    nullinds <- switch(side,1+(nwi+1):(2*nwi),NULL,1:nwi)
    wi[nullinds] <- 0     
    if (standard) wi <- wi/sum(wi)
    y <- filter(x,filter=wi,method="convolution",sides=2)
    return(list(y = y, wi = wi))
}

# Kernel functions for argument Kernel
Kuniform <- function(x){
    return(1/2 * (abs(x)<=1))         
}

Ktriangular <- function(x){
    return((1-abs(x)) * (abs(x)<=1))
}

Kepanechnikov <- function(x){
    return((3/4)*(1-x^2)*(abs(x)<=1))
}

Kbiweight <- function(x){
    return((abs(x)<=1)*(15/16)*(1-x^2)^2)
}

Ktriweight <- function(x){
    return((abs(x)<=1)*(35/32)*(1-x^2)^3)
}

Ktricube <- function(x){
    return((abs(x)<=1)*(70/81)*(1-abs(x)^3)^3)
}

Kgaussian <- function(x){
    return(dnorm(x))    
}

Kcosine <- function(x){
    return((pi/4)*cos(x*pi/2)*(abs(x)<=1))
}

Klogistic <- function(x){
    return(1/(exp(x)+2+exp(-x)))      
}

Ksigmoid <- function(x){
    return((2/pi)/(exp(x)+exp(-x)))
}

Ksilverman <- function(x){
    return(exp(-abs(x)/sqrt(2))*sin((abs(x)/sqrt(2))+(pi/4))/2)
}

Kmichels <- function(x){
    return((1+x-x^2-x^3)*(3/4)*(abs(x)<=1))
}

##############################################################
#
#  tsRegression(formula, data, arima.params = NULL, period = 1, ...)
#
#  Arguments
#      formula: the formula for the model. Used as in the 'lm' function. Can include a crossbasis for dlnm
#      data: data.frame containing the variables in the model
#      arima.params: vector of length 6. Orders c(p,d,q,P,D,Q) for the arima model (see function 'arima'). If NULL (the default), chosen automatically using auto.arima from the package forecast
#      period: numeric. The seasonality if a seasonal component for the arima model is expected (default to 1)
#      ... : other arguments to be passed to lm or auto.arima methods
#
#  Value
#      a list containing the following elements :
#          call: the model call
#          data: the variables used in the model
#          n: the number of observations used for fitting the model
#          lm.model: design matrix for linear model
#          Kalman.model: the state-space model fitted for arima
#          coef: fitted coefficients of the final model (includes both coefficients from the linear model and the arima of residuals)
#          order: the order of the arima fitted on residuals
#          yhat: fitted values
#          residuals: the residuals of the model
#          var.coef: covariance matrix of the coefficients
#          sigma2: the variance of residuals
#          loglik: log-likelihood of the final model
#          aic: AIC value of the final model          

tsRegression <- function(formula, data, arima.params = NULL, period = 1, ...)
{
   dots <- list(...)
   n <- nrow(data) 
   argsLm <- dots[names(dots) %in% names(formals(lm))]
   argsLm$formula <- formula
   argsLm$data <- data
   prel <- do.call("lm",argsLm)
   res <- residuals(prel)
   res <- ts(res,frequency=period)
   if (is.null(arima.params)){
      argsAutoArima <- dots[names(dots) %in% unique(c(names(formals(arima)),names(formals(auto.arima))))]
      argsAutoArima$y <- res
      res.fit <- do.call("auto.arima",argsAutoArima)
      arima.params <- res.fit$arma[c(1,6,2,3,7,4)]
   } else {
      stopifnot(length(arima.params)==6)
   }
   mf <- model.frame(formula=formula, data=data)
   mt <- attr(mf, "terms")
   argsArima <- dots[names(dots) %in% names(formals(arima))]
   argsArima$x <- model.response(mf)
   argsArima$xreg <- model.matrix(mt,mf)[,-1]
   argsArima$order <- arima.params[1:3]
   argsArima$seasonal <- list(order = arima.params[4:6], period = period)
   final.model <- do.call("arima",argsArima)
   if (!is.null(attributes(mf)$na.action)){
      yhat <- eps <- rep(NA,n)
      yhat[setdiff(1:n,attributes(mf)$na.action)-1] <- fitted(final.model)
      eps[setdiff(1:n,attributes(mf)$na.action)] <- residuals(final.model)
   } else {
      yhat <- c(NA,as.vector(fitted(final.model)))
      eps <- as.vector(residuals(final.model))
   }
   output <- list(call = match.call(), data = data, n = final.model$nobs, lm.model = mf, Kalman.model = final.model$model, coef = coef(final.model), order = arima.params, yhat = yhat, residuals = eps, var.coef = final.model$var.coef, sigma2 = final.model$sigma2, loglik = final.model$loglik, aic = final.model$aic)
   class(output) <- "tsreg"
   return(output)
}

##############################################################
#
#  predict.tsreg(object, newdata, na.action)
#
#  Arguments
#      object: output object from tsRegression
#      newdata: data.frame containing new values for prediction. If missing, calibration data are used.
#      na.action: vector of length 6. Orders c(p,d,q,P,D,Q) for the arima model (see function 'arima'). If NULL (the default), chosen automatically using auto.arima from the package forecast
#      period: how to deal with missing values. Default to na.remove
#
#  Value
#      a vector containing the predicted time series

predict.tsreg <- function(object, newdata, na.action)
{
   stopifnot(inherits(object, "tsreg"))
   mt <- attr(object$lm.model, "terms")
   if (missing(newdata) || is.null(newdata)){
      #newX <- model.matrix(mt,object$lm.model)[,-1]
      pred <- object$yhat
   } else {
      Terms <- delete.response(mt)
      m <- model.frame(Terms, newdata, na.action = na.action)
      newX <- model.matrix(Terms, m)
      n <- nrow(newX)
      narma <- sum(object$order[c(1,3,4,6)])
      if (narma > 0){
         xcoef <- object$coef[-(1:narma)]
      } else {
         xcoef <- object$coef
      }      
      xpred <- newX %*% xcoef
      epred <- KalmanForecast(n, object$Kalman.model)
      pred <- epred$pred + xpred 
   }
   return(as.vector(pred))
}
