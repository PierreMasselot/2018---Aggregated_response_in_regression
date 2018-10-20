# Regression with aggregated response
***
R functions implementing the methodology of time series regression when the response series is temporally aggregated, as proposed in the paper:

> Masselot, P., Chebana, F., Bélanger, D., St-Hilaire, A., Abdous, B., Gosselin, P., Ouarda, T.B.M.J., 2018. Aggregating the response in time series regression models, applied to weather-related cardiovascular mortality. *Science of The Total Environment* **628–629**, 217–225. [https://doi.org/10.1016/j.scitotenv.2018.02.014](https://doi.org/10.1016/j.scitotenv.2018.02.014)

## Description

In this paper, it is proposed to temporally aggregate the health response series in order to study low magnitude cumulative effects. A two-step methodology is proposed for this purpose:
1. Temporally aggregates the health response series, for instance through Nadaraya-Watson kernel smoothing;
2. Perform a time series regression model on the aggregated response, *i.e.* by modelling residuals through an ARIMA model.

The software is composed of three main functions (see examples below):
1. `NWsmooth`: aggregates a data series through Nadaraya-Watson kernel smoothing;
2. `tsRegression`: applies time series regression;
3. `predict.tsreg`: predict method for results from `tsRegression`. 
In addition, several functions for Kernels are provided for use with NWsmooth.

Description of input parameters and output of each function can be found in the header of each function in the file tsRegression.R

Dependency: R package `forecast`.

## Example

Simple example of a distributed lag nonlinear model (DLNM) with aggregated response.

```
library(forecast)
library(dlnm)
library(splines)

data(chicagoNMMAPS)

# Aggregating the Cardiovascular death series using the Epanechnikov kernel on future values with H = 7
ytilde <- NWsmooth(chicagoNMMAPS$cvd, Kernel = Kepanechnikov, h = 7, side = 3)

# Construct crossbasis for DLNM
Tcb <- crossbasis(chicagoNMMAPS$temp, lag = 21, argvar = list(knots = c(-4, 20, 24)), arglag = list(knots = c(3, 8, 13)))

# Include smooth time component in the model
mod.data <- data.frame(ytilde = log(ytilde$y), Tcb, time = ns(chicagoNMMAPS$time,df=8*14))

# Fit model
fit <- tsRegression(ytilde ~ ., data = mod.data, period = 365)

# Obtain fitted values
yhat <- exp(predict(fit))

# Plot the surface
Tcp <- crosspred(Tcb, coef = fit$coef[colnames(Tcb)], vcov = fit$var.coef[colnames(Tcb),colnames(Tcb)], model.link = "log")
plot(Tcp, ptype="contour", xlab="Temperature (°C)", ylab="Lag (Day)")
```