library(readxl)
library(astsa)
library(tseries)

#Data Cleaning
# spy data 
spy_raw <- read_excel("SPY_2020-2024.xlsx")
spy_ts <- ts(spy_raw$Close, start = c(2019, 8), frequency = 12)

# hpi data
hpi_raw <- read.csv("Case-Shiller Home Price Index.csv")
hpi_ts <- ts(hpi_raw$CSUSHPINSA[-1], start = c(2019, 8), frequency = 12)


#EDA
par(mfrow = c(2, 1))
plot(spy_ts, main = "SPY Price", xlab = "Year", ylab = "Price ($)", col = "blue", lwd = 2)
plot(hpi_ts, main = "House Price Index", xlab = "Year", ylab = "Index Jan2000=100", col = "red", lwd = 2)
par(mfrow = c(1, 1))

ccf2(spy_ts, hpi_ts)

library(fGarch)

plot(diff(new_spy))
acf2(diff(new_spy)) # -> white noise

dif_spy <- diff(diff(log(new_spy)))
dif_hpi <- diff(hpi_ts)
acf2(dif_spy, max.lag = 50)
sarima(new_spy, p=0, d=2, q=1)

#fit an ARCH model
diff_spy_centered <- dif_spy - mean(dif_spy)
acf2(diff_spy_centered)  # -> white noise
acf2(diff_spy_centered^2, max.lag = 50)  # -> seems to be AR(1)

ARCH_1_spy <- garchFit(~garch(1,2), new_spy, cond.dist = "std", trace=FALSE)
summary(ARCH_1_spy) #<- #fail to pass Shapiro test
plot(ARCH_1_spy, which=13) #<- #not really normal distribution 

ARCH_11_spy <- garchFit(~garch(1,1), dif_spy, trace=FALSE)
summary(ARCH_11_spy)  # <- Garch(1,1) is not better

ARCH_2_spy <- garchFit(~garch(2,0), dif_spy, trace=FALSE)
summary(ARCH_2_spy)  #<- not better, keep the garch(1,0)


spy_model <- garchFit(~garch(1,0), dif_spy, cond.dist = "std")
summary(spy_model)
plot(spy_model, which=13)

# fitted value
plot(dif_spy, type = "l", main = "Time Series Plot of Dif SPY")
time_difspy <- c(time(dif_spy))
fitted_value <- spy_model@fitted
lines(time_difspy,fitted_value, col="red")

# the conditional standard deviation
plot(dif_spy, type = "l", main = "Time Series Plot of Dif SPY")
fitted_sd <- spy_model@sigma.t
lines(time_difspy, fitted_sd, col="blue")
plot(spy_model, which = 3)

##############Prediction HERE###############
predict(ARCH_1_spy_std, n.ahead = 5)

#     meanForecast meanError standardDeviation
#1     7.323045  15.96550          15.96550
#2     7.323045  19.69929          19.69929
#3     7.323045  22.30748          22.30748
#4     7.323045  24.24388          24.24388
#5     7.323045  25.72801          25.72801

#using 95% confidence interval, next 5 value

#     Lower Bound	Upper Bound
#1		-23.976335	38.622425
#2		-31.388755	45.034845
#3		-36.399605	51.045695
#4  	-40.186955	54.833045
#5 		-43.094855	57.740945

#Aug 2024 SPY is 541.3495
#Garch model show that the monthly average of Dif_Spy is 7.323045
#The mean value of Sep 2024 SPY is 548.672545
#And the range of Sep 2024 SPY is between 517.373165 and 579.971925

#The mean value of Oct 2024 SPY is 555.99559
#and so on

ar1_fit <- sarima(new_spy, 1, 0, 0)
res <- ar1_fit$fit$residuals
plot(res)
acf2(res)
acf2(res^2)

#############################################################################

#Building the regression
ccf2(dif_spy, dif_hpi)  #>- we need pre-whitening

#pre-whiten using large AR model
acf2(dif_spy)
sarima(dif_spy, 5,0,0)
spy_AR5 <- sarima(dif_spy, 5,0,0, no.constant = T)
spy_AR5_residuals <- spy_AR5$fit$residuals
spy_AR5$ttable
hpi_filtered <-filter(dif_hpi, filter = c(1, -spy_AR5$ttable[, 1]), sides=1)
ccf2(spy_AR5_residuals, hpi_filtered, na.action=na.omit)
# lag 2 is significant, spy leads hpi by two monthes ahead

#fit the regression model
df <- ts.intersect(y = dif_hpi, 
                   xlag2 = lag(dif_spy,-2),
                   dframe = TRUE)
lm_model <- lm(y ~ xlag2, data=df, na.action = na.omit)
summary(lm_model)
acf2(resid(lm_model))
sarima(df$y, p = 1, d = 0, q = 0, xreg = df$xlag2)
# the diagnostic test looks weird, and the p-value of intercept and xlag2 are too large
# maybe the error terms do not follow ARIMA

#Coefficients:
#                 Estimate   Std. Error t value Pr(>|t|)    
#(Intercept)    1.86419    0.32134   5.801 3.54e-07 ***
#  xlag2        0.03664    0.01594   2.299   0.0254 *  

# formula: dif_hpi = 1.86419 + 0.03664*dif_spy(t-2) + error


#predicting
spy_future = data.frame(xlag2 = c(19.1784668, -2.3228760, 7.323045, 7.323045, 7.323045))
res <- predict(lm_model, newdata = spy_future)
print(res)
# 1        2        3        4        5 
# 2.566818 1.779089 2.132480 2.132480 2.132480 





##### Build ARIMAX Model #####

install.packages("TSA")
library(TSA)

ccf2(dif_spy, dif_hpi)
prewhiten(dif_spy, dif_hpi) # -2, -3, -4 are significant
prewhiten(spy_ts, hpi_ts)
curr_df <- ts.intersect(y = dif_hpi,
                        xlag2 = lag(dif_spy, -2),
                        xlag3 = lag(dif_spy, -3),
                        xlag4 = lag(dif_spy, -4),
                        dframe = TRUE)
new_model <- lm(y ~ xlag2 + xlag3 + xlag4, data=curr_df, na.action = na.omit)

summary(new_model) # all sig
acf(resid(new_model))
pacf(resid(new_model))
acf2(resid(new_model), max.lag = 50)

sarima(curr_df$y, p = 2, d = 0, q = 0,
       xreg = cbind(curr_df$xlag2, curr_df$xlag3, curr_df$xlag4))
#passes diagnostics




##### Build SARIMA Model #####

install.packages("TSA")
library(TSA)

ccf2(dif_spy, dif_hpi)
prewhiten(dif_spy, dif_hpi) # -2, -3, -4 are significant

curr_df <- ts.intersect(y = dif_hpi,
                        xlag2 = lag(dif_spy, -2),
                        xlag3 = lag(dif_spy, -3),
                        xlag4 = lag(dif_spy, -4),
                        dframe = TRUE)
new_model <- lm(y ~ xlag2 + xlag3 + xlag4, data=curr_df, na.action = na.omit)

summary(new_model) # all sig
acf(resid(new_model))
pacf(resid(new_model))
acf2(resid(new_model), max.lag = 50)
library(forecast)
auto.arima(resid(new_model), seasonal = TRUE)
sarima(curr_df$y, p = 2, d = 0, q = 0,
                  P = 1, D = 0, Q = 0, S = 12,
                  xreg = cbind(curr_df$xlag2, curr_df$xlag3, curr_df$xlag4), no.constant = T)
#passes diagnostics





##### USING LOG WITH DIFFERENCING #####

logdifspy <- diff(log(spy_ts))
logdifhpi <- diff(log(hpi_ts))


ccf2(logdifspy, logdifhpi)
prewhiten(logdifspy, logdifhpi)

curr_df <- ts.intersect(y = logdifhpi,
                        ylag = lag(logdifhpi, -1),
                        ylag2 = lag(logdifhpi, -2),
                        xlag1 = lag(logdifspy, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + ylag2 + xlag1, data=curr_df, na.action = na.omit)
summary(new_model)
par(mfrow = c(2, 1))
acf(resid(new_model)) 
pacf(resid(new_model))
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=0,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))



##### USING LOG NO DIFFERENCING #####

loghpi <- log(hpi_ts)
logspy <- log(spy_ts)

curr_df <- ts.intersect(y = loghpi,
                        ylag = lag(loghpi, -1),
                        ylag2 = lag(loghpi, -2),
                        xlag1 = lag(logspy, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + ylag2 + xlag1, data=curr_df, na.action = na.omit)
summary(new_model) # all sig
par(mfrow = c(2, 1))
acf(resid(new_model)) # looks like ar(2)
pacf(resid(new_model)) # looks like ar(2)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))




loghpi <- log(hpi_ts)
logspy <- log(spy_ts)

curr_df <- ts.intersect(y = loghpi,
                        ylag1 = lag(loghpi, -1),
                        xlag3 = lag(logspy, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag1 + xlag3, data=curr_df, na.action = na.omit)
summary(new_model) # all sig
par(mfrow = c(2, 1))
acf(resid(new_model)) # looks like ar(2)
pacf(resid(new_model)) # looks like ar(2)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))

model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag1, curr_df$xlag3))




curr_df <- ts.intersect(y = hpi_ts,
                        ylag = lag(hpi_ts, -1),
                        ylag2 = lag(hpi_ts, -2),
                        xlag1 = lag(spy_ts, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + ylag2 + xlag1, data=curr_df, na.action = na.omit)
summary(new_model) # all sig
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))


curr_df <- ts.intersect(y = loghpi,
                        ylag = lag(loghpi, -1),
                        ylag2 = lag(loghpi, -2),
                        xlag1 = lag(logspy, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + ylag2 + xlag1, data=curr_df, na.action = na.omit)
summary(new_model)
par(mfrow = c(2, 1))
acf(resid(new_model)) 
pacf(resid(new_model))
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=0,d=1,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))


curr_df <- ts.intersect(y = hpi_ts,
                        ylag = lag(hpi_ts, -1),
                        ylag2 = lag(hpi_ts, -2),
                        xlag1 = lag(spy_ts, -1),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + ylag2 + xlag1, data=curr_df, na.action = na.omit)
summary(new_model)
par(mfrow = c(2, 1))
acf(resid(new_model)) 
pacf(resid(new_model))
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=0,d=1,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))


curr_df <- ts.intersect(y = hpi_ts,
                        ylag = lag(hpi_ts, -1),
                        xlag3 = lag(spy_ts, -3),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + xlag3, data=curr_df, na.action = na.omit)
summary(new_model)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$xlag3))


curr_df <- ts.intersect(y = dif_hpi,
                        ylag = lag(dif_spy, -1),
                        xlag2 = lag(dif_spy, -2),
                        xlag3 = lag(dif_spy, -3),
                        xlag4 = lag(dif_spy, -4),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + xlag2 + xlag3 + xlag4, data=curr_df, na.action = na.omit)
summary(new_model)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$xlag2, curr_df$xlag3, curr_df$xlag4))


new_model <- lm(y ~ xlag3, data=curr_df, na.action = na.omit)
summary(new_model)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=2,D=1,Q=1,S=12, xreg = curr_df$xlag3)
library(forecast)
auto.arima(resid(new_model))

log_spy <- diff(log(spy_ts))
log_hpi <- diff(diff(log(hpi_ts), 12))
plot(spy_ts)
plot(diff(spy_ts))
plot(log_spy)
plot(l)

plot(hpi_ts)
plot(log(hpi_ts))
plot(diff(hpi_ts, 12))

ccf2(diff(spy_ts), diff(hpi_ts, 12))
prewhiten(diff(spy_ts), diff(hpi_ts, 12))
acf2(diff(spy_ts))
model <- sarima(diff(spy_ts), p=8,d=0,q=0, no.constant = T)
x_residuals <- model$fit$residuals
y_filtered <- filter(diff(hpi_ts, 12), filter = c(1, -model$ttable[, 1]), sides=1)
ccf2(x_residuals, y_filtered)
curr_df <- ts.intersect(y = diff(hpi_ts, 12),
                        xlag2 = lag(diff(spy_ts), 2),
                        xlag3 = lag(diff(spy_ts), 3),
                        dframe = TRUE)
aligned_data <- na.omit(ts.intersect(
  y = diff(hpi_ts, 12),
  xlag2 = lag(diff(spy_ts), 2),
  xlag3 = lag(diff(spy_ts), 3)
))
model <- lm(y ~ xlag2 + xlag3, data=aligned_data, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 40)
sarima(aligned_data[, "y"], p=3,d=0,q=1, P=3,D=0,Q=0,S=12, xreg = cbind(aligned_data[, "xlag2"], aligned_data[, "xlag3"]))

plot(log_hpi)
acf2(log_hpi)
ccf2(log_spy, log_hpi)
acf2(log_spy, max.lag = 50)
prewhiten(log_spy, log_hpi)
model <- sarima(log_spy, p=8,d=0,q=0, P=0,D=0,Q=0,S=12)
acf2(model$fit$residuals)
x_resids <- model$fit$residuals
y_filtered <- filter(log_hpi, filter = c(1, -model$ttable[, 1]), sides=1)
ccf2(x_resids, y_filtered)

curr_df <- ts.intersect(y = log_hpi,
                        ylag = lag(log_hpi, -1),
                        xlag1 = lag(log_spy, -1),
                        xlag2 = lag(log_spy, -2),
                        xlag3 = lag(log_spy, -3),
                        xlag4 = lag(log_spy, -4),
                        dframe = TRUE)
new_model <- lm(y ~ xlag1 + xlag2 + xlag3, data=curr_df, na.action = na.omit)
summary(new_model)
acf2(resid(new_model), max.lag = 50)
plot(resid(new_model))
model <- sarima(curr_df$y, p=1,d=0,q=0, P=1,D=1,Q=1,S=12, xreg = curr_df$xlag3)


acf2(diff(hpi_ts))
sarima(hpi_ts, p=3,d=1,q=0, P=1,D=0,Q=0,S=12)
sarima.for(hpi_ts, n.ahead = 3, p=3,d=1,q=0, P=1,D=0,Q=0,S=12)
triplees <- HoltWinters(hpi_ts)
pred <- forecast(triplees, 3)
plot(pred)

# Create newxreg with constant values (last observed values)
newxreg <- cbind(
  ylag = tail(curr_df$ylag, 1),    # Use last observed ylag
  xlag2 = tail(curr_df$xlag2, 1),  # Use last observed xlag2
  xlag3 = tail(curr_df$xlag3, 1),  # Use last observed xlag3
  xlag4 = tail(curr_df$xlag4, 1)   # Use last observed xlag4
)

# Repeat for n.ahead rows
newxreg <- matrix(rep(newxreg, each = 3), nrow = 3, byrow = TRUE)

sarima.for(curr_df$y, n.ahead = 1, p=2,d=0,q=0, P=1,D=0,Q=0,S=12,
           xreg = cbind(curr_df$ylag, curr_df$xlag2, curr_df$xlag3, curr_df$xlag4),
           newxreg = newxreg)

sarima(spy_ts, p=1,d=0,q=0, P=0,D=0,Q=0,S=12)
hpi_ts

# Example time series
n <- length(hpi_ts)
train_data <- ts(hpi_ts[1:(n - 5)], start = c(2019, 8), frequency = 12)  # Training data (excluding the last 5 points)
test_data <- ts(hpi_ts[(n - 4):n], start = c(2024, 3), frequency = 12)  # Last 5 data points (to evaluate the models)

model1 <- sarima(train_data, p=2,d=1,q=0, P=1,D=0,Q=0,S=12)
model2 <- HoltWinters(train_data)
forecast1 <- sarima.for(train_data, n.ahead = 5, p=2,d=1,q=0, P=1,D=0,Q=0,S=12)$pred
forecast2 <- forecast(model2, 5)
forecast1
forecast2
forecast2 <- ts(c(319.1786, 320.3689, 321.4619, 322.8275, 324.7552), start = c(2024, 3), frequency = 12)
forecast2
plot(forecast2)
forecast3 <- c(320.5539, 323.5403, 326.8730, 329.3700, 331.6809)
forecasts <- cbind(forecast1, forecast2, forecast3)
forecasts
# Compute error metrics
mae <- colMeans(abs(forecasts - test_data))  # Mean Absolute Error
rmse <- sqrt(colMeans((forecasts - test_data)^2))  # Root Mean Squared Error

# Print results
results <- data.frame(Model = c("Model 1", "Model 2", "Model 3"), MAE = mae, RMSE = rmse)
print(results)

acf2(dif_spy, max.lag = 50)
sarima(spy_ts, p=0,d=1,q=0, P=0,D=0,Q=0,S=12)
sarima(spy_ts, p=0,d=1,q=0)

tail(hpi_ts,5)

spy_train <- ts(spy_ts[1:(n - 5)], start = c(2019, 8), frequency = 12)

curr_df <- ts.intersect(y = diff(train_data),
                        ylag = lag(diff(spy_train), -1),
                        xlag2 = lag(diff(spy_train), -2),
                        xlag3 = lag(diff(spy_train), -3),
                        xlag4 = lag(diff(spy_train), -4),
                        dframe = TRUE)
new_model <- lm(y ~ ylag + xlag2 + xlag3 + xlag4, data=curr_df, na.action = na.omit)
summary(new_model)
acf2(resid(new_model), max.lag = 45)
plot(resid(new_model))
model <- sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$xlag2, curr_df$xlag3, curr_df$xlag4))

forecast_ylag <- c(tail(curr_df$ylag, 1), predict(arima(spy_train, order = c(0, 1, 0)), n.ahead = 4)$pred)
forecast_xlag2 <- c(tail(curr_df$xlag2, 2), predict(arima(spy_train, order = c(0, 1, 0)), n.ahead = 3)$pred)
forecast_xlag3 <- c(tail(curr_df$xlag3, 3), predict(arima(spy_train, order = c(0, 1, 0)), n.ahead = 2)$pred)
forecast_xlag4 <- c(tail(curr_df$xlag4, 4), predict(arima(spy_train, order = c(0, 1, 0)), n.ahead = 1)$pred)
forecast_xlag4
curr_df$xlag4
newxreg <- cbind(
  ylag = forecast_ylag,
  xlag2 = forecast_xlag2,
  xlag3 = forecast_xlag3,
  xlag4 = forecast_xlag4
)

forecast3 <- sarima.for(curr_df$y, n.ahead = 5,
                        p=2,d=0,q=0, P=1,D=0,Q=0,S=12,
                        xreg = cbind(curr_df$ylag, curr_df$xlag2, curr_df$xlag3, curr_df$xlag4),
                        newxreg = newxreg)
forecast3
hpi_ts
last_val <- ts(rep(316.972,5))
undiff <- last_val + forecast3




curr_df <- ts.intersect(y = hpi_ts,
                        ylag = lag(hpi_ts, -1),
                        ylag2 = lag(hpi_ts, -2),
                        xlag1 = lag(spy_ts, -1),
                        xlag2 = lag(spy_ts, -2),
                        xlag3 = lag(spy_ts, -3),
                        xlag4 = lag(spy_ts, -4),
                        xlag5 = lag(spy_ts, -5),
                        xlag6 = lag(spy_ts, -6),
                        dframe = TRUE)
model <- lm(y ~ ylag + ylag2 + xlag3, data=curr_df, na.action = na.omit)
summary(model)
acf2(diff(diff(resid(model)), 12), max.lag = 40)
plot(diff(resid(model)))
auto.arima(resid(model))
sarima(curr_df$y, p=2,d=0,q=0, P=1,D=0,Q=0,S=12, xreg = cbind(curr_df$ylag, curr_df$ylag2, curr_df$xlag1))

hpi_ts
new_spy_raw <- read.csv('spy_monthly.csv')
new_spy <- ts(new_spy_raw$Close, start = c(2019, 8), frequency = 12)
new_spy

par(mfrow = c(1,1))
ccf2(diff(new_spy), diff(hpi_ts))
prewhiten(diff(new_spy), diff(hpi_ts))
acf2(diff(diff(new_spy, 12)), max.lag = 40)
model <- sarima(new_spy, p=0,d=1,q=0, P=0,D=0,Q=0,S=12, no.constant = T)
new_res <- model$fit$residuals
acf2(new_res, max.lag = 50)
y_filter = filter(diff(hpi_ts), filter = c(1, -model$ttable[, 1]), sides=1)
ccf2(new_res, y_filter)
df <- ts.intersect(y = diff(hpi_ts),
                        ylag = lag(diff(hpi_ts), -1),
                        ylag2 = lag(diff(hpi_ts), -2),
                        xlag1 = lag(diff(new_spy), -1),
                        xlag2 = lag(diff(new_spy), -2),
                        xlag3 = lag(diff(new_spy), -3),
                        xlag4 = lag(diff(new_spy), -4),
                        xlag5 = lag(diff(new_spy), -5),
                        xlag6 = lag(diff(new_spy), -6),
                        dframe = TRUE)
df
model <- lm(y ~ xlag2 + xlag3, data=df, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 50)
plot(resid(model))
auto.arima(resid(model))
sarima(df$y, p=0,d=0,q=4, P=1,D=1,Q=0,S=12, xreg = cbind(df$xlag2, df$xlag3))


ccf2(diff(log(new_spy)), diff(log(hpi_ts)))
plot(log(new_spy))
acf2(log(new_spy), max.lag = 50)
model <- sarima(new_spy, p=1,d=0,q=0, P=0,D=0,Q=0,S=12)
log_res = model$fit$residuals
hpi_filtered <- filter(hpi_ts, filter = c(1, -model$ttable[, 1]), sides=1)
ccf2(log_res, hpi_filtered)
ccf2()

df <- ts.intersect(y = diff(hpi_ts),
                   ylag = lag(diff(hpi_ts), -1),
                   ylag2 = lag(diff(hpi_ts), -2),
                   xlag1 = lag(diff(new_spy), -1),
                   xlag2 = lag(diff(new_spy), -2),
                   xlag3 = lag(diff(new_spy), -3),
                   xlag4 = lag(diff(new_spy), -4),
                   xlag5 = lag(diff(new_spy), -5),
                   xlag6 = lag(diff(new_spy), -6),
                   dframe = TRUE)
model <- lm(y ~ xlag2 + xlag3, data=df, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 50)
plot(resid(model))
auto.arima(resid(model))
sarima(df$y, p=0,d=0,q=4, P=1,D=1,Q=0,S=12, xreg = cbind(df$xlag2, df$xlag3))

ccf2(diff(log(new_spy)), diff(log(hpi_ts)))
df <- ts.intersect(y = diff(log(hpi_ts)),
                   xlag2 = lag(diff(log(new_spy)), -2),
                   #xlag3 = lag(diff(log(new_spy)), -3),
                   dframe = TRUE)
model <- lm(y ~ xlag2, data=df, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 50)
plot(resid(model))
auto.arima(resid(model))
sarima(df$y, p=1,d=0,q=4, P=1,D=0,Q=0,S=12, xreg = cbind(df$xlag2))


df <- ts.intersect(y = diff(hpi_ts),
                   xlag2 = lag(diff(new_spy), -2),
                   xlag3 = lag(diff(new_spy), -3),
                   dframe = TRUE)
model <- lm(y ~ xlag2, data=df, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 50)
plot(resid(model))
auto.arima(resid(model))
sarima(df$y, p=0,d=0,q=5, P=1,D=1,Q=0,S=12, xreg = cbind(df$xlag2))


model <- lm(y ~ xlag2 + xlag3, data=df, na.action = na.omit)
summary(model)
acf2(resid(model), max.lag = 50)
plot(resid(model))
auto.arima(resid(model))
sarima(df$y, p=1,d=0,q=3, P=1,D=1,Q=0,S=12, xreg = cbind(df$xlag2, df$xlag3))


sarima(new_spy, p=1,d=0,q=0, P=0,D=0,Q=0,S=12)
sarima(new_spy, p=0,d=1,q=0, P=0,D=0,Q=0,S=12)




######### COMPARE MODELS ############

n <- length(hpi_ts)
train_hpi <- ts(hpi_ts[1:(n - 5)], start = c(2019, 8), frequency = 12)
test_hpi <- ts(hpi_ts[(n - 4):n], start = c(2024, 3), frequency = 12)
train_spy <- ts(new_spy[1:(n - 5)], start = c(2019, 8), frequency = 12)
test_spy <- ts(new_spy[(n - 4):n], start = c(2024, 3), frequency = 12)

#SARIMA
spy_forecast_sar <- sarima.for(train_spy, n.ahead = 5, p=0,d=1,q=0)
predicted_spy_sar <- ts(c(train_spy, spy_forecast_sar$pred), start = c(2019, 8), frequency = 12)

#ETS
spy_model_ets <- HoltWinters(train_spy, gamma = FALSE)
spy_forecast_ets <- forecast(spy_model_ets, 5)
predicted_spy_ets <- ts(c(train_spy, spy_forecast_ets$mean), start = c(2019, 8), frequency = 12)

spy_forecasts <- ts.intersect(true = test_spy,
                              sar = spy_forecast_sar$pred,
                              ets = spy_forecast_ets$mean)

mae_sar <- mae(spy_forecasts[, "true"], spy_forecasts[, "sar"])
rmse_sar <- rmse(spy_forecasts[, "true"], spy_forecasts[, "sar"])

mae_ets <- mae(spy_forecasts[, "true"], spy_forecasts[, "ets"])
rmse_ets <- rmse(spy_forecasts[, "true"], spy_forecasts[, "ets"])

spy_metrics_table <- data.frame(
  Method = c("sar", "ets"),
  MAE = c(mae_sar, mae_ets),
  RMSE = c(rmse_sar, rmse_ets)
)
spy_metrics_table

# Get SPY lags
diff_spy_etc <- diff(predicted_spy_sar)
spylag2 <- lag(diff_spy_etc, -2)
spylag3 <- lag(diff_spy_etc, -3)

# Split data into train and test sets
train_spy2 <- ts(spylag2[1:(n - 7)], start = c(2019, 10), frequency = 12)
test_spy2 <- ts(spylag2[(n - 6):n], start = c(2024, 3), frequency = 12)
train_spy3 <- ts(spylag3[1:(n - 8)], start = c(2019, 11), frequency = 12)
test_spy3 <- ts(spylag3[(n - 7):n], start = c(2024, 3), frequency = 12)

train_data <- ts.intersect(train_hpi = diff(train_hpi),
                           train_spy2 = train_spy2,
                           train_spy3 = train_spy3)
newxreg2 <- head(ts.intersect(spylag2 = test_spy2,
                              spylag3 = test_spy3),
                              5)
newxreg3 <- head(test_spy2, 5)


# SARIMA
model_sar <- arima(train_hpi, order = c(2, 1, 0),
                   seasonal = c(c(1, 0, 0), period = 12))
forecast_sar <- predict(model_sar, n.ahead = 5) 

# ETS
model_ets <- HoltWinters(train_hpi)
forecast_ets <- forecast(model_ets, 5)

# Our Model 2
model_olr2 <- arima(train_data[, "train_hpi"], order = c(1, 0, 3), 
                   seasonal = c(c(1, 1, 0), period = 12),
                   xreg = cbind(train_data[, "train_spy2"], train_data[, "train_spy3"]))
forecast_olr2 <- predict(model_olr2, n.ahead = 5, newxreg = newxreg2)
last_hpi <- tail(train_hpi, 1)
nondiff_preds2 <- diffinv(forecast_olr2$pred, xi = last_hpi)
forecast_olr2 <- nondiff_preds[-1]

# Our Model 3
model_olr3 <- arima(train_data[, "train_hpi"], order = c(0, 0, 5), 
                    seasonal = c(c(1, 1, 0), period = 12),
                    xreg = train_data[, "train_spy2"])
forecast_olr3 <- predict(model_olr3, n.ahead = 5, newxreg = newxreg3)
nondiff_preds3 <- diffinv(forecast_olr3$pred, xi = last_hpi)
forecast_olr3 <- nondiff_preds3[-1]

# All
all_forecasts <- ts.intersect(true = test_hpi,
                              sar = forecast1$pred,
                              ets = forecast2$mean,
                              olr2 = forecast_olr2,
                              olr3 = forecast_olr3)
all_forecasts

# Metrics
library(Metrics)

mae_sar <- mae(all_forecasts[, "true"], all_forecasts[, "sar"])
rmse_sar <- rmse(all_forecasts[, "true"], all_forecasts[, "sar"])

mae_ets <- mae(all_forecasts[, "true"], all_forecasts[, "ets"])
rmse_ets <- rmse(all_forecasts[, "true"], all_forecasts[, "ets"])

mae_olr2 <- mae(all_forecasts[, "true"], all_forecasts[, "olr2"])
rmse_olr2 <- rmse(all_forecasts[, "true"], all_forecasts[, "olr2"])

mae_olr3 <- mae(all_forecasts[, "true"], all_forecasts[, "olr3"])
rmse_olr3 <- rmse(all_forecasts[, "true"], all_forecasts[, "olr3"])

metrics_table <- data.frame(
  Method = c("SARIMA", "ETS", "Model 2", "Model 3"),
  MAE = c(mae_sar, mae_ets, mae_olr2, mae_olr2),
  RMSE = c(rmse_sar, rmse_ets, rmse_olr3, rmse_olr3)
)

print(metrics_table)

