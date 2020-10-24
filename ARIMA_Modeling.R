# Import packages
library(readr)
library(ggplot2)
library(dplyr)
library(tseries)
library(forecast)
library(haven)
library(fma)
library(expsmooth)
library(lmtest)
library(zoo)
library(seasonal)
library(lubridate)
library(tidyr) 
library(ggthemes)

# import data and convert Date column from character to date 
data <- read_csv('~/MSA 21/AA 502/Time Series/Ozone_Raleigh2.csv', 
                 col_types = cols(Date = col_date(format = '%m/%d/%Y')))

# Change Daily Max column name so it is shorter and easier to handle
names(data)[5] <- "maxOzoneC"

# Create a time plot of the mean monthly max 8 hour ozone concentration and upload it.
monthly_data <- data %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(mean_ozone = mean(maxOzoneC, na.rm = T))

# Create monthly ts object
monthly_ts <- ts(monthly_data$mean_ozone, start = c(2014, 1), end = c(2020, 5), frequency = 12)

# Training data:
trainingTS <- subset(monthly_ts, end = length(monthly_ts) - 17)

# Validation data:
validTS <- subset(monthly_ts, start = length(monthly_ts) - 16, end = length(monthly_ts) - 5)

# Test data
testTS <- subset(monthly_ts, start = length(monthly_ts) - 4)

# Print out 10 decimal places
options(digits=10)

# Change plot size output in jupyter notebook
options(repr.plot.width=4, repr.plot.height=3)

# Visualize data to see what you have-- definitely seasonality!
plot(trainingTS)


#### ORIGINAL DATA PLOT ######
trainingTS %>% autoplot() + 
  labs(title = 'Average Monthly Ozone Concentration Over Time', 
       y = 'Average Monthly Ozone Concentration')+
  theme_tufte()+
  theme(plot.title = element_text(hjust =.5))

#########################################################################################################

### Dummy Variable Method: better for shorter season lengths ###
# Since we have monthly data, we can fit 11 dummy variables (12 - 1)
# Set the month variable in training set to a factor (60 obs in training data); factor is needed for model.matrix()
month.factor <- factor(monthly_data$month[1:60])

# Create a model matrix on this month variable: returns 11 dummy variables
# Make design matrix with one-sided formula bc we don't need to know the response to be able to determine the model matrix
# two-sided formula will still work but left-hand side is ignored
reg.ozone <- model.matrix(~month.factor)
reg.ozone <- reg.ozone[,-1] # delete the first column (y-int)

# Fit a seasonal ARIMA model on time series (training), using dummy variables as external regressors
reg.dummies <- Arima(trainingTS, xreg = reg.ozone)
summary(reg.dummies) # MAE = 0.00214, MAPE = 5.736

### Sin/Cosine Method: better for longer season lengths ###
reg.trig <- Arima(trainingTS, xreg = fourier(trainingTS, K=4))
summary(reg.trig) # MAE = 0.00218, MAPE = 5.791

##########################################################################################################

# Check residuals on dummies model
options(repr.plot.width=4, repr.plot.height=3) # adjust plot size
# checkresiduals(arima.dummies); fancier function to display residuals
# plot(reg.dummies$residuals, main = 'Dummy Method')

####### RESIDUAL PLOTS #########

autoplot(reg.dummies$residuals, size = .8) + 
  theme_tufte() + 
  labs(title = 'Residuals of Dummy Variable Model', 
       y = 'Residuals') +
  theme(plot.title = element_text(hjust =.5))

ggsave('~/MSA 21/AA 502/Time Series/ResidualsDummies.png', width = 11, height = 3)


autoplot(reg.trig$residuals, size = .8) + 
  theme_tufte() + 
  labs(title = 'Residuals of Fourier Model', 
       y = 'Residuals') +
  theme(plot.title = element_text(hjust =.5))
  
ggsave('~/MSA 21/AA 502/Time Series/ResidualsTrig.png', width = 11, height = 3)


# Check residuals on trig model
# checkresiduals(arima.trig); fancier function to display residuals
# plot(reg.trig$residuals, main = 'Trig Method')

# See if residuals are stationary
ndiffs(reg.dummies$residuals)
ndiffs(reg.trig$residuals)

# Both dummies and trig methods have stationary residuals since 0 differences need to be taken.

# # Code for plotting ACF and PACF plots
# options(repr.plot.width=4, repr.plot.height=3)
# Acf(arima.trig$residuals, lag = 10)$acf
# Pacf(arima.trig$residuals, lag = 10)$acf

##########################################################################################################

# Get best AR and MA terms for both models
dummies_fit <- auto.arima(reg.dummies$residuals)
dummies_fit # 1 nonseasonal AR term

trig_fit <- auto.arima(reg.trig$residuals)
trig_fit # 0 nonseasonal terms

##########################################################################################################

# Models with AR and MA terms
arima.dummies <- Arima(trainingTS, order=c(1,0,0), xreg = reg.ozone)
summary(arima.dummies)

arima.trig <- Arima(trainingTS, order=c(0,0,0), xreg = fourier(trainingTS, K=4))
summary(arima.trig)

##########################################################################################################

#------------------------------------- White Noise for Dummies Model -----------------------------------#
# Perform Ljung-Box test up to 20 lags on residuals then store p-values in white.LB_resid
# Number of lags observed has increased from HW 3 to capture more seasonality 
# fitdf=0 because we did not fit any nonseasonal AR or MA terms

White.LB_resid <- rep(NA, 20)

for(i in 1:20){
  White.LB_resid[i] <- Box.test(arima.dummies$residuals, lag = i, type = "Lj", fitdf = 0)$p.value
}


# p-values >= 0.2 are recorded as 0.2 (for plotting purposes)
White.LB_resid <- pmin(White.LB_resid, 0.2)
LBDummy = White.LB_resid

# Let's look at a plot of these p-values (lags 1,2,...,10)
# The horizontal lines let us see which lags have p-values <0.05 and <0.01
# barplot(White.LB_resid, 
#         main = "Ljung-Box Test for Dummies Residuals", 
#         ylab = "Probabilities", 
#         xlab = "Lags", 
#         ylim = c(0, 0.2),
#         names.arg = seq(0,19))
# 
# abline(h = 0.01, lty = "dashed", col = "black")
# abline(h = 0.05, lty = "dashed", col = "black")


#individual LB plot for Dummy Model
data.frame(cbind(Probabilities = White.LB_resid, Lags = seq(1,20))) %>% 
  ggplot( aes(x = Lags, y = Probabilities)) + 
  geom_col() +
  geom_hline(yintercept = .05, linetype = 'dashed')+
  geom_hline(yintercept = .01, linetype = 'dashed')+
  labs(title = 'Ljung-Box Test for Dummy Variable Model Residuals')+
  theme_tufte()+
  theme(plot.title = element_text(hjust = .5))
  
#----------------------------------------- White Noise for Trig Model ---------------------------------- #
White.LB_resid <- rep(NA, 20)

for(i in 1:20){
  White.LB_resid[i] <- Box.test(arima.trig$residuals, lag = i, type = "Lj", fitdf = 0)$p.value
}

# p-values >= 0.2 are recorded as 0.2 (for plotting purposes)
White.LB_resid <- pmin(White.LB_resid, 0.2)
LBTrig = White.LB_resid

# Let's look at a plot of these p-values (lags 1,2,...,10)
# The horizontal lines let us see which lags have p-values <0.05 and <0.01
# barplot(White.LB_resid, 
#         main = "Ljung-Box Test for Trig Residuals", 
#         ylab = "Probabilities", 
#         xlab = "Lags", 
#         ylim = c(0, 0.2),
#         names.arg = seq(0,19))
# 
# abline(h = 0.01, lty = "dashed", col = "black")
# abline(h = 0.05, lty = "dashed", col = "black")


#Individual LB plot for Fourier
data.frame(cbind(Probabilities = White.LB_resid, Lags = seq(1,20))) %>% 
  ggplot( aes(x = Lags, y = Probabilities)) + 
  geom_col() +
  geom_hline(yintercept = .05, linetype = 'dashed')+
  geom_hline(yintercept = .01, linetype = 'dashed')+
  labs(title = 'Ljung-Box Test for Fourier Model Residuals')+
  theme_tufte()+
  theme(plot.title = element_text(hjust = .5))


# Side by side LB plots
data.frame(cbind(`Dummy Variable Model` = LBDummy, Lags = seq(1,20), `Fourier Model` = LBTrig)) %>% 
  gather(key = 'Model', value = 'Probabilities', -Lags) %>% 
  ggplot( aes(x = Lags, y = Probabilities)) + 
  geom_col() +
  geom_hline(yintercept = .05, linetype = 'dashed')+
  geom_hline(yintercept = .01, linetype = 'dashed')+
  facet_wrap(~Model, nrow = 1)+
  labs(title = 'Ljung-Box Test for Model Residuals')+
  theme_tufte()+
  theme(plot.title = element_text(hjust = .5))

#########################################################################################################
#------------------------------ Forecast on validation data using Dummies--------------------------------#
dummies_pred <- forecast(arima.dummies, xreg = reg.ozone)
dummies_pred_v <- dummies_pred$mean[1:12] # predictions for validation

# Calculate validation data error
dummies_pred_v_error <- validTS - dummies_pred_v
# Validation MAE = 0.00242
MAE_dummies_v <- mean(abs(dummies_pred_v_error))
# Validation MAPE = 5.789
MAPE_dummies_v <- mean(abs(dummies_pred_v_error)/abs(validTS))*100


#------------------------------ Forecast on validation data using Trig ----------------------------------#
trig_pred <- forecast(arima.trig, xreg = fourier(trainingTS, K=4))
trig_pred_v <- trig_pred$mean[1:12] # predictions for validation

# Calculate validation data error
trig_pred_v_error <- validTS - trig_pred_v
# Dummies validation MAE = 0.00221
MAE_trig_v <- mean(abs(trig_pred_v_error))
# Trig validation MAPE = 5.401
MAPE_trig_v <- mean(abs(trig_pred_v_error)/abs(validTS))*100

######################################### Choose final model ############################################

#---------------------- Forecast predictions using Trig model on testing (5 months)---------------------#

trig_pred_test <- trig_pred$mean[13:17] # predictions for testing data
# Calculate test data error
trig_pred_test_error <- testTS - trig_pred_test
# Trig test MAE = 0.00374
MAE_trig_test <- mean(abs(trig_pred_test_error))
# Trig test MAPE = 9.789
MAPE_trig_test <- mean(abs(trig_pred_test_error)/abs(testTS))*100

#------------------------------- Plot actual vs predicted for validation data --------------------------#

# Create df with predicted values, actual values (test data), and date for validation set
results_df <- data.frame(trig_pred_v = trig_pred_v,
                         dummies_pred_v = dummies_pred_v,
                         actual = validTS,
                         date = seq(as.Date('2019-01-01'), as.Date('2019-12-01'), by="months"))

# Validation Predicted vs. Actual
results_df %>% 
  select(date, actual, trig_pred_v, dummies_pred_v) %>% 
  gather(key = 'Model', value = value, -date) %>% 
  mutate(Model = case_when(Model == 'actual' ~ 'Actual Value',
                           Model == 'dummies_pred_v' ~ 'Dummy Variable Model Prediction',
                           Model == 'trig_pred_v' ~ 'Fourier Model Prediction')) %>% 
  ggplot(aes(x = date, y = value, color = Model)) +
  geom_line(size = 1, aes(linetype = Model)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype_manual(values = c('solid', 'dotted', 'dashed'))+
  theme_tufte()+ 
  labs(x = 'Date', 
       y = 'Monthly Average Ozone Concentration',
       title = 'Ozone Concentration Model Comparison',
       subtitle = 'Validation Data') +
  theme(legend.position = c(.77,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'white', size = .5, linetype = 'solid'),
        plot.title = element_text( hjust =.5), 
        plot.subtitle = element_text( hjust =.5))

#------------------------------- Plot actual vs predicted for the test data ----------------------------#

# Create df with predicted values, actual values (test data), and date for the 5 month forecast (test)
results_df <- data.frame(trig_pred = trig_pred_test,
                         dummies_pred = dummies_pred$mean[13:17],
                         actual = testTS,
                         date = seq(as.Date('2020-01-01'), as.Date('2020-05-01'), by="months"))
# Color test data
results_df %>% 
  select(date, actual, trig_pred) %>% 
  gather(key = 'Model', value = value, -date) %>% 
  mutate(Model = case_when(Model == 'actual' ~ 'Actual Value',
                           Model == 'dummies_pred' ~ 'Dummy Variable Model Prediction',
                           Model == 'trig_pred' ~ 'Fourier Model Prediction')) %>% 
  ggplot(aes(x = date, y = value, color = Model)) +
  geom_line(size = 1, aes(linetype = Model)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype_manual(values = c('solid', 'dashed'))+
  theme_tufte()+ 
  labs(x = 'Date', 
       y = 'Monthly Average Ozone Concentration',
       title = 'Ozone Concentration Model Comparison',
       subtitle = 'Test Data') +
  theme(legend.position = c(.10,.95),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'white', size = .5, linetype = 'solid'),
        plot.title = element_text( hjust =.5), 
        plot.subtitle = element_text( hjust =.5)) 


# Black and White only
results_df %>% 
  select(date, actual, trig_pred) %>% 
  gather(key = 'Model', value = value, -date) %>% 
  mutate(Model = case_when(Model == 'actual' ~ 'Actual Value',
                           Model == 'dummies_pred' ~ 'Dummy Variable Model Prediction',
                           Model == 'trig_pred' ~ 'Fourier Model Prediction')) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_line(size = 1, aes(linetype = Model)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype_manual(values = c('solid', 'dashed'))+
  theme_tufte()+ 
  labs(x = 'Date', 
       y = 'Monthly Average Ozone Concentration',
       title = 'Ozone Concentration Model Comparison',
       subtitle = 'Test Data') +
  theme(legend.position = c(.10,.95),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'white', size = .5, linetype = 'solid'),
        plot.title = element_text( hjust =.5), 
        plot.subtitle = element_text( hjust =.5))