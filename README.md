# Ozone Air Quality: Project Overview

We want to explore the average monthly concentration of ozone (ppm) in a local area using a seasonal ARIMA model. We will address how seasonality is incorporated into the model for ozone and how this decision was made.

### Part 1: Data Preparation 
* The data recorded daily max 8-hour ozone concentration per observation.
* Rolled up data to the monthly level by using the average of the observations present
* Evaluated time plot for potential trends or seasonality captured

### Part 2: Model Building
* Split the data into training (withholding last 17 months), validation (next 12 months), test (last 5 months) data sets
* Explored dummy variables and trigonometric functions for capturing seasonality
* Determined Autoregressive Moving Average terms and incorporated into a monthly ARIMA model with the training data set

### Part 3: Testing Stationarity and White Noise
* Checked for stationarity of residuals including any potential trend and/or random walks
* Tested whether the stationary series exhibited white noise using the Ljung-Box test

### Part 4: Forecast Ozone Data
* Made data visualizations of the following:
  * Time Plot of the predicted versus actual for the validation and test data
* Calculated MAPE and MAE values for evaluating accuracy of forecasts with test data
