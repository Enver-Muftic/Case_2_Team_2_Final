# Load required libraries
library(ggplot2)
library(tseries)
library(forecast)
library(vars)

# Get the list of CSV files in the current working directory
csv_files <- list.files(pattern = "\\.csv$")

# Create a list to store the data frames
data_list <- list()

# Loop through each CSV file
for (file in csv_files) {
  # Read the CSV file and store it in the list
  data_list[[length(data_list) + 1]] <- read.csv(file)
}

# Extract relevant data from different data frames
# Read the data representing the date, industrial production, consumer price index, and federal funds rate.
# Adjust vector lengths and handle blank entries.

# Extract data from the first data frame
first_data_frame <- data_list[[1]]
date <- as.Date(first_data_frame$sasdate, format = "%m/%d/%Y")
date <- date[-(1:2)]
date <- date[1:(length(date) - 1)]
ip <- first_data_frame$INDPRO
ip <- ip[1:(length(ip) - 2)]
ip <- ip[-1]

# Extract data from the seventh data frame
seventh_data_frame <- data_list[[7]]
cpi <- seventh_data_frame$CPIAUCSL
cpi <- cpi[1:(length(cpi) - 8)]
cpi <- cpi[-1]

# Extract data from the sixth data frame
sixth_data_frame <- data_list[[6]]
ffr <- sixth_data_frame$FEDFUNDS
ffr <- ffr[1:(length(ffr) - 7)]
ffr <- ffr[-1]

# Create data frames for each variable
df1 <- data.frame(date, ip)
df2 <- data.frame(date, cpi)
df3 <- data.frame(date, ffr)

# Plots of the variables using ggplot
ggplot(df1, aes(x = date, y = ip)) +
  geom_line() +
  labs(title = "Industrial Production Over Time",
       x = "Date",
       y = "Industrial Production")

# Plot of Consumer Price Index over time
ggplot(df2, aes(x = date, y = cpi)) +
  geom_line() +
  labs(title = "Consumer Price Index Over Time",
       x = "Date",
       y = "Consumer Price Index")

# Plot of Federal Funds Rate over time
ggplot(df3, aes(x = date, y = ffr)) +
  geom_line() +
  labs(title = "Federal Funds Over Time",
       x = "Date",
       y = "Federal Funds")

# Conduct Augmented Dickey-Fuller (ADF) tests and perform differencing for industrial production
# Industrial production is initially non-stationary, check ADF test result and autocorrelation plot.
adf_test_ip <- adf.test(ip)
print(adf_test_ip)
acf_ip <- acf(ip)

# Perform differencing on industrial production and check ADF test again
diff_ip <- diff(log(ip))
adf_test_diff_ip <- adf.test(diff_ip)
print(adf_test_diff_ip)
df4 <- data.frame(date[-1], diff_ip)
# Plot the differentiated industrial production
ggplot(df4, aes(x = date[-1], y = diff_ip)) +
  geom_line() +
  labs(title = "Differentiated Industrial Production",
       x = "Date",
       y = "Differentiated Industrial Production")
acf_diff_ip <- acf(diff_ip)

# Similar process for Consumer Price Index: ADF test, autocorrelation plot, differencing, and ADF test again
adf_test_cpi <- adf.test(cpi)
print(adf_test_cpi)
acf_cpi <- acf(cpi)
diff_cpi <- diff(diff(log(cpi)))
acf_diff_cpi <- acf(diff_cpi)
adf_test_diff_cpi <- adf.test(diff_cpi)
df5 <- data.frame(date[-(1:2)], diff_cpi)
# Plot the differentiated consumer price index
ggplot(df5, aes(x = date[-(1:2)], y = diff_cpi)) +
  geom_line() +
  labs(title = "Differentiated Consumer Price Index",
       x = "Date",
       y = "Differentiated Consumer Price Index")
print(adf_test_diff_cpi)

# Similar process for Federal Funds Rate: ADF test, autocorrelation plot, and differencing
adf_test_ffr <- adf.test(ffr)
print(adf_test_ffr)
acf_ffr <- acf(ffr)
diff_ffr <- diff(ffr)

# Perform Augmented Dickey-Fuller (ADF) test on the differenced Federal Funds Rate
adf_test_diff_ffr <- adf.test(diff_ffr)
print(adf_test_diff_ffr)

# Create a data frame for the differenced Federal Funds Rate
df6 <- data.frame(date[-1], diff_ffr)

# Plot the differenced Federal Funds Rate
ggplot(df6, aes(x = date[-1], y = diff_ffr)) +
  geom_line() +
  labs(title = "Differentiated Federal Funds Rate",
       x = "Date",
       y = "Differentiated Federal Funds Rate")

# Calculate autocorrelation function for differenced Federal Funds Rate
acf_diff_ffr <- acf(diff_ffr, lag.max = 20, plot = FALSE)

# Determine VAR order using VARselect
# Select the lag order based on the AIC criterion
df7 <- data.frame(diff_ip[-1], diff_cpi, diff_ffr[-1])
lag_selection <- VARselect(df7, lag.max = 24, type = "const")
print(lag_selection)

# Set lag order based on the AIC criterion
lag_order <- 3

# Fit VAR model
var_model <- VAR(df7, p = lag_order, type = "const")
print(var_model)

# Examine autocorrelation and residuals of the VAR model
plot(residuals(var_model))
acf(residuals(var_model))

# Examine normality of residuals using Shapiro-Wilk test and Q-Q plot
shapiro.test(residuals(var_model))
qqnorm(residuals(var_model))
qqline(residuals(var_model))

# Examine model stability using the stability function
stability(var_model)

# Perform Granger causality test
# Function to perform Granger causality test
granger_causality_test <- function(y, x, max_lag) {
  n <- length(y)
  y_lagged <- c(rep(NA, max_lag), y[1:(n - max_lag)])
  
  # Fit the model without the potential cause variable
  lm_without_x <- lm(y ~ y_lagged)
  rss_without_x <- sum(resid(lm_without_x)^2)
  
  # Fit the model with the potential cause variable
  lm_with_x <- lm(y ~ y_lagged + x)
  rss_with_x <- sum(resid(lm_with_x)^2)
  
  # Calculate F-statistic
  df_diff = length(lm_with_x$coefficients) - length(lm_without_x$coefficients)
  f_statistic <- ((rss_without_x - rss_with_x) / df_diff) / (rss_with_x / (n - length(lm_with_x$coefficients)))
  
  # Calculate p-value
  p_value <- 1 - pf(f_statistic, df_diff, n - length(lm_with_x$coefficients))
  
  # Display results
  cat("F-Statistic:", f_statistic, "\n")
  cat("P-Value:", p_value, "\n\n")
}

# Perform Granger causality tests for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    if (i != j) {
      cat("Granger causality test for", colnames(df7[i]), "on", colnames(df7[j]), ":\n")
      granger_causality_test(df7[[i]], df7[[j]], 3)  # Assuming a maximum lag of 3
    }
  }
}

# Calculate impulse response functions (irf) for every pair of variables and plot them. Plot 15 periods ahead.
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(var_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 15) 
    plot(irf, main = paste("Impulse variable: ", colnames(var_model$y)[i], " Response variable: ", colnames(var_model$y)[j]), xlab = "", ylab = "")
    
    # Plot leveled irfs. For each variable undo transformations done to ensure stationarity.
    if (i == 1) {
      leveled_irf <- exp(cumsum(irf$irf$diff_ip..1.))
      plot(leveled_irf, main = paste("Leveled IRF, Impulse variable: ", colnames(var_model$y)[i], " Response variable: ", colnames(var_model$y)[j]), xlab = "", ylab = "", type = "l")
    }
    if (i == 2) {
      leveled_irf <- exp(cumsum(cumsum(irf$irf$diff_cpi)))
      plot(leveled_irf, main = paste("Leveled IRF, Impulse variable: ", colnames(var_model$y)[i], " Response variable: ", colnames(var_model$y)[j]), xlab = "", ylab = "", type = "l")
    }
    if (i == 3) {
      leveled_irf <- cumsum(irf$irf$diff_ffr..1.)
      plot(leveled_irf, main = paste("Leveled IRF, Impulse variable: ", colnames(var_model$y)[i], " Response variable: ", colnames(var_model$y)[j]), xlab = "", ylab = "", type = "l")
    } 
  }
}

# Create an identity matrix with NAs on the diagonal
amat <- diag(3)
diag(amat) <- NA 

# Fit structural VAR (svar) model using scoring algorithm
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for every pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}

# Create a data frame combining different orderings of variables for VAR and SVAR models

# Order: diff_ip, diff_ffr, diff_cpi
df8 <- data.frame(diff_ip[-1], diff_ffr[-1], diff_cpi)
var_model <- VAR(df8, p = lag_order, type = "const")
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}

# Repeat the process for different variable orderings

# Order: diff_cpi, diff_ffr, diff_ip
df9 <- data.frame(diff_cpi, diff_ffr[-1], diff_ip[-1])
var_model <- VAR(df9, p = lag_order, type = "const")
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}

# Repeat the process for different variable orderings

# Order: diff_cpi, diff_ip, diff_ffr
df10 <- data.frame(diff_cpi, diff_ip[-1], diff_ffr[-1])
var_model <- VAR(df10, p = lag_order, type = "const")
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}

# Repeat the process for different variable orderings

# Order: diff_ffr, diff_ip, diff_cpi
df11 <- data.frame(diff_ffr[-1], diff_ip[-1], diff_cpi)
var_model <- VAR(df11, p = lag_order, type = "const")
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}

# Repeat the process for different variable orderings

# Order: diff_ffr, diff_cpi, diff_ip
df12 <- data.frame(diff_ffr[-1], diff_cpi, diff_ip[-1])
var_model <- VAR(df12, p = lag_order, type = "const")
svar_model <- SVAR(var_model, estmethod = "scoring", Amat = amat, Bmat = NULL, max.iter = 435)

# Plot impulse response functions for each pair of variables
for (i in 1:3) {
  for (j in 1:3) {
    irf <- irf(svar_model, impulse = colnames(var_model$y)[i], response = colnames(var_model$y)[j], n.ahead = 10) 
    plot(irf, main = paste("Impulse Response Function, impulse variable: ", colnames(var_model$y)[i], " and response variable: ", colnames(var_model$y)[j]))
  }
}