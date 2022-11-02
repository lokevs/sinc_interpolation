
###                              von Schmalensee 2022                         ###
##                                                                             ##
##     Script containing the sinc function, the sinc interpolation function,   ##
##                    and a tutorial for how to use them.                      ##
##                                                                             ##
##  Author: Loke von Schmalensee, Department of Zoology, Stockholm University  ##
##                  Email: loke.von.schmalensee@zoologi.su.se                  ##
##                                                                             ##
###                                                                           ###



#### The functions ####

## The sinc function
sinc <- function(x) {
  result <- sin(pi * x) / (pi * x)
  result[is.nan(result)] <- 1
  result
}

## The sinc interpolation function:  returns a vector of interpolated values (between measurements).
# measurements is a vector of measurements (e.g. temperature) sampled at even time intervals ...
# which are assumed to be of length 1.
# interpolation_times is a vector of the time points with even intervals at which interpolated values ...
# should be calculated. For example, if the interval between these points is 0.5, one value is ...
# interpolated between each measurement.
sinc_interpolation <- function(measurements, interpolation_times) {
  mat <- array(                                                 # Create a matrix ...
    rep(interpolation_times, each = length(measurements)),      # of the interpolation points sequence repeated on every row ...
    c(length(measurements), length(interpolation_times)))       # with as many rows as the number of data points.
  mat = measurements * sinc(mat - seq(1, length(measurements))) # Apply sinc function and multiply by y-value.
  colSums(mat)                                                  # Sum each column.
}

############



#### Reconstructing a thermal regime: a tutorial ####

## First, let's specify at what time intervals we will sample temperature (let's start with 0.7)
sampling_interval <- 0.85 # This can be changed to a sampling interval of choice



### A simple example

## Assume that temperature is fluctuating in a very simple way, ...
## like this sine wave with a wavelength of 2, for example: 
## 1.5 + sin(pi * x + 1.2)

## Let's start with plotting our artificial, simple, thermal regime from time point 1 to time point "end":
end <- 11.5 # This can be changed to represent a measurement period of choice
curve(1.5 + sin(pi * x + 1.2), 
      1, end, 
      ylim = c(-0.5, 3), xlim = c(0.5, end + 0.5), 
      lwd = 3, cex.lab = 1.5, 
      ylab = "Temperature", xlab = "Time",  xaxt = "n", yaxt = "n", bty = "l")
xtick <- seq(from = 1, to = end, by = 2)
ytick<- c(0)
axis(side = 1, at = xtick, labels = FALSE)
text(x = xtick,  par("usr")[3], labels = xtick, pos = 1, xpd = TRUE, cex = 1.5)
axis(side = 2, at = ytick, labels = FALSE)
text(y = ytick,  par("usr")[1], labels = ytick, pos = 2, xpd = TRUE, cex = 1.5)

## Now, let's sample this thermal regime at the specified sampling intervals.
measurement_time_points <- seq(from = 1, to = end, by = sampling_interval)
measured_temperature <- 1.5 + sin(pi * measurement_time_points + 1.2)

## Let's plot the measurements on top of the underlying thermal regime.
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Now, we want to perform sinc interpolation to reconstruct the underlying thermal regime.
## This is done by summing separate sinc functions over the time period of interest.
## Each sinc function is scaled on the y-axis and shifted on the x-axis to one unique measurement.
## Let's start by looking at one separate sinc function to get an intuition for what it does.
curve(measured_temperature[5] * # Scaling the temperature value to match the measurement
        sinc((x - measurement_time_points[5]) / sampling_interval), # Divide by the sampling interval to normalize on the x-axis
      -1, end + 2, add = T, lty = 3, lwd = 2, col = "grey")

## Notice how at each measurement time, the sinc function evaluates to zero.
## Let's draw a horizontal line at zero to demonstrate this
abline(0, 0, col = "red")

## Repeat the sinc function fitting for each measurement.
for(i in 1:length(measurement_time_points)){
  curve(measured_temperature[i] * sinc((x - measurement_time_points[i]) / sampling_interval),
        -1, end + 2, add = T, lty = 3, lwd = 2, col = "grey")
}

## Because each function evaluates to 0 at all other measurement time points ...
## the data speaks for itself at each measurement. But between the measurements we must interpolate.
## This is done by summing all the sinc functions, which we can do with the sinc_interpolation function.
## Let's use the sinc_interpolation function to generate a data frame of interpolated temperature.
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 0.1 * sampling_interval))

## Let's now plot the interpolated temperatures
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")

## Notice how the interpolation is most accurate in the middle of the sampling. 
## This is because perfect reconstruction can theoretically only be done over infinite time series.
## However, for most ecological applications, interpolated data converge sufficiently with the true data ...
## after a few sampling intervals (given that temperature is sampled frequently enough)



### A more complex example

## To see that the interpolation works for complex wave forms as well, we can repeat ...
## the exercise on one. Let's specify an arbitrary complex waveform (three sine waves):
## 2.5 + sin(pi * x + 1.2) + 0.3 * sin(pi * x * 1.3) + 1.2 * sin(pi * x * 0.4)

## Plot it
curve(2.5 + sin(pi * x + 1.2) + 0.3 * sin(pi * x * 1.3) + 1.2 * sin(pi * x * 0.4),
      1, end, 
      ylim = c(-0.5, 5), xlim = c(0.5, end + 0.5), 
      lwd = 3, cex.lab = 1.5, 
      ylab = "Temperature", xlab = "Time",  xaxt = "n", yaxt = "n", bty = "l")

## Sample it again
measurement_time_points <- seq(from = 1, to = end, by = sampling_interval)
measured_temperature <- 2.5 + sin(pi * measurement_time_points + 1.2) + 0.3 * sin(pi * measurement_time_points * 1.3) + 1.2 * sin(pi * measurement_time_points * 0.4)

## Plot the measurements
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Generate data frame with interpolated data
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 0.1 * sampling_interval))

## Let's now plot the interpolated temperatures
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")
## Even if the second sine wave fluctuates at frequencies more then twice above the sampling rate, the interpolation ...
## does a pretty good job of retrieving the variation.



### A realistic example

## Read a real temperature time series dataset (a part of one of the thermal regimes from ...
## von Schmalensee et al. 2022 Ecol. Lett.)
temperature <- read.delim("./temperature.txt")

## Here, temperatures are measured every 15 minutes. We begin by assuming as the "true" thermal regime.
plot(temp ~ time,
     data = temperature,
     ylim = c(-3, 45), xlim = c(0.5, nrow(temperature) + 0.5), 
     type = "l", lwd = 3, cex.lab = 1.5, 
     ylab = "Temperature (deg. C)", xlab = "Time (day)", yaxt = "n", xaxt = "n", bty = "l")
xtick <- seq(from = 1, to = nrow(temperature), by = 4 * 24 * 2)
ytick<- c(0, 20, 40)
axis(side = 1, at = xtick, labels = FALSE)
text(x = xtick,  par("usr")[3], labels = (1:length(xtick) - 1) * 2 + 1, pos = 1, xpd = TRUE, cex = 1.5)
axis(side = 2, at = ytick, labels = FALSE)
text(y = ytick,  par("usr")[1], labels = ytick, pos = 2, xpd = TRUE, cex = 1.5)

## Artificially thin the dataset to reflect a sparser sampling scheme
measurement_time_points <- seq(from = 1, to = nrow(temperature), by = 4 * 11) # Let's sample every 11 hours
measured_temperature <- temperature$temp[measurement_time_points]

## Plot the measurements
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Generate data frame with interpolated data
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 15 minutes.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 1 / (4 * 11))), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 1))

## Plot the 15-minute temperatures interpolated using 11-hour measurement intervals
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")

## Finally, let's look at the correlation between measured and interpolated 15-minute temperatures
plot(temperature$temp[1:nrow(interpolated_data)] ~ interpolated_data$interpolated_temperature,
     data = temperature,
     ylim = c(-3, 45), xlim = c(-3, 45), 
     pch = 21, cex = 1.5, cex.lab = 1.5, bg = "grey",
     ylab = "Observed temperature (deg. C)", xlab = "Interpolated temperature (deg. C)", yaxt = "n", xaxt = "n", bty = "l")
tick <- c(0, 20, 40)
axis(side = 1, at = tick, labels = FALSE)
text(x = tick,  par("usr")[3], labels = tick, pos = 1, xpd = TRUE, cex = 1.5)
axis(side = 2, at = tick, labels = FALSE)
text(y = tick,  par("usr")[1], labels = tick, pos = 2, xpd = TRUE, cex = 1.5)

# Plot the 1:1 slope, representing a perfect match between interpolated and observed temperatures
abline(0, 1, lwd = 2, col = "red")

## Calculate R2 based on this 1:1 slope. This represents how much of the observed time-specific 
## variation in temperature that you will actually capture if you believe that 15-minute temperatures
## can be reconstructed using 11-hour sampling intervals.
resid_ss <- sum((resid(lm(temperature$temp[1:nrow(interpolated_data)] ~ 0 + offset(interpolated_data$interpolated_temperature))))^2) # Residual variation from the 1:1 model
total_ss <- sum(resid(lm(temperature$temp ~ 1))^2)
R2 <- round(1 - resid_ss/total_ss, 3)

## Put it on the figure
legend("topleft", paste0("R-squared = ", R2), cex = 1.5, bty = "n")
## As can be seen, 11-hour sampling intervals can retrieve 74% (!) of the observed, time-specific, 15-minute temperature variation!

##########