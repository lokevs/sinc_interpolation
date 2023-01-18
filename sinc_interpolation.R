
###                       von Schmalensee 2023 Methods Ecol Evol                         ###
##  (use following code as you want, but cite the above study if used in a publication)   ##
##                                                                                        ##
##          Script containing the sinc function, the sinc interpolation functions,        ##
##                           and a tutorial for how to use them.                          ##
##                                                                                        ##
##       Author: Loke von Schmalensee, Department of Zoology, Stockholm University        ##
##                      Email: loke.von.schmalensee@zoologi.su.se                         ##
##                                                                                        ##
###                                                                                      ###


#### The functions ####

## The sinc function
sinc <- function(x) {
  result <- sin(pi * x) / (pi * x)
  result[is.nan(result)] <- 1
  result
}

## The sinc interpolation function:  returns a vector of interpolated values (between measurements).
## measurements is a vector of measurements (e.g. temperature) sampled at even time intervals ...
## which are assumed to be of length 1.
## interpolation_times is a vector of the time points with even intervals at which interpolated values ...
## should be calculated. For example, if the interval between these points is 0.5, one value is ...
## interpolated between each measurement.
sinc_interpolation <- function(measurements, interpolation_times) {
  mat <- array(                                                 # Create a matrix ...
    rep(interpolation_times, each = length(measurements)),      # of the interpolation points sequence repeated on every row ...
    c(length(measurements), length(interpolation_times)))       # with as many rows as the number of data points.
  mat = measurements * sinc(mat - seq(1, length(measurements))) # Apply sinc function and multiply by y-value.
  colSums(mat)                                                  # Sum each column.
}

## This is a slightly different, and less intuitive version of the interpolation function, which ...
## uses a method called "zero padding" in the frequency domain. It is much faster, and therefore ...
## better if analyzing very large datasets. It is ever so slightly less exact, but this effect is ...
## for all intents and purposes negligible. Compare results at the end of the script.
fast_sinc_interpolation <- function(measurements, interpolation_times) {
  
  n_input <- length(measurements) # Number of measurements
  n_input_even <- 2^(ceiling(log2(n_input))) # Next number that is a power of 2
  resolution <- 1/interpolation_times %% 1 # Pick out the interpolation factors for the desired interpolation times
  
  # Add zeros to the data to make the sample size a power of 2
  if(n_input != n_input_even){
    measurements <- c(measurements, rep(0, n_input_even - n_input)) 
  }
  
  # If interpolation should take place, run this.
  if(is.finite(min(resolution)) == T){
    
    resolution <- max(resolution[is.finite(resolution)]) # Interpolation resolution (number of times it enhances the resolution of the data)
    n_interpolation <- n_input_even * resolution # Calculate required number of interpolated values
    n_zeros <- n_interpolation - n_input_even # Number of zeros that should be padded in the frequency domain
    measurements_ft <- fft(measurements)  # Calculate the discrete Fourier transform of the measurements
    
    # Zero pad the transformed measurements (add zeros to the middle of the FFT'ed sequence)
    zero_padded_measurements_ft <- c(measurements_ft[1:(n_input_even/2)], rep(0, n_zeros),measurements_ft[((n_input_even/2)+1):(n_input_even)])
    
    # Calculate the inverse Fourier transform of the padded measurements and return values at the desired times
    return(Re(fft(zero_padded_measurements_ft, inverse = T) / n_input_even)[(interpolation_times-1) * resolution + 1])
  }
  
  # Otherwise return the data at the desired time points
  else{
    return(measurements[interpolation_times])
  }
}

## And finally, set the working directory to that of the script ...
## (this assumes that the script and the data are in the same folder ...
## and the code only works with R-studio).
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############


##########################################################
###                                                    ###
###   Reconstructing a thermal regime: a tutorial      ###
###                                                    ###
##########################################################


## Below follows some sections showing how sinc interpolation can be used to reconstruct ...
## thermal regimes at a higher resolution than they were sampled.


## The sections are:

# 1. DECONSTRUCTING THE SINC INTERPOLATION FUNCTION USING A SIMPLE EXAMPLE 
# 2. A MORE COMPLEX EXAMPLE
# 3. A REALISTIC EXAMPLE AND VISUALIZING THE FREQUENCY DOMAIN
# 4. A COMPARISON BETWEEN THE ORIGINAL AND THE FAST INTERPOLATION METHODS
# 5. IMPUTING MISSING MEASUREMENTS



#### 1. DECONSTRUCTING THE SINC INTERPOLATION FUNCTION USING A SIMPLE EXAMPLE  ####

## First, let's specify at what time intervals we want to sample the temperature (let's start with 0.85)
## I encourage the reader to play around with different values.
sampling_interval <- 0.85 # This can be changed to a sampling interval of choice 

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

##########

#### 2. A MORE COMPLEX EXAMPLE ####

## To see that the interpolation works for complex wave forms as well, we can repeat ...
## the previous exercise on one. Let's specify an arbitrary complex waveform (three sine waves):
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

##########

#### 3. A REALISTIC EXAMPLE AND VISUALIZING THE FREQUENCY DOMAIN ####

## Read a real temperature time series dataset (a part of one of the thermal regimes from ...
## von Schmalensee et al. 2022 Ecol. Lett.)
temperature <- read.delim("./temperature.txt")

## Here, temperatures are measured every 15 minutes. We begin by assuming this as the "true" thermal regime.
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

## Let's see if these results are in accordance with the Nyquist-Shannon sampling theorem by looking at ...
## the temperature data in the frequency domain:

## Make an empty plot
plot(-100,-100, ylim = (c(-8, 13)), xlim = c(-11.2, -2), xlab = "Fluctuation", ylab = "Amplitude (log scale)", xaxt = "n", yaxt = "n", bty = "n", cex.lab = 1.5, cex.axis = 1)

## It can be nice to plot the frequencies on the log-scale. Here, I have just calculated what ...
## log-transformed values correspond to fluctuation frequencies at familiar time-scales.
xtick<- c(log2((0.25/24)/30),log2((0.25/24)/7),log2((0.25/24)), log2((0.25/24)/(1/24)))
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3],
     labels = c("Monthly", "Weekly", "Daily", "Hourly"), pos = 1, xpd = TRUE, cex=1.5)
box(bty = "l")

## I have skipped plotting the y-axis, because we are not interested in the absolute values here.
## However, also the amplitudes will be plotted on a logarithmic scale for clarity (highlighting ...
## proportional differences instead of absolute differences)

## Extract the frequency spectrum using the function 'spectrum'
frequencies <- spectrum(temperature$temp, plot = F)$freq
amplitudes <- spectrum(temperature$temp, plot = F)$spec

## Plot it!
lines(log(amplitudes) ~ log2(frequencies), lwd = 1.5)

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

#### 4. A COMPARISON BETWEEN THE ORIGINAL AND THE FAST INTERPOLATION METHODS ####

## Continue using the real world example from von Schmalensee et al. 2022 Ecol. Lett.
temperature <- read.delim("./temperature.txt")

## Plot it again
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

## Generate data frame with interpolated data using the original function
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 15 minutes.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 1 / (4 * 11))), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 1))

## Plot the 15-minute temperatures interpolated using the original function
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")




## Generate data frame with interpolated data using the faster function
fast_interpolated_data <- data.frame(interpolated_temperature = fast_sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 15 minutes.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 1 / (4 * 11))), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 1))

## Again, plot the 15-minute temperatures interpolated using the fast function
lines(interpolated_temperature ~ time, data = fast_interpolated_data, lty = 1, lwd = 2, col = "red")

## Not a big difference.

## Let's look at the correlation between the interpolated temperatures from the two different functions
plot(y = fast_interpolated_data$interpolated_temperature, x = interpolated_data$interpolated_temperature,
     ylim = c(-3, 45), xlim = c(-3, 45), 
     pch = 21, cex = 1.5, cex.lab = 1.5, bg = "grey",
     ylab = "Fast interpolated temperature (deg. C)", xlab = "Original interpolated temperature (deg. C)", yaxt = "n", xaxt = "n", bty = "l")
tick <- c(0, 20, 40)
axis(side = 1, at = tick, labels = FALSE)
text(x = tick,  par("usr")[3], labels = tick, pos = 1, xpd = TRUE, cex = 1.5)
axis(side = 2, at = tick, labels = FALSE)
text(y = tick,  par("usr")[1], labels = tick, pos = 2, xpd = TRUE, cex = 1.5)

# Plot the 1:1 slope, representing a perfect match between interpolated and observed temperatures
abline(0, 1, lwd = 2, col = "red")

## Calculate R2 based on this 1:1 slope. 
resid_ss <- sum((resid(lm(fast_interpolated_data$interpolated_temperature ~ 0 + offset(interpolated_data$interpolated_temperature))))^2) # Residual variation from the 1:1 model
total_ss <- sum(resid(lm(temperature$temp ~ 1))^2)
R2 <- round(1 - resid_ss/total_ss, 3)

## Put it on the figure
legend("topleft", paste0("R-squared = ", R2), cex = 1.5, bty = "n")
## Basically identical

##########

#### 5. IMPUTING MISSING MEASUREMENTS ####

## Finally, we can attempt to impute missing data, since temperature time series ...
## sometimes have holes in it.


## Let's start by seeing WHAT WE CAN'T DO.
## Recall that the sinc function is normalized on the x-axis to match the sampling intervals ...
## leading to there being no interference with actual data points in the interpolation. Therefore, ...
## it is clear that we cannot use sinc-interpolation when samples are unevenly spaced. I visualize ...
## it below:

## Here is the previous artificial, simple, thermal regime:
sampling_interval <- 0.85 
end <- 11.5
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

## Now, let's sample this thermal regime at the previous sampling intervals ...
## but remove the middle measurement.
measurement_time_points <- seq(from = 1, to = end, by = sampling_interval)[-7]
measured_temperature <- 1.5 + sin(pi * measurement_time_points + 1.2)
measured_temperature[-7]


## Let's plot the measurements on top of the underlying thermal regime.
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Recall that sinc interpolation is performed by summing separate sinc functions over the time- ...
## period of interest. Let's plot the sinc functions again.
abline(0, 0, col = "red")
for(i in 1:length(measurement_time_points)){
  curve(measured_temperature[i] * sinc((x - measurement_time_points[i]) / sampling_interval),
        -1, end + 2, add = T, lty = 3, lwd = 2, col = "grey")
}

## Since one data point is missing, there is one less sinc function plotted here as compared to before.
## Another way to phrase it is that there is a sinc function with an amplitude of 0 fitted to the ...
## the time point where the measurement is missing. This is equivalent to a missing sinc function, ...
## since the sinc functions are summed in the final step of the interpolation, and summing with zero ...
## does nothing. This is also how we hack our data frame so we do not have to re-specify the function:
measurement_time_points <- seq(from = 1, to = end, by = sampling_interval)
measured_temperature <- 1.5 + sin(pi * measurement_time_points + 1.2)
measured_temperature[7] <- 0

## Plot the data points again (adding the 0 at the 7th time point)
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Do the interpolation
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 0.1 * sampling_interval))

## Plot it
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")

## A-ha! It does not work.


## Now let's focus on what WE CAN POTENTIALLY DO!
## Imagine if we have a relatively well-sampled thermal regime with just a single or a few missing ...
## data points - these missing data points would ruin the efficiency of the sinc interpolation. Can we ... 
## salvage this? Well, potentially. If the sampling is done at least at twice the Nyquist frequency, ...
## we could interpolate the temperature where it is missing using sinc interpolation, and add that ...
## interpolated value to the data. Then, we can run the rest of the interpolation with a higher resolution ...
## and efficiency. Like this:

## Here is the previous artificial, simple, thermal regime:
sampling_interval <- 0.42 # Increase sampling interval
end <- 11.5
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

## Remove the middle measurement (row 13) and plot them
measurement_time_points <- seq(from = 1, to = end, by = sampling_interval)
measured_temperature <- 1.5 + sin(pi * measurement_time_points + 1.2)
measured_temperature[13] <- 0
points(measured_temperature[-13] ~ measurement_time_points[-13], 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")
abline(0,0, col = "red")

## Still, of course, the interpolation fails
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 0.1 * sampling_interval))
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")

## However, if we interpolate using only half the samples (red dots and line):
measurement_time_points_2 <- measurement_time_points[seq(from = 2, to = length(measurement_time_points), by = 2)]
measured_temperature_2 <- measured_temperature[seq(from = 2, to = length(measurement_time_points), by = 2)]
points(measured_temperature_2 ~ measurement_time_points_2, 
       cex = 1, pch = 21, 
       bg = "red")
interpolated_data_2 <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature_2,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature_2),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points_2[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points_2[length(measurement_time_points_2)],
                                           by = 0.2 * sampling_interval))
lines(interpolated_temperature ~ time, data = interpolated_data_2, lty = 1, lwd = 2, col = "red")

## Not surprising, since the Nyquist-criterion is met, this interpolation already does a very good job.
## But can we make it do an even better job by imputing the value that was missing, and then use ...
## all the other data points?
measured_temperature[13] <- interpolated_data_2$interpolated_temperature[which(interpolated_data_2$time == measurement_time_points[13])]

## Plot the data points again (adding the imputed data at the 13th time point instead of 0)
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
points(measured_temperature ~ measurement_time_points, 
       cex = 1.5, pch = 21, 
       bg = "lightgrey")

## Do the interpolation
interpolated_data <- data.frame(interpolated_temperature = sinc_interpolation(measurements = measured_temperature,  # Input vector of measured temperature (sampling interval assumed to be 1).
                                                                              interpolation_times = seq(from = 1,   # Interpolate temperature every 0.1 sampling intervals.
                                                                                                        to = length(measured_temperature),
                                                                                                        by = 0.1)), 
                                time = seq(from = measurement_time_points[1], # Make another column with the time points on the original unit scale.
                                           to = measurement_time_points[length(measurement_time_points)],
                                           by = 0.1 * sampling_interval))

## Plot the interpolation using only half of the data points (red lines)
lines(interpolated_temperature ~ time, data = interpolated_data_2, lty = 1, lwd = 2, col = "red")

## Plot the interpolation using all the data points and the imputed data (orange line)
lines(interpolated_temperature ~ time, data = interpolated_data, lty = 2, lwd = 2, col = "darkorange")

## Add legend
legend("topleft", c("Without imputation", "With imputation"), lty = c(1, 2), col = c("red", "darkorange"), lwd = 2, cex = 1.5, bty = "n")

## Clearly, it does better than any of the other interpolations without the imputed data.
## This will be even more important when interpolating complex waveforms. By clever implementation of ...
## sinc interpolation to impute missing data (e.g. locally in long time-series), interpolation can ...
## be improved!

##########
