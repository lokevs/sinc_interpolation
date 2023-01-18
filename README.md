# sinc_interpolation

An R-script (sinc_interpolation.R) with the sinc function and the sinc interpolation function. The sinc interpolation function can be used to generate high-resolved temperature data from sparse measurements.

The script contains a step-by-step guide for employing sinc interpolation over environmental temperature (or other) data, with figures along the way, including how to Fourier transform temperature data and visualize it in the frequency domain, and how to impute missing temperatures. The script has no dependencies but utilizes some data from temperature.txt (a temperature time-series).

Contents:

0. THE SINC AND SINC INTERPOLATION FUNCTIONS
1. DECONSTRUCTING THE SINC INTERPOLATION FUNCTION USING A SIMPLE EXAMPLE 
2. A MORE COMPLEX EXAMPLE
3. A REALISTIC EXAMPLE AND VISUALIZING THE FREQUENCY DOMAIN
4. A COMPARISON BETWEEN THE ORIGINAL AND THE FAST INTERPOLATION METHODS
5. IMPUTING MISSING MEASUREMENTS
