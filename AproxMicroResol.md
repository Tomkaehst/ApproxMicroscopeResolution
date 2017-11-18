Approximating Microscopic Resolution
================

Used Data
---------

This document aims at approximating the resolution of an microscopic image by measureing the intensity profile of a structure. The data used are images of EGFP-RAD51-overexpressing U2OS cells after induction of a DNA damage by bleomycin. RAD51 is the recombinase, which binds a damaged DNA site and is involved in strand probing and exchange in the process of homologous recombination during the S and G2 phases of the cell cycle.

Theoretically, the resolution according to the Rayleigh criterion
$$d = \\frac{0.61 \\lambda}{NA}$$
 would be ~ **330** nm using a 0.95 NA water immersion objective and an EGFP emission of **514 nm**.

![My Figure](rad51.jpg)

The data used is the plotted intensity profile perpendicular to a RAD51 nucleoprotein filament.

``` r
dat = read.csv("Values.csv")
print(dat)
```

    ##           X       Y
    ## 1  0.000000 24.0000
    ## 2  0.094512 24.4492
    ## 3  0.189024 29.4746
    ## 4  0.283536 34.7160
    ## 5  0.378048 35.3546
    ## 6  0.472560 33.5960
    ## 7  0.567071 29.3889
    ## 8  0.661583 26.7037
    ## 9  0.756095 24.6468
    ## 10 0.850607 24.5000

``` r
plot(dat, xlab = "location [microns]", ylab = "Absolute Intensity")
```

![](AproxMicroResol_files/figure-markdown_github/unnamed-chunk-1-1.png)

One can get an idea about the resolution of the microscope by using the Gauss normal distribution with its equation

$$f(x) = k e ^ \\left( \\frac{(-x - \\mu\_0)^2}{(2 \\sigma)^2} \\right) $$

with

$$ k = \\frac{1}{\\sqrt{2 \\pi \\sigma^2}} $$

We will now fit the intensity profile to the provided data and determine the resolution approximation by calculating the Full Width at Half Maximum value (FWHM). This is the distance between the two points of the x-axis, where the function value of f is at its half maximum. By doing some algebra with the Gauss normal distribution, we can find out that this is the case at
$$2 \\sqrt{2 ln(2)} \\sigma \\approx 2.355 \\sigma$$

See this example

``` r
curve(1/(sqrt(2*pi*0.5^2))*exp(-(x)^2/(2*0.5^2)), from = -1.5, to = 1.5, xlab = "Location [microns]", ylab = "Intensity")

abline(v = 0.58875, lty = 2, lwd = 2)
abline(v = -0.58875, lty = 2, lwd = 2)
abline(h = 0.3989, lty = 3, lwd = 2)
text(0, 0.1, " = 0.5")
```

![](AproxMicroResol_files/figure-markdown_github/unnamed-chunk-2-1.png)

In this case sigma was set to 0.5. The FWHM value is 1.1775 accordingly.

Determining the Resolution
--------------------------

We use the data, that have been imported previously to the dataframe **dat**. First, we define the Gauss normal distribution as a function in R. *(Note: We can use k instead of the actual term in order to division by 0)* It is also beneficial to normalise the intensity values, because we want our fitted function to actually cross the y-axis at the lowest intensity values. A general approach to normalise data from 0 to 1 is
$$ X\_n = \\frac{(x\_i - X\_\\min)}{(X\_\\max - X\_\\min)}$$
.

``` r
gauss = function(x, x0, sigma){
  1/sqrt(2*pi*sigma^2)*exp(-(x - x0)^2 / (2 * sigma)^2)
}

curve(gauss(x, 0, 1), from = -4, to = 4)
```

![](AproxMicroResol_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
i = 1
while(i <= length(dat$Y)){
  dat$y.norm[i] = (dat$Y[i] - min(dat$Y))/(max(dat$Y) - min(dat$Y))
  i = i + 1
}
```

The loaded data is now fitted to the gauss()-function using the non-linear least square (nls) method, which is in the standard implementation of R. Non-linear fitting can easily go in the wrong direction and the process can get stuck at a local minimum of the error function. Therefore, we define some guesses for starting parameters by just looking at the plot of our data points. The parameter x0 for example has to be a value of around 0.4. We can also get a guess of the sigma and k value from the definition of the Rayleigh criteron and the gaussian distribution. (*Be aware that you have to convert nanometers to microns!*)

$$\\sigma\_{estimate}=\\frac{\\lambda\_{emission}}{2.355} = 0.140 $$
$$k\_{estimate} = \\frac{1}{\\sqrt{2 \\pi \\sigma\_{estimate}^2}} $$

``` r
sigma_estimate = 330 / 2.355
# k_estimate = (1/sqrt(2 * pi * sigma_estimate^2))

fit = nls(dat$y.norm ~ gauss(dat$X, x0, sigma),data = dat, start = list(x0 = 0.4, sigma = 0.5), algorithm = "default")

summary(fit)
```

    ## 
    ## Formula: dat$y.norm ~ gauss(dat$X, x0, sigma)
    ## 
    ## Parameters:
    ##       Estimate Std. Error t value Pr(>|t|)  
    ## x0    -0.03624    1.22795  -0.030   0.9772  
    ## sigma  0.83560    0.42763   1.954   0.0865 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4119 on 8 degrees of freedom
    ## 
    ## Number of iterations to convergence: 7 
    ## Achieved convergence tolerance: 4.103e-06

``` r
plot(dat$X, predict(fit), type = "l")
points(dat$X, dat$y.norm, pty = 2)

FWHM = summary(fit)[[10]][2] * 2 * sqrt(2 * log(2))

abline(v = (summary(fit)[[10]][1] + FWHM/2))
abline(v = (summary(fit)[[10]][1] - FWHM/2))
abline(h = 0.5)
```

![](AproxMicroResol_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
plot(fitted(fit), residuals(fit), main = "Residuals")
```

![](AproxMicroResol_files/figure-markdown_github/unnamed-chunk-4-2.png)

**NOTE: WHAT THE HECK IS WRONG HERE?**
