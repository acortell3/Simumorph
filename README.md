
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Simumorph

<!-- badges: start -->
<!-- badges: end -->

The goal of Simumorph is to simulate discrete morphometric evolution of
shapes, using outline-based geometric morphometric and harmonic
mathematical properties. For a more detailed explanation, see vignette
[here](https://acortell3.github.io/Simumorph/simumorph_vignette_1.html).
More functions, models and info will be added, but currently it supports
the following:

-   `amplitude`: It computes the amplitude of $n$ given harmonics given
    a series of Fourier parameters $an$, $bn$, $cn$ and $dn$.
-   `phase`: It computes the phase of $n$ given harmonics given a series
    of Fourier parameters $an$, $bn$, $cn$ and $dn$.
-   `morphospace`: It prepares the parameters of a morphospace
    (e.g. observed shapes) to use for self-correlated simulation.
-   `interpo_s`: Interpolates one or several “mean” shapes, given a set
    of Fourier parameters.
-   `proc_dist`: It computes the Procrustes distances of two shapes
    after GPA. It can be used with $x,y$ coordinates.
-   `build_s`: It rebuilds a shape given a set of amplitude and phase
    values.
-   `simumorph`: It combines the functions above to perform a full
    step-wise simulation.

## Installation

You can install the version version of Simumorph from
[GitHub](https://github.com/acortell3/Simumorph) with:

``` r
 install.packages("devtools")
devtools::install_github("acortell3/Simumorph")
```

# 

## Example

# 

\#This is a basic example which shows you how to solve a common problem:
\# \#`{r example} #library(Simumorph) ## basic example code #` \# \#What
is special about using `README.Rmd` instead of just `README.md`? You can
include R chunks like so: \# \#`{r cars} #summary(cars) #` \# \#You’ll
still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. \# \#You can
also embed plots, for example: \#
\#`{r pressure, echo = FALSE} #plot(pressure) #` \# \#In that case,
don’t forget to commit and push the resulting figure files, so they
display on GitHub and CRAN.
